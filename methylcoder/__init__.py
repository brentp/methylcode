"""
generate methylation data given a reference genome and a set of bisulfite
treated reads.
"""

__version__ = "0.2.5"

from pyfasta import Fasta
import numpy as np
import bsddb
import sys
import os.path as op
import os
from subprocess import Popen
from calculate_methylation_points import calc_methylation
from cbowtie import _update_conversions
from fastqindex import FastQIndex, FastQEntry
import string
import glob
import datetime
np.seterr(divide="ignore")


def revcomp(s, _comp=string.maketrans('ATCG', 'TAGC')):
    return s.translate(_comp)[::-1]

CPU_COUNT = 4
try:
    import multiprocessing
    CPU_COUNT = multiprocessing.cpu_count()
except ImportError:
    import processing
    CPU_COUNT = processing.cpuCount()

def write_c2t(fasta_name):
    """
    given a fasta file, write a new file:
        `some.fr.c2t.fasta` which contains:
          + the same headers prefixed with 'f' with all C's converted to T
          + headers prefixed with 'r' reverse complemented with
                                 all C's converted to T.
    """
    d = op.join(op.dirname(fasta_name), "bowtie_index")
    if not op.exists(d): os.mkdir(d)

    p, ext = op.splitext(op.basename(fasta_name)) # some.fasta -> some, fasta
    fname = "%s/%s.fr.c2t%s" % (d, p, ext)
    if op.exists(fname): return fname
    fasta = Fasta(fasta_name)

    c2t_fh = open(fname, 'w')
    print >>sys.stderr, "writing forward and reverse c2t to: %s" % (fname,)

    try:
        for header in fasta.iterkeys():
            seq = str(fasta[header]).upper()
            assert not ">" in seq
            # c2t, prefix header with f and write
            print >>c2t_fh, ">f%s" % header
            print >>c2t_fh, seq.replace('C', 'T')
            # then r-c, c2t, prefix header with r and write
            print >>c2t_fh, ">r%s" % header
            print >>c2t_fh, revcomp(seq).replace('C', 'T')

        c2t_fh.close()
    except:
        try: c2t_fh.close()
        except: pass
        os.unlink(fname)
        raise

    return fname

def is_up_to_date_b(a, b):
    return op.exists(b) and os.stat(b).st_mtime >= os.stat(a).st_mtime


def run_bowtie_builder(bowtie_path, fasta_path):
    d = os.path.dirname(fasta_path)
    p, ext = op.splitext(op.basename(fasta_path)) # some.fasta -> some, fasta
    cmd = '%s/bowtie-build -q -f %s %s/%s > %s/%s.bowtie-build.log' % \
                (bowtie_path, fasta_path, d, p, d, p)

    if is_up_to_date_b(fasta_path, "%s/%s.1.ebwt" % (d, p)):
        return None
    print >>sys.stderr, "running: %s" % cmd
    process = Popen(cmd, shell=True)
    process.wait()

def is_same_cmd(cmd, prev_file):
    return op.exists(prev_file) and \
        cmd.strip() == open(prev_file).read().strip()

def run_bowtie(opts, ref_path, reads_c2t, bowtie_args, threads=CPU_COUNT):
    out_dir = opts.out_dir
    bowtie_path = opts.bowtie
    sam_out_file = out_dir + "/" + op.basename(ref_path) + ".sam"
    cmd = ("%(bowtie_path)s/bowtie %(bowtie_args)s --sam --sam-nohead " + \
           "--chunkmbs 1024 --norc " + \
           "--best -p %(threads)d %(ref_path)s -q %(reads_c2t)s") % locals()

    if opts.k is not None: cmd += " -k %i" % opts.k
    if opts.m is not None: cmd += " -m %i" % opts.m
    if opts.mismatches: cmd += " -v  %i" % opts.mismatches

    cmd += " %(sam_out_file)s 2>&1 | tee %(out_dir)s/bowtie.log" % locals()
    print >>sys.stderr, cmd.replace("//", "/")

    if is_up_to_date_b(ref_path + ".1.ebwt", sam_out_file) and \
       is_up_to_date_b(reads_c2t, sam_out_file) and \
       is_same_cmd(cmd, sam_out_file + ".bowtie.sh"):

        print >>sys.stderr, "^ up to date, not running ^"
        return sam_out_file

    process = Popen(cmd, shell=True)
    print >> open(sam_out_file + ".bowtie.sh", "w"), cmd
    process.wait()
    return sam_out_file


# http://bowtie-bio.sourceforge.net/manual.shtml#sam-bowtie-output
def parse_sam(sam_aln_file, chr_lengths, get_records):

    for sline in open(sam_aln_file):
        # comment.
        if sline[0] == "@": continue
        # it was excluded because of -M
        line = sline.split("\t")
        # no reported alignments.
        if line[1] == '4': continue
        # extra found via -M
        if line[4] == '0': continue
        assert line[1] == '0', line

        read_id = line[0]
        seqid = line[2]
        direction = seqid[0]
        assert direction in 'fr'

        seqid = seqid[1:]
        line[2] = seqid

        pos0 = int(line[3]) - 1
        converted_seq = line[9]

        # we want to include the orginal, non converted reads
        # in the output file to view the alignment.
        # read_id is the line in the file.
        #fh_raw_reads.seek((read_id * read_len) + read_id)
        #raw_seq = fh_raw_reads.read(read_len)
        raw_fastq, converted_fastq = get_records(read_id)
        read_len = len(converted_seq)
        raw_seq = raw_fastq.seq

        if direction == 'f':
            line[9] = raw_seq
        else:
            pos0 = chr_lengths[seqid] - pos0 - read_len
            line[3] = str(pos0 + 1)
            # since the read matched the flipped genome. we flip it here.
            line[9] = raw_seq = revcomp(raw_seq)
            # flip the quality as well.
            line[10] = line[10][::-1]
            line[1] = '16' # alignment on reverse strand.
            converted_seq = revcomp(converted_fastq.seq)

        # NM:i:2
        NM = [x for x in line[11:] if x[0] == 'N' and x[1] == 'M'][0].rstrip()
        nmiss = int(NM[-1])
        line[-1] = line[-1].rstrip()
        yield dict(
            read_id=read_id,
            seqid=line[2],
            pos0=pos0,
            mapq=line[4],
            nmiss=nmiss,
            read_sequence=converted_seq,
            raw_read=raw_seq,
        ), line, read_len, direction

def bin_paths_from_fasta(original_fasta, out_dir='', pattern_only=False):
    """
    given the fasta, return the paths to the binary
    files that will be created
    """
    opath = op.splitext(op.basename(original_fasta))[0]
    if pattern_only:
        return ((out_dir + "/") if out_dir else "") + opath + ".%s.*.bin"

    paths = [ out_dir + "/" + opath + ".%s.c.bin",
              out_dir + "/" + opath + ".%s.t.bin",
              out_dir + "/" + opath + ".%s.methyltype.bin" ]
    if out_dir == '':
        return [p.lstrip("/") for p in paths]
    return paths

def get_raw_and_c2t(header, fqidx, fh_raw_reads, fh_c2t_reads):
    """
    since we're sharing the same index for the reads and the c2t,
    we take the header and return each record
    """
    fpos = fqidx.get_position(header)
    fh_raw_reads.seek(fpos)
    fh_c2t_reads.seek(fpos)
    return FastQEntry(fh_raw_reads), FastQEntry(fh_c2t_reads)

def count_conversions(original_fasta, sam_file, raw_reads, out_dir,
                      allowed_mismatches):
    # direction is either 'f'orward or 'r'everse. if reverse, need to subtract
    # from length of chromsome.
    print >>sys.stderr, "tabulating methylation for %s" % (original_fasta,)

    fqidx = FastQIndex(raw_reads + ".c2t")
    fa = Fasta(original_fasta)
    fh_raw_reads = open(raw_reads, 'r')
    fh_c2t_reads = open(raw_reads + ".c2t", 'r')

    def get_records(header):
        return get_raw_and_c2t(header, fqidx, fh_raw_reads, fh_c2t_reads)


    chr_lengths = dict((k, len(fa[k])) for k in fa.iterkeys())

    out_sam = open(out_dir + "/methylcoded.sam", 'w')

    # get the 3 possible binary files for each chr
    fc, ft, fmethyltype = \
            bin_paths_from_fasta(original_fasta, out_dir)

    counts = {}
    for k in fa.iterkeys():
        # so this will be a dict of position => conv
        # here we add to fc and ft np.fromfile() from the forward,
        # and add to it in the reverse. otherwise, just overwriting
        # below.
        counts[k] = {'t': np.zeros((len(fa[k]),), dtype=np.uint32),
                     # total reads in which this c changed to t
                     'c': np.zeros((len(fa[k]),), dtype=np.uint32)}
        assert len(fa[k]) == len(counts[k]['t']) == len(counts[k]['c'])

    skipped = 0
    align_count = 0
    for p, sam_line, read_len, direction in parse_sam(sam_file, chr_lengths, get_records):
        # read_id is also the line number from the original file.
        pairs = "CT" if direction == "f" else "GA" #
        read_id = p['read_id']
        pos0 = p['pos0']
        align_count += 1
        raw = p['raw_read']

        # the position is the line num * the read_id + read_id (where the +
        # is to account for the newlines.
        genomic_ref = str(fa[p['seqid']][pos0:pos0 + read_len])
        DEBUG = False
        if DEBUG:
            araw, ac2t = get_records(read_id)
            print >>sys.stderr, "f['%s'][%i:%i + %i]" % (p['seqid'], pos0,
                                                         pos0, read_len)
            #fh_c2t_reads.seek((read_id * read_len) + read_id)
            print >>sys.stderr, "mismatches:", p['nmiss']
            print >>sys.stderr, "ref        :", genomic_ref
            if direction == 'r':
                print >>sys.stderr, "raw_read(r):", raw
                c2t = ac2t.seq
                c2t = revcomp(c2t)
                assert c2t == p['read_sequence']

            else:
                print >>sys.stderr, "raw_read(f):", raw
                c2t = ac2t.seq
                assert c2t == p['read_sequence']
            print >>sys.stderr, "c2t        :",  c2t, "\n"

        # have to keep the ref in forward here to get the correct bp
        # positions. look for CT when forward and GA when back.
        current_mismatches = p['nmiss']
        # we send in the current mismatches and allowed mismatches so that in
        # the case where a C in the ref seq has becomes a T in the align seq
        # (but wasnt calc'ed as a mismatch because we had converted C=>T. if
        # these errors cause the number of mismatches to exceed the number of
        # allowed mismatches, we dont include the stats for this read.
        remaining_mismatches = allowed_mismatches - current_mismatches
        this_skipped = _update_conversions(genomic_ref, raw, pos0, pairs,
                                       counts[p['seqid']]['c'],
                                       counts[p['seqid']]['t'],
                                      remaining_mismatches, read_len)
        if this_skipped == 0:
            # only print the line to the sam file if we use it in our calcs.
            print >>out_sam, "\t".join(sam_line)
        skipped += this_skipped
        if DEBUG:
            raw_input("press any key...\n")

    print >>sys.stderr, "total alignments: %i" % align_count
    print >>sys.stderr, \
            "skipped %i alignments (%.3f%%) where read T mapped to genomic C" % \
                  (skipped, 100.0 * skipped / align_count)

    out = open(op.dirname(fmethyltype) + "/methyl-data-%s.txt" \
                    % datetime.date.today(), 'w')
    print >>out, make_header()
    print >>out, "#seqid\tmt\tbp\tc\tt"

    f_pat = bin_paths_from_fasta(original_fasta, out_dir, pattern_only=True)
    print >>sys.stderr, "#> writing:", f_pat

    summary_counts = dict.fromkeys(('CHG', 'CHH', 'CG'))
    for ctx in summary_counts.keys():
        summary_counts[ctx] = {'cs': 0, 'ts': 0}

    for seqid in sorted(counts.keys()):
        cs = counts[seqid]['c']
        ts = counts[seqid]['t']
        csum = float(cs.sum())
        tsum = float(ts.sum())
        mask = (cs + ts) > 0

        cs.tofile(fc % seqid)
        ts.tofile(ft % seqid)


        seq = str(fa[seqid])
        mtype = calc_methylation(seq)

        # print a quick summary of methylation for each context.
        summary = "chr: %-12s cs: %-12i ts: %-12i [methylation]" \
                  % (seqid, csum, tsum)
        contexts = [('CG', (1, 4)), ('CHG', (2, 5)), ('CHH', (3, 6))]
        for context, (mp, mm) in contexts:
            ctx = (mtype == mp) | (mtype == mm)
            ctx_mask = mask & ctx
            summary_counts[context]['cs'] += cs[ctx_mask].sum()
            summary_counts[context]['ts'] += ts[ctx_mask].sum()
            meth = (cs[ctx_mask].astype('f') / (cs[ctx_mask] + ts[ctx_mask]))
            summary += (" %s: %.4f" % (context, meth.mean()))
        print >>sys.stderr, summary.rstrip(",")

        mtype.tofile(fmethyltype % seqid)
        to_text_file(cs, ts, mtype, seqid, out)

    print >>sys.stderr, "genome-wide: methylation (cs, ts)"
    for ctx in sorted(summary_counts.keys()):
        d = summary_counts[ctx]
        print >>sys.stderr, "%s: %.5f (%i %i)" % \
                (ctx, d['cs'] / float(d['cs'] + d['ts']), d['cs'], d['ts'])

def to_text_file(cs, ts, methylation_type, seqid, out=sys.stdout):
    """
    convert the numpy arrays to a file of format:
    seqid [TAB] methylation_type [TAB] bp(0) [TAB] cs [TAB] ts
    """
    idxs, = np.where(cs + ts)
    for bp, mt, c, t in np.column_stack((idxs, methylation_type[idxs],
                                           cs[idxs], ts[idxs])):
        print >>out, "\t".join(map(str, (seqid, mt, bp, c, t)))

def write_sam_commands(out_dir, fasta):
    fh_lens = open("%s/chr.lengths.txt" % out_dir, "w")
    for k in fasta.keys():
        print >>fh_lens, "%s\t%i" % (k, len(fasta[k]))
    fh_lens.close()
    out = open("%s/commands.sam.sh" % out_dir, "w")
    print >> out, """\
SAMTOOLS=/usr/local/src/samtools/samtools

$SAMTOOLS view -b -t %(odir)s/chr.lengths.txt %(odir)s/methylcoded.sam \
        -o %(odir)s/methylcoded.unsorted.bam
$SAMTOOLS sort %(odir)s/methylcoded.unsorted.bam %(odir)s/methylcoded # suffix added
$SAMTOOLS index %(odir)s/methylcoded.bam
# TO view:
# $SAMTOOLS tview %(odir)s/methylcoded.bam %(fapath)s
    """ % dict(odir=out_dir, fapath=fasta.fasta_name)
    out.close()

def convert_reads_c2t(reads_path):
    """
    assumes all capitals returns the new path and creates and index.
    """
    c2t = reads_path + ".c2t"
    idx = c2t + FastQIndex.ext

    if is_up_to_date_b(reads_path, c2t) and is_up_to_date_b(c2t, idx):
        return c2t, FastQIndex(c2t)
    print >>sys.stderr, "converting C to T in %s" % (reads_path)

    try:
        out = open(c2t, 'wb')
        db = bsddb.btopen(idx, 'n', cachesize=32768*2, pgsize=512)

        fh_fq = open(reads_path)
        tell = fh_fq.tell
        next_line = fh_fq.readline
        while True:
            pos = tell()
            header = next_line().rstrip()
            if not header: break

            db[header[1:]] = str(pos)
            seq = next_line()
            out.write(header + '\n')
            out.write(seq.replace('C', 'T'))
            out.write(next_line())
            out.write(next_line())
        out.close()
        print >>sys.stderr, "opening index"
        db.close()
    except:
        os.unlink(c2t)
        os.unlink(idx)
        raise
    return c2t, FastQIndex(c2t)

def get_fasta(opts, args):
    "all the stuff to get the fasta from cmd line in single spot"
    fasta = opts.reference
    if fasta is None:
        assert len(args) == 1, "must specify path to fasta file"
        fasta = args[0]
        assert os.path.exists(fasta), "fasta: %s does not exist" % fasta
    if glob.glob("%s/*.bin" % opts.out_dir):
        print >>sys.stderr, "PLEASE use an empty out directory or move "\
                "the existing .bin files from %s" % opts.out_dir
        sys.exit(1)
    return fasta

def get_out_dir(out_dir, reads):
    out_dir = out_dir or op.splitext(reads)[0] + "_methylcoder"
    if not op.exists(out_dir):
        os.makedirs(out_dir)
    print >>sys.stderr, "using %s for writing output" % out_dir
    return out_dir

def make_header():
    return """\
#created by: %s
#on: %s
#from: %s
#version: %s""" % (" ".join(sys.argv), datetime.date.today(),
                   op.abspath("."), __version__)

def main():
    import optparse
    p = optparse.OptionParser(__doc__)

    p.add_option("--bowtie", dest="bowtie", help="path to bowtie directory")
    p.add_option("--reads", dest="reads", help="path to fastq reads file")
    p.add_option("--outdir", dest="out_dir", help="path to a directory in "
                 "which to write the files", default=None)

    p.add_option("--bowtie_args", dest="bowtie_args", default="",
                 help="any extra arguments to pass directly to bowtie. must "
                 " be specified inside a string. e.g.: "
                 "--bowtie_args '--strata --solexa-quals'")

    p.add_option("-v", "--mismatches", dest="mismatches", type="int",
             help="number of mismatches allowed. sent to bowtie executable")
    p.add_option("--reference", dest="reference",
             help="path to reference fasta file to which to align reads")
    p.add_option("-k", dest="k", type='int', help="bowtie's -k parameter", default=None)
    p.add_option("-m", dest="m", type='int', help="bowtie's -m parameter", default=None)

    opts, args = p.parse_args()

    if not (opts.reads and opts.bowtie):
        sys.exit(p.print_help())

    out_dir = opts.out_dir = get_out_dir(opts.out_dir, opts.reads)
    fasta = get_fasta(opts, args)
    fr_c2t = write_c2t(fasta)


    run_bowtie_builder(opts.bowtie, fr_c2t)

    raw_reads = opts.reads
    c2t_reads, c2t_index = convert_reads_c2t(raw_reads)
    ref_base = op.splitext(fr_c2t)[0]

    try:

        sam = run_bowtie(opts, ref_base, c2t_reads, opts.bowtie_args)
        # if they didn't specify mismatches (used quality scores, then allow 2 more
        # from c2t coversions.
        count_conversions(fasta, sam, raw_reads, opts.out_dir, opts.mismatches or 2)
    except:
        files = bin_paths_from_fasta(fasta, opts.out_dir, pattern_only=True)
        for f in glob.glob(files):
            print >>sys.stderr, "deleting:", f
            try: os.unlink(f)
            except OSError: pass
        print >>sys.stderr, "ERROR: don't use .bin or text files"
        raise
    finally:
        cmd = open(opts.out_dir +"/cmd.ran", "w")
        print >>cmd, "#date:", str(datetime.date.today())
        print >>cmd, "#path:", op.abspath(".")
        print >>cmd, " ".join(sys.argv)
        write_sam_commands(opts.out_dir, Fasta(fasta))

    print >>sys.stderr, "SUCCESS"

if __name__ == "__main__":
    main()
