"""
generate methylation data given a reference genome and a set of bisulfite
treated reads. usage:

    methylcoder [options] reads.fastq (or fasta)

  or for pair-end reads:

    methylcoder [options] reads1.fastq reads2.fastq

any arguments to bowtie should be sent as a string via --bowtie_args.
e.g. --bowtie_args "--phred64-quals -m 1 --fr "

NOTE: it is highly recommended to use -m 1 to ensure unique mappings.

"""
from pyfasta import Fasta
import numpy as np
import bsddb
import sys
import os.path as op
import os
from subprocess import Popen
from calculate_methylation_points import calc_methylation
from cbowtie import _update_conversions
from fastqindex import guess_index_class
import string
import glob
import datetime
np.seterr(divide="ignore")

import version
__version__ = version.version

def revcomp(s, _comp=string.maketrans('ATCG', 'TAGC')):
    return s.translate(_comp)[::-1]

CPU_COUNT = 4
try:
    import multiprocessing
    CPU_COUNT = multiprocessing.cpu_count()
except ImportError:
    import processing
    CPU_COUNT = processing.cpuCount()

def write_c2t(fasta_name, unconverted):
    """
    given a fasta file, write a new file:
        `some.fr.c2t.fasta` which contains:
          + the same headers prefixed with 'f' with all C's converted to T
          + headers prefixed with 'r' reverse complemented with
                                 all C's converted to T.

    if unconverted is false, then also save a file with the forward and reverse
    without conversion.
    """
    d = op.join(op.dirname(fasta_name), "bowtie_index")
    if not op.exists(d): os.mkdir(d)

    p, ext = op.splitext(op.basename(fasta_name)) # some.fasta -> some, fasta
    fname = "%s/%s.fr.c2t%s" % (d, p, ext)
        # no conversion, just copy the file into the index dir.
    unconverted_fname = "%s/%s.fr%s" % (d, p, ext)
    if op.exists(fname):
        if not unconverted: return fname, unconverted_fname
        elif op.exists(unconverted_fname): return fname, unconverted_fname

    fasta = Fasta(fasta_name)

    c2t_fh = open(fname, 'w')
    unc_fh = open(unconverted_fname, 'w') if unconverted else None

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
            rseq = revcomp(seq)
            print >>c2t_fh, rseq.replace('C', 'T')
            if unc_fh is not None:
                print >>unc_fh, ">f%s\n%s" % (header, seq)
                print >>unc_fh, ">r%s\n%s" % (header, rseq)

        c2t_fh.close()
    except:
        os.unlink(fname)
        os.unlink(unconverted_fname)
        raise

    return fname, unconverted_fname

def is_up_to_date_b(a, b):
    return op.exists(b) and os.stat(b).st_mtime >= os.stat(a).st_mtime


def run_bowtie_builder(bowtie_path, fasta_path):
    d = os.path.dirname(fasta_path)
    p, ext = op.splitext(op.basename(fasta_path)) # some.fasta -> some, fasta
    cmd = '%s/bowtie-build -q -f %s %s/%s > %s/%s.bowtie-build.log' % \
                (bowtie_path, fasta_path, d, p, d, p)

    if is_up_to_date_b(fasta_path, "%s/%s.1.ebwt" % (d, p)):
        print >>sys.stderr, "^ not running: %s ^" % cmd
        return None
    print >>sys.stderr, "running: %s" % cmd
    process = Popen(cmd, shell=True)
    return process

def is_same_cmd(cmd, prev_file):
    return op.exists(prev_file) and \
        cmd.strip() == open(prev_file).read().strip()

def run_bowtie(opts, ref_path, reads_list_c2t, bowtie_args, bowtie_sequence_flag,
               threads=CPU_COUNT):
    out_dir = opts.out_dir
    bowtie_path = opts.bowtie
    sam_out_file = out_dir + "/" + op.basename(ref_path) + ".sam"
    cmd = ("%(bowtie_path)s/bowtie %(bowtie_args)s --sam --sam-nohead " + \
           "--chunkmbs 1024 --norc -p %(threads)d %(ref_path)s " + \
           "-%(bowtie_sequence_flag)s ") % locals()
    # single end reads
    if len(reads_list_c2t) == 1:
        cmd +=" --best %s" % reads_list_c2t[0]
    else:
        cmd += " -1 %s -2 %s" % tuple(reads_list_c2t)

    cmd += " %(sam_out_file)s 2>&1 | tee %(out_dir)s/bowtie.log" % locals()
    print >>sys.stderr, cmd.replace("//", "/")

    if is_up_to_date_b(ref_path + ".1.ebwt", sam_out_file) and \
       is_up_to_date_b(reads_list_c2t[0], sam_out_file) and \
       is_same_cmd(cmd, sam_out_file + ".bowtie.sh"):

        print >>sys.stderr, "^ up to date, not running ^"
        return sam_out_file

    process = Popen(cmd, shell=True)
    print >> open(sam_out_file + ".bowtie.sh", "w"), cmd
    if process.wait() != 0:
        print >>sys.stderr, "errors running bowtie"
        sys.exit(1)
    return sam_out_file


# http://bowtie-bio.sourceforge.net/manual.shtml#sam-bowtie-output
def parse_sam(sam_aln_file, chr_lengths, get_records, unmapped_name):
    unmapped = open(unmapped_name, "w")
    print >>sys.stderr, "writing unmapped reads to %s" % (unmapped.name, )
    idx = 0
    for sline in open(sam_aln_file):
        # comment.
        if sline[0] == "@": continue
        # it was excluded because of -M
        line = sline.split("\t")
        read_id = line[0]
        # no reported alignments.
        # extra via -m
        if line[1] == '4':
            if not "XM:i:0" in sline:
                # write stuff that was excluded because of too many mappings.
                raw_fastq, converted_fastq = get_records(read_id, 0)
                print >> unmapped, str(raw_fastq)
            continue
        # extra found via -M
        if line[4] == '0' and line[1] == '0':
            raw_fastq, converted_fastq = get_records(read_id, 0)
            print >> unmapped, str(raw_fastq)
            continue

        pair_end = line[1] != '0'
        if pair_end:
            # if the pair doesn't map to same place, skip.
            if line[6] != "=": continue
            # flags are (1 | 2 | 32 | 64) or (1 | 2 | 16 | 128)
            idx = 0 if pair_end == '99' else 0
            # bowtie prints the alignment without the pair end info.
            # add back /0 or /1 here.
            read_id = read_id + "/" + str(idx + 1)

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
        raw_fastq, converted_fastq = get_records(read_id, idx)
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
            line[1] = str(int(line[1]) + 16) # alignment on reverse strand.
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
    IndexEntry = fqidx.entry_class
    fpos = fqidx.get_position(header)
    fh_raw_reads.seek(fpos)
    fh_c2t_reads.seek(fpos)
    return IndexEntry(fh_raw_reads), IndexEntry(fh_c2t_reads)

def methylation_summary(cs, ts, mt):
    contexts = [('CG', (1, 4)), ('CHG', (2, 5)), ('CHH', (3, 6))]
    m = {}
    for context, (mp, mm) in contexts:
        mask = (mt == mp) | (mt == mm)
        mask &= ((cs + ts) > 0)
        m[context] = {}
        m[context]['cs'] = css = cs[mask].sum()
        m[context]['ts'] = tss = ts[mask].sum()
        m[context]['methylation'] = css / float(css + tss)
    return m

def print_genome_summary(summary_counts):
    """
    like print_summary() except does the genome-wide averages"
    """
    genome_dict = dict(seqid='genome-wide', total_cs=0, total_ts=0)
    #print >>sys.stderr, "genome-wide: methylation (cs, ts)"
    for ctx in summary_counts.keys():
        d = summary_counts[ctx]
        genome_dict["%s_cs" % ctx] = d['cs']
        genome_dict["%s_ts" % ctx] = d['ts']
        genome_dict[ctx] = d['cs'] / float(d['cs'] + d['ts'])
        genome_dict['total_cs'] += d['cs']
        genome_dict['total_ts'] += d['ts']

    print >>sys.stderr, print_summary.format % genome_dict

def print_summary(seqid, cs, ts, mtype, summary_counts, print_header=False):
    ms = methylation_summary(cs, ts, mtype)
    if print_header:
        header = print_summary.format.replace("-12i", "-12s").replace(".6f", "8s")
        labels = "seqid total_cs total_ts CG CG_cs CG_ts CHG CHG_cs CHG_ts CHH CHH_cs CHH_ts".split()
        print >>sys.stderr, header % dict(zip(labels, labels))
    d = {'seqid': seqid, 'total_cs': cs.sum(), 'total_ts': ts.sum()}
    # print a quick summary of methylation for each context.
    for context in sorted(summary_counts.keys()):
        d["%s_cs" % context] = ms[context]['cs']
        d["%s_ts" % context] = ms[context]['ts']
        d[context]           = ms[context]['methylation']

        summary_counts[context]['cs'] += ms[context]['cs']
        summary_counts[context]['ts'] += ms[context]['ts']
        #summary += (" %s: %.4f" % (context, ms[context]['methylation']))

    print >>sys.stderr, print_summary.format % d
print_summary.format = "%(seqid)-12s %(total_cs)-12i %(total_ts)-12i %(CG)-.6f %(CG_cs)-12i %(CG_ts)-12i %(CHG)-.6f %(CHG_cs)-12i %(CHG_ts)-12i %(CHH)-.6f %(CHH_cs)-12i %(CHH_ts)-12i "


def get_counts(fcpatt, ftpatt, fa):

    counts = {}
    for seqid in fa.iterkeys():
        seq = fa[seqid]
        fc = fcpatt % seqid
        ft = ftpatt % seqid
        tc = {}
        if op.exists(fc) and op.exists(ft):
            print >>sys.stderr, "using existing files: %s, %s" % (fc, ft)
            tc = {'t': np.fromfile(ft, dtype=np.uint32),
                  'c': np.fromfile(fc, dtype=np.uint32)}
        else:
            tc = {'t': np.zeros((len(seq),), dtype=np.uint32),
                   # total reads in which this c changed to t
                  'c': np.zeros((len(seq),), dtype=np.uint32)}
        assert len(seq) == len(tc['t']) == len(tc['c'])
        counts[seqid] = tc
    return counts

def count_conversions(original_fasta, sam_file, raw_reads_list, c2t_reads_list, index_class, out_dir,
                      allowed_mismatches, mode='w', counts=None):
    print >>sys.stderr, "tabulating methylation for %s" % (original_fasta,)

    fh_raw_reads_list = [open(raw_reads) for raw_reads in raw_reads_list]
    fh_c2t_reads_list = [open(c2t_reads) for c2t_reads in c2t_reads_list]
    fqidx_list = [index_class(c2t_reads) for c2t_reads in c2t_reads_list]
    fa = Fasta(original_fasta)
    #fh_raw_reads = open(raw_reads, 'r')
    #fh_c2t_reads = open(c2t_reads, 'r')

    def get_records(header, idx):
        return get_raw_and_c2t(header, fqidx_list[idx], fh_raw_reads_list[idx], fh_c2t_reads_list[idx])

    chr_lengths = dict((k, len(fa[k])) for k in fa.iterkeys())

    out_sam = open(out_dir + "/methylcoded.sam", mode)

    # get the 3 possible binary files for each chr
    fc, ft, fmethyltype = \
            bin_paths_from_fasta(original_fasta, out_dir)

    # if we're re-running without conversion, it was passed via kwargs.
    # otherwise, we get a new set.
    retry = counts is not None
    unmapped_name = op.dirname(sam_file) +  "/unmapped.reads"
    if retry:
        unmapped_name += ".try2"
    else:
        counts = get_counts(fc, ft, fa)

    skipped = 0
    align_count = 0
    for p, sam_line, read_len, direction in parse_sam(sam_file, chr_lengths,
                                                      get_records, unmapped_name):
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
        # if we're retrying, the reads are treated, but the genome is not ct-converted, so
        # dont check the mismtaches since it was done correctly by bowtie.
        remaining_mismatches = -1 if retry else (allowed_mismatches - current_mismatches)

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
    return counts, unmapped_name

def write_files(original_fasta, out_dir, counts):
    fa = Fasta(original_fasta)
    fc, ft, fmethyltype = \
            bin_paths_from_fasta(original_fasta, out_dir)
    out = open(op.dirname(fmethyltype) + "/methyl-data-%s.txt" \
                    % datetime.date.today(), 'w')
    print >>out, make_header()
    print >>out, "#seqid\tmt\tbp\tc\tt"

    f_pat = bin_paths_from_fasta(original_fasta, out_dir, pattern_only=True)
    print >>sys.stderr, "#> writing:", f_pat

    summary_counts = dict.fromkeys(('CHG', 'CHH', 'CG'))
    for ctx in summary_counts.keys():
        summary_counts[ctx] = {'cs': 0, 'ts': 0}

    for i, seqid in enumerate(sorted(counts.keys())):
        cs = counts[seqid]['c']
        ts = counts[seqid]['t']

        seq = str(fa[seqid])
        mtype = calc_methylation(seq)

        print_summary(seqid, cs, ts, mtype, summary_counts, print_header=(i==0))

        cs.tofile(fc % seqid)
        ts.tofile(ft % seqid)
        mtype.tofile(fmethyltype % seqid)

        to_text_file(cs, ts, mtype, seqid, out)

    print_genome_summary(summary_counts)

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

    # the index can be either for Fasta or FastQ file.
    # determine that here by the start of the header > or @
    IndexClass = guess_index_class(reads_path)
    # if convert is false, then we just return the original file.
    c2t = reads_path + ".c2t"
    idx = c2t + IndexClass.ext

    if is_up_to_date_b(reads_path, c2t) and is_up_to_date_b(c2t, idx):
        return c2t, IndexClass(c2t)

    print >>sys.stderr, "converting C to T in %s" % (reads_path)

    try:
        out = open(c2t, 'wb')
        db = bsddb.btopen(idx, 'n', cachesize=32768*2, pgsize=512)

        fh_fq = open(reads_path)
        tell = fh_fq.tell
        next_line = fh_fq.readline
        lines = range(2, IndexClass.entry_class.lines)

        while True:
            pos = tell()
            header = next_line().rstrip('\n')
            if not header: break

            db[header[1:]] = str(pos)
            seq = next_line()
            out.write(header + '\n')
            out.write(seq.replace('C', 'T'))
            for i in lines: out.write(next_line())

        out.close()
        print >>sys.stderr, "opening index"
        db.close()
    except:
        os.unlink(c2t)
        os.unlink(idx)
        raise
    return c2t, IndexClass(c2t)

def get_fasta(opts):
    "all the stuff to get the fasta from cmd line in single spot"
    fasta = opts.reference
    if fasta is None:
        raise Exception("must specify path to fasta file")
    assert os.path.exists(fasta), "fasta: %s does not exist" % fasta
    if glob.glob("%s/*.bin" % opts.out_dir):
        cmd = raw_input("""%s already contains .bin files. do you want to:
                        (U)pate them with this data.
                        (D)elete them.
                        (A)bort ? """ % opts.out_dir).rstrip().upper()[0]
        if cmd == "A":
            sys.exit(1)
        elif cmd == "D":
            for f in glob.glob("%s/*.bin" % opts.out_dir): os.unlink(f)
    return fasta

def get_out_dir(out_dir, reads):
    out_dir = out_dir or op.splitext(reads[0])[0] + "_methylcoder"
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
    p.add_option("--outdir", dest="out_dir", help="path to a directory in "
                 "which to write the files", default=None)
    p.add_option("--unconverted", dest="unconverted", action="store_true",
                 default=False, help="map unconverted reads against "
                 "unconverted genome in addtion to mapping c2t reads against"
                 " c2t genome. only works for un-paired reads.")

    p.add_option("--bowtie_args", dest="bowtie_args", default="",
                 help="any extra arguments to pass directly to bowtie. must "
                 " be specified inside a string. e.g.: "
                 "--bowtie_args '--strata --solexa-quals'")

    p.add_option("--mismatches", dest="mismatches", type="int", default=0,
             help="number of mismatches allowed. NOT sent to bowtie executable"
                " used to improve mapping accuracy after seeing reads aligned"
                " with original genome. use bowtie_args to send to bowtie."
                " default: %default")
    p.add_option("--reference", dest="reference",
             help="path to reference fasta file to which to align reads")

    opts, read_paths = p.parse_args()

    if not (len(read_paths) in (1, 2) and opts.bowtie):
        sys.exit(p.print_help())

    unconverted = opts.unconverted and len(read_paths) == 1

    out_dir = opts.out_dir = get_out_dir(opts.out_dir, read_paths)
    fasta = get_fasta(opts)

    fr_c2t, fr_unc = write_c2t(fasta, unconverted)

    pc2t = run_bowtie_builder(opts.bowtie, fr_c2t)
    punc = run_bowtie_builder(opts.bowtie, fr_unc) if unconverted else None
    if pc2t:
        if pc2t.wait() != 0: sys.exit(1)
    if punc:
        if punc.wait() != 0: sys.exit(1)

    IndexClass = None
    c2t_reads_list = []
    for read_path in read_paths:
        c2t_reads_path, c2t_reads_index = convert_reads_c2t(read_path)
        IndexClass = c2t_reads_index.__class__
        c2t_reads_list.append(c2t_reads_path)


    ref_base = op.splitext(fr_c2t)[0]

    try:
        bowtie_reads_flag = IndexClass.entry_class.bowtie_flag
        sam = run_bowtie(opts, ref_base, c2t_reads_list, opts.bowtie_args, bowtie_reads_flag)
        counts, unmatched = count_conversions(fasta, sam, read_paths, c2t_reads_list, IndexClass, opts.out_dir, opts.mismatches)
        if unconverted and len(c2t_reads_list) == 1:
            sam = run_bowtie(opts, ref_base, unmatched, opts.bowtie_args, bowtie_reads_flag)
            counts, _ = count_conversions(fasta, sam, unmatched, unmatched,  IndexClass,
                                          opts.out_dir, opts.mismatches, mode='a', counts=counts)
        write_files(fasta, opts.out_dir, counts)
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
