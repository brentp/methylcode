"""
write the forward and reverse mappings for use by bowtie.
much of the code in this file benefited from a read of the 
code contained in bsseeker.
"""

__version__ = "0.2.1"

from pyfasta import Fasta
import numpy as np
import sys
import os.path as op
import os
from subprocess import Popen
from calculate_methylation_points import calc_methylation
from cbowtie import _update_conversions
import string
import glob
import datetime


def revcomp(s, _comp=string.maketrans('ATCG', 'TAGC')):
    return s.translate(_comp)[::-1]

CPU_COUNT = 4
try:
    import multiprocessing
    CPU_COUNT = multiprocessing.cpu_count()
except:
    import processing
    CPU_COUNT = processing.cpuCount()

def write_c2t(fasta_name):
    """
    given a fasta file and write 2 new files. given some.fasta:
        `some.forward.c2t.fasta` contains the same headers but all C's 
                                 converted to T
        `some.reverse.c2t.fasta` contains the reverse-complemented sequece
                                 with all C's converted to T.
    """
    d = op.join(op.dirname(fasta_name), "bowtie_index")
    if not op.exists(d): os.mkdir(d)

    p, ext = op.splitext(op.basename(fasta_name)) # some.fasta -> some, fasta
    revname = "%s/%s.reverse.c2t%s" % (d, p, ext)
    forname = "%s/%s.forward.c2t%s" % (d, p, ext)
    if op.exists(revname) and op.exists(forname): return forname, revname
    fasta = Fasta(fasta_name)

    reverse_fh = open(revname, 'w')
    forward_fh = open(forname, 'w')
    print >>sys.stderr, "writing: %s, %s" % (revname, forname)

    try:
        for header in fasta.iterkeys():
            seq = str(fasta[header]).upper()
            assert not ">" in seq
            print >>reverse_fh, ">%s" % header
            print >>forward_fh, ">%s" % header

            print >>reverse_fh, revcomp(seq).replace('C', 'T')
            print >>forward_fh, seq.replace('C', 'T')

        reverse_fh.close(); forward_fh.close()
    except:
        try: reverse_fh.close(); forward_fh.close()
        except: pass
        os.unlink(revname)
        os.unlink(forname)
        raise

    return forname, revname

def is_up_to_date_b(a, b):
    return op.exists(b) and os.stat(b).st_mtime >= os.stat(a).st_mtime


def run_bowtie_builder(bowtie_path, fasta_path):
    d = os.path.dirname(fasta_path)
    p, ext = op.splitext(op.basename(fasta_path)) # some.fasta -> some, fasta
    cmd = '%s/bowtie-build -f %s %s/%s > %s/%s.log' % \
                (bowtie_path, fasta_path, d, p, d, p)
    cmd = '%s/bowtie-build -f %s %s/%s > %s/%s.log' % \
                (bowtie_path, fasta_path, d, p, d, p)

    if is_up_to_date_b(fasta_path, "%s/%s.1.ebwt" % (d, p)):
        return None
    print >>sys.stderr, "running: %s" % cmd
    process = Popen(cmd, shell=True)
    return process


def run_bowtie(bowtie_path, ref_path, out_dir, reads_c2t, mismatches, threads=CPU_COUNT, **kwargs):
    sam_out_file = out_dir + "/" + op.basename(ref_path) + ".sam"
    cmd = ("%(bowtie_path)s/bowtie --sam --sam-nohead " + \
           " -v %(mismatches)d --norc -k1 --best " + \
          "--mm -p %(threads)d %(ref_path)s  -r %(reads_c2t)s " + \
          "%(sam_out_file)s") % locals() 
    # TODO: can add a sort here by the first column so that later access to the
    # raw reads file will have better caching behavior. but pretty fast for now
    print >>sys.stderr, cmd.replace("//", "/")

    if is_up_to_date_b(ref_path + ".1.ebwt", sam_out_file) and \
       is_up_to_date_b(reads_c2t, sam_out_file):
        print >>sys.stderr, "^ up to date, not running ^"
        return sam_out_file

    process = Popen(cmd, shell=True)
    process.wait()
    return sam_out_file


# <QNAME> <FLAG> <RNAME> <POS> <MAPQ> <CIGAR> <MRNM> <MPOS> \
#        <ISIZE> <SEQ> <QUAL> \ [<TAG>:<VTYPE>:<VALUE> [...]]
"""
['0', '4', '*', '0', '0', '*', '*', '0', '0', 'NATTATTTTTTTAATTGTAAATTTTTNATTGAAATTTTGTAAATA', 'IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII', 'XM:i:0\n']
"""
# http://bowtie-bio.sourceforge.net/manual.shtml#sam-bowtie-output
def parse_sam(sam_aln_file):
    #print >>sys.stderr, "parsing %s" % sam_aln_file
    for sline in open(sam_aln_file):
        if sline[0] == "@": continue
        line = sline.split("\t")
        # no reported alignments.
        if line[1] == '4': continue 
        assert line[1] == '0', line
        #print >>sys.stderr, line
        seq = line[9]
        # NM:i:2
        NM = [x for x in line[11:] if x[0] == 'N' and x[1] == 'M'][0].rstrip()
        nmiss = int(NM[-1])
        yield dict(
            read_id=int(line[0]),
            seqid=line[2],
            pos0=int(line[3]) - 1,
            mapq=line[4],
            nmiss=nmiss,
            read_sequence=seq
        )


def bin_paths_from_fasta(original_fasta, out_dir='', pattern_only=False):
    """
    given the fasta, return the paths to the binary
    files that will be created
    """
    opath = op.splitext(op.basename(original_fasta))[0]
    if pattern_only:
        return ((out_dir + "/") if out_dir else "") + opath + ".%s.*.bin"

    paths = [ out_dir + "/" + opath + ".%s.converted.bin",
              out_dir + "/" + opath + ".%s.total.bin",
              out_dir + "/" + opath + ".%s.methyl.bin",
              out_dir + "/" + opath + ".%s.methyltype.bin" ]
    if out_dir == '':
        return [p.lstrip("/") for p in paths]
    return paths

def count_conversions(original_fasta, direction, sam_file, raw_reads, out_dir,
                      allowed_mismatches,
                      read_len=76, use_existing=False, write_text_file=True):
    # direction is either 'f'orward or 'r'everse. if reverse, need to subtract
    # from length of chromsome.
    assert direction in 'fr'
    print >>sys.stderr, "tabulating %s methylation for %s" % \
            (direction, original_fasta)

    fa = Fasta(original_fasta)
    fh_raw_reads = open(raw_reads, 'r')
    fh_c2t_reads = open(raw_reads + ".c2t", 'r')

    if direction == 'r':
        chr_lengths = dict((k, len(fa[k])) for k in fa.iterkeys())
 
    # get the 4 possible binary files for each chr
    fconverted, ftotal, fmethyl, fmethyltype = \
            bin_paths_from_fasta(original_fasta, out_dir)

    counts = {}
    for k in fa.iterkeys():
        # so this will be a dict of position => conv
        # the keys in 'n', 't' are the integer positions.
        #            total reads covering a particular c
        # here we add to total and converted np.fromfile() from the forward,
        # and add to it in the reverse. otherwise, just overwriting
        # below.
        if use_existing:
            counts[k] = {'total': np.fromfile(ftotal % k, dtype=np.uint32),
                         'converted': np.fromfile(fconverted % k,
                                                  dtype=np.uint32)}
        else:
            counts[k] = {'total': np.zeros((len(fa[k]),), dtype=np.uint32),
                         # total reads in which this c changed to t 
                         'converted': np.zeros((len(fa[k]),), dtype=np.uint32)}
    assert len(fa[k]) == len(counts[k]['total']) == len(counts[k]['converted'])

    skipped = 0
    align_count = 0
    pairs = "CT" if direction == "f" else "GA" # 
    for p in parse_sam(sam_file):
        # this is also the line number from the original file.
        read_id = p['read_id']
        pos0 = p['pos0']
        seqid = p['seqid']
        align_count += 1

        fh_raw_reads.seek((read_id * read_len) + read_id)
        raw = fh_raw_reads.read(read_len)

        if direction == 'r':
            # reason out that this pos is right
            pos0 = chr_lengths[seqid] - pos0 - read_len
            # since the read matched the flipped genome. we flip it here.
            raw = revcomp(raw)


        print >>sys.stderr, "f['%s'][%i:%i + %i]" % (seqid, pos0, pos0, read_len)
        genomic_ref = str(fa[seqid][pos0:pos0 + read_len])

        # so the position is the line num * the read_id + read_id (where the +
        # is to account for the newlines.
        DEBUG = True
        if DEBUG:
            #print >>sys.stderr, p
            fh_c2t_reads.seek((read_id * read_len) + read_id)
            print >>sys.stderr, p['nmiss']
            print >>sys.stderr, "ref        :", genomic_ref
            if direction == 'r':
                print >>sys.stderr, "raw_read(r):", raw
                c2t = fh_c2t_reads.read(read_len)
                assert c2t == p['read_sequence']
                c2t = revcomp(c2t)

            else:
                print >>sys.stderr, "raw_read(f):", raw
                c2t = fh_c2t_reads.read(read_len)
                assert c2t == p['read_sequence']
            print >>sys.stderr, "c2t        :",  c2t, "\n"
            raw_input("press any key...\n")
        # so we have to keep the ref in forward here to get the correct bp
        # positions. so we look for CT when forward and GA when back.
        current_mismatches = p['nmiss']
        # we send in the current mismatches and allowed mismatches so that in
        # the case where a C in the ref seq has becomes a T in the align seq
        # (but wasnt calc'ed as a mismatch because we had converted C=>T. if
        # these errors cause the number of mismatches to exceed the number of
        # allowed mismatches, we dont include the stats for this read.
        remaining_mismatches = allowed_mismatches - current_mismatches
        skipped += _update_conversions(genomic_ref, raw, pos0, pairs, 
                                       counts[seqid]['total'], 
                                       counts[seqid]['converted'],
                                      remaining_mismatches, read_len)
    print >>sys.stderr, \
            "skipped %i alignments where read T mapped to genomic C" % skipped
    print >>sys.stderr, "total alignments: %i" % align_count

    for seqid in sorted(counts.keys()):
        total     = counts[seqid]['total']
        converted = counts[seqid]['converted']

        total.tofile(ftotal % seqid)
        converted.tofile(fconverted % seqid)

        if write_text_file:
            file_pat = bin_paths_from_fasta(original_fasta, out_dir, 
                                            pattern_only=True)

            print >>sys.stderr, "#> writing:", file_pat % seqid
            np.seterr(divide="ignore")
            assert np.all(converted <= total)

            meth = (1.0 - (converted/total.astype(np.float32))).astype(np.float32)
            meth[np.isnan(meth) | np.isinf(meth)] = 0
            meth.tofile(fmethyl % seqid)
            del meth


            seq = str(fa[seqid])
            mtype = calc_methylation(seq)
            mtype.tofile(fmethyltype % seqid)
            out = open(op.dirname(fmethyltype) + "/methyl-data-%s.txt" \
                    % datetime.date.today(), 'w')
            to_text_file(total, converted, mtype, seqid, out)


def to_text_file(total, converted, methylation_type, seqid, out=sys.stdout):
    """
    convert the numpy arrays to a file of format:
    seqid [TAB] bp(0-based) [TAB] total-reads covering this bp [TAB] reads where c => t
    """
    print >>out, make_header()
    idxs, = np.where(total)
    for bp, mt, tt, cc in np.column_stack((idxs, methylation_type[idxs],
                                           total[idxs], converted[idxs])):
        print >>out, "\t".join(map(str, (seqid, mt, bp, tt, cc)))


def convert_reads_c2t(reads_path):
    """
    assumes all capitals returns the new path and the read_length.
    """
    c2t = reads_path + ".c2t"
    if is_up_to_date_b(reads_path, c2t): 

        return c2t, len(open(reads_path).readline().strip())
    print >>sys.stderr, "converting C to T in %s" % (reads_path)
    out = open(c2t, 'wb')
    for line in open(reads_path):
        out.write(line.replace('C', 'T'))
    out.close()
    return c2t, len(line) 

def make_header():
    return """
#created by: %s
#on: %s
#from: %s
#version: %s
""" % (" ".join(sys.argv), datetime.date.today(), op.abspath("."), __version__)


if __name__ == "__main__":
    import optparse
    p = optparse.OptionParser(__doc__)

    p.add_option("--bowtie", dest="bowtie", help="path to bowtie directory")
    p.add_option("--reads", dest="reads", help="path to reads file containing only sequence")
    p.add_option("--outdir", dest="out_dir", help="path to a directory in "
                 "which to write the files", default=None)

    p.add_option("--mismatches", dest="mismatches", default=2, type="int",
             help="number of mismatches allowed. sent to bowtie executable")
    p.add_option("--fasta", dest="fasta",
             help="path to fasta file to which to align reads")

    opts, args = p.parse_args()

    if not (opts.reads and opts.bowtie):
        sys.exit(p.print_help())

    fasta = opts.fasta
    if fasta is None:
        assert len(args) == 1, "must specify path to fasta file"
        fasta = args[0]
        assert os.path.exists(fasta), "fasta: %s does not exist" % fasta
    if glob.glob("%s/*.bin" % opts.out_dir):
        print >>sys.stderr, "PLEASE use an empty out directory or move "\
                "the existing .bin files from %s" % opts.out_dir
        sys.exit(1)

    fforward_c2t, freverse_c2t = write_c2t(fasta)
    pforward = run_bowtie_builder(opts.bowtie, fforward_c2t)
    preverse = run_bowtie_builder(opts.bowtie, freverse_c2t)
    if preverse is not None: preverse.wait()
    if pforward is not None: pforward.wait()

    raw_reads = opts.reads
    c2t_reads, read_len = convert_reads_c2t(raw_reads)  
    ref_forward = op.splitext(fforward_c2t)[0]
    ref_reverse = op.splitext(freverse_c2t)[0]
    forward_sam = run_bowtie(opts.bowtie, ref_forward, opts.out_dir,
                             c2t_reads, opts.mismatches)
    reverse_sam = run_bowtie(opts.bowtie, ref_reverse, opts.out_dir,
                             c2t_reads, opts.mismatches)
    try:
        count_conversions(fasta, 'f', forward_sam, raw_reads, opts.out_dir,
                          opts.mismatches,
                          read_len, use_existing=False, write_text_file=False)
        count_conversions(fasta, 'r', reverse_sam, raw_reads, opts.out_dir,
                      opts.mismatches,
                      read_len, use_existing=True, write_text_file=True)
    except:
        files = bin_paths_from_fasta(fasta, opts.out_dir, pattern_only=False)
        for f in files:
            try: os.unlink(f)
            except OSError: pass
        print >>sys.stderr, "ERROR: don't use .bin or text files"
        raise
    finally:
        cmd = open(opts.out_dir +"/cmd.ran", "w")
        print >>cmd, "#date:", str(datetime.date.today())
        print >>cmd, "#path:", op.abspath(".")
        print >>cmd, " ".join(sys.argv)

    print >>sys.stderr, "SUCCESS"
