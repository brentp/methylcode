import sys
import os.path as op
sys.path.insert(0, op.dirname(__file__))
import os
from __init__ import is_up_to_date_b, CPU_COUNT, bin_paths_from_fasta, \
        get_counts, _update_conversions, write_files, write_sam_commands, \
        is_same_cmd
from subprocess import Popen
from pyfasta import Fasta

# /opt/src/gmap/gmap-2010-03-09/util/gmap_setup -D ../ -d hg19e_gmap  ../hg19e.f
def gmap_setup(gsnap_dir, out_dir, ref_fasta):
    ref_base = op.splitext(ref_fasta)[0]
    ref_dir = op.dirname(ref_fasta)
    ref_name = op.basename(ref_base)
    make_path = "%(out_dir)s/Makefile.%(ref_name)s" % locals()
    gsnap_setup = op.join(gsnap_dir, "gmap_setup")
    # have to cd to the out_dir because gsnap writes to cwd.
    cmd = "\n cd %(ref_dir)s && \n"
    cmd += gsnap_setup
    cmd += " -D %(ref_dir)s -o %(make_path)s -d %(ref_base)s %(ref_fasta)s > %(out_dir)s/gmap_setup.log && "
    cmd += "\n make -f %(make_path)s coords > gmap_coords.log && "
    cmd += "\n make -f %(make_path)s gmapdb > gmap_gmapdb.log &&"
    cmd += "\n %(gsnap_dir)s/cmetindex -d %(ref_name)s -F %(ref_dir)s > gmap_cmetindex.log 2> gmap_cmetindex.error.log"
    cmd %= locals()
    print >>sys.stderr, "[ command ] $", cmd
    cmd_last = op.join(out_dir, "ran_gsnap_setup.sh")
    if not op.exists(cmd_last) or not is_up_to_date_b(ref_fasta, cmd_last) or not is_same_cmd(cmd, cmd_last):
        fh = open(cmd_last, "w")
        print >>fh, cmd
        fh.close()
    elif is_up_to_date_b(ref_fasta, cmd_last) and not is_same_cmd(cmd, cmd_last):
        fh = open(cmd_last, "w")
        print >>fh, cmd
        fh.close()
    if not op.exists(op.join(out_dir, make_path)) or not (is_up_to_date_b(op.join(out_dir, make_path), cmd_last)):
        # TODO ^ add check for cmetindex files as well.
        p = Popen(cmd.replace('\n', ' '), shell=True)
        print >>sys.stderr, "^ executing gmap/gsnap setup^"
        if p.wait() != 0:
            pclean = Popen("cd %(out_dir)s; make -f %(make_path)s clean" % locals())
            pclean.wait()
            pass
    else:
        print >>sys.stderr, "^ NOT executing gmap/gsnap setup. everything is up to date.^"
    return ref_base

def run_gsnap(gsnap_dir, gsnap_args, out_dir, ref_fasta, reads_paths, cpu_count):
    #/opt/src/gmap/gmap-2010-03-09/src/gsnap --npaths 1 --quiet-if-excessive -A sam --nofails --nthreads 4 -D ./ -d hg19e_gmap --cmet bs_reads.fasta > bs.align.sam
    ref_base = op.splitext(ref_fasta)[0]
    ref_name = op.basename(ref_base)
    ref_dir = op.dirname(ref_fasta)
    reads_paths = [op.abspath(r) for r in reads_paths]
    log = op.join(out_dir, "gsnap_run.log")

    reads_paths_str = " ".join(reads_paths)
    out_sam = op.abspath(op.join(out_dir, "methylcoded.gsnap.sam"))
    cmd = "%(gsnap_dir)s/gsnap --quiet-if-excessive -A sam"
    cmd += " --nofails --nthreads %(cpu_count)i -D %(ref_dir)s %(gsnap_args)s"
    cmd += " -d %(ref_name)s --use-cmet %(reads_paths_str)s > %(out_sam)s 2> %(log)s"
    cmd %= locals()
    cmd_path = op.join(out_dir, "ran_gsnap.sh")
    new_cmd = False
    if not is_same_cmd(cmd, cmd_path):
        fh = open(cmd_path, "w")
        print >>fh, cmd
        fh.close()
        new_cmd = True

    print >>sys.stderr, "\n" + cmd
    if not new_cmd and all(is_up_to_date_b(r, out_sam) for r in reads_paths) \
        and all(is_up_to_date_b(r, cmd_path) for r in reads_paths):
        print >>sys.stderr, "^ NOT executing gsnap. everything is up to date.^"
    else:
        fh = open(cmd_path, "w")
        print >>fh, cmd
        fh.close()
        print >>sys.stderr, "^ executing gsnap. ^"
        p = Popen(cmd.replace('\n', ' '), shell=True)
        if p.wait() != 0:
            os.unlink(out_sam)
            print >>sys.stderr, "ERROR:\n", open(log).read()
            sys.exit(1)
    return out_sam

def fastx_to_gsnap_fasta(reads, out_fa=sys.stdout):

    fh1 = open(reads[0])
    fh2 = open(reads[1]) if len(reads) > 1 else None
    read1 = fh1.readline
    read2 = fh2.readline if fh2 is not None else None
    while True:
        header1 = read1().strip()
        if header1 == "": break
        assert header1[0] in "@>"
        seq1 = read1()
        # for fastx, discard the 2 "extra" rows
        if header1[0] == "@":
            read1(); read1()
        print >>out_fa, ">" + header1[1:]
        print >>out_fa, seq1,

        if fh2 is not None:
            header2 = read2()
            assert header2.startswith(header1[:-4])
            seq2 = read2()
            print >>out_fa, seq2,
            if header1[0] == "@":
                read2(); read2()
    return out_fa.name

def parse_gsnap_sam(gsnap_f, ref_path, out_dir, paired_end, write_bin):
    fa = Fasta(ref_path)

    fc, ft, fmethyltype = \
            bin_paths_from_fasta(fa.fasta_name, out_dir)
    counts = get_counts(fc, ft, fa)
    #chr_lengths = dict((k, len(fa[k])) for k in fa.iterkeys())


    print >>sys.stderr, "tabulating methylation"

    for sline in open(gsnap_f):
        if sline.startswith("@SQ"):
            continue

        # the ends didn't map to same spot.
        line = sline.split("\t")
        sam_flag = int(line[1])
        if paired_end:
            if line[6] != "=": continue
        else:
            # no reported alignments.
            if sam_flag == 4: continue

        seqid = line[2]
        aln_seq = line[9]
        read_length = len(aln_seq)
        bp0 = int(line[3]) - 1
        ga = ((sam_flag & 16) != 0) ^ (sam_flag & 128 != 0)
        insert_length = int(line[8])
            #line[9] = aln_seq
            #line[10] = line[10][:len(aln_seq)]

        # both ends start at exactly the same place.
        if insert_length == 0: continue
        # handle overlapping reads. one side has + insert, the other is -
        if -read_length < insert_length < 0:
            insert_length = abs(insert_length)
            aln_seq = aln_seq[:-(read_length - insert_length)]
            read_length = len(aln_seq)
        if line[7] == '0': continue

        bp1 = bp0 + read_length
        ref_seq = (fa[seqid][bp0:bp1]).upper()


        letters = 'GA' if ga else 'CT'
        read_length = len(ref_seq)
        assert read_length > 0, (bp0, bp1)
        _update_conversions(ref_seq, aln_seq, bp0, letters,
                            counts[seqid]['c'], counts[seqid]['t'],
                            50, read_length, line[5])

    write_files(fa.fasta_name, out_dir, counts, write_bin)

    cmd = open(out_dir +"/cmd.ran", "w")
    import datetime
    print >>cmd, "#date:", str(datetime.date.today())
    print >>cmd, "#path:", op.abspath(".")
    print >>cmd, " ".join(sys.argv)
    write_sam_commands(out_dir, fa, "methylcoded.gsnap")


def is_fastq(f):
    fh = open(f)
    ifastq = (fh.readline()[0] == "@")
    fh.close()
    return ifastq

def main(out_dir, ref_fasta, reads, gsnap_path, gsnap_args, write_bin):
    fa_reads = out_dir + "/" + op.basename(reads[0]).rstrip("_1") + ".fasta"
    if all(is_fastq(r) for r in reads):
        print >>sys.stderr, "using existing reads files"
        gsnap_reads = reads
    elif all(is_up_to_date_b(r, fa_reads) for r in reads):
        # fasta reads. up to date.
        gsnap_reads = reads
    else:
        # its a pair of fasta files. put into gsnap format.
        if len(reads) > 1:
            out_fa = open(fa_reads, "w")
            fastx_to_gsnap_fasta(reads, out_fa)
            out_fa.close()
            gsnap_reads = (fa_reads,)
        else:
            # if it's a single fasta, don't need to write any new files.
            gsnap_reads = (reads[0],)
    gmap_setup(gsnap_path, out_dir, ref_fasta)
    gsnap_sam = run_gsnap(gsnap_path, gsnap_args, out_dir, ref_fasta, gsnap_reads, cpu_count=CPU_COUNT)
    paired_end = len(reads) > 1

    parse_gsnap_sam(gsnap_sam, ref_fasta, out_dir, paired_end, write_bin)

if __name__ == "__main__":
    import optparse
    # NOTE: the actual entry point to this script is from __init__.py
    # which calls the main() function.
    p = optparse.OptionParser( """
convert fastq or fasta files, (paired-end or single) to fasta format used by gsnap.

      paired end fastq usage:
            %prog pair_1.fastq pair_2.fastq > some.fasta

      paired end fasta usage:
            %prog --fasta pair_1.fasta pair_2.fasta > gsnap.fasta

      single end fastq usage:
            %prog reads.fastq > gsnap.reads.fasta""")
    opts, reads = p.parse_args()
    if not len(reads) in (1, 2):
        sys.exit(p.print_help())
    fastx_to_gsnap_fasta(reads)
