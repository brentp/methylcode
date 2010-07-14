"""
convert bisulfite treated reads generated with Cokus protocol to format suitable for methylcoder.py

%prog [options] > newfile.fastq

then newfile.fastq is suitable to send to methylcoder
"""
import sys
from pyfasta import complement

def gen_record_from_raw(reads, na="I"):
    for i, line in enumerate(open(reads)):
        line = line.strip()
        yield [">%i" % i,
               "%s" % line]

def gen_record_from_fastq(reads):
    fh = open(reads)
    rl = fh.readline
    while True:
        header = rl().strip()
        if not header: break
        assert header[0] == "@"
        yield [header, rl().strip(), rl().strip(), rl().strip()]


def main(reads, tags, fmt):
    gen_record = gen_record_from_raw if fmt == "raw" else gen_record_from_fastq

    for record in gen_record(reads):
        # TODO: support case where there are no tags and just want to add both normal
        # and rc read to new file.
        if tags:
            seq = record[1]
            if seq.startswith(tags[0]):
                seq = seq[len(tags[0]):]
            elif seq.startswith(tags[1]):
                seq = complement(seq[len(tags[1]):])[::-1]
            else:
                print >>sys.stderr, "warning:%s does not start with specified tags!"
            record[1] = seq
            if len(record) > 2:
                record[3] = record[3][-len(seq):][::-1] # strip the stuff for tags.
        print "\n".join(record)
        """
        if not tags:
            record[0] += "r"
            record[3] = complement(record[3])[::-1]
            print "\n".join(record)
        """

if __name__ == "__main__":
    import optparse
    p = optparse.OptionParser(__doc__)
    p.add_option("--tags", dest="tags", help=\
                 "if reads have tags indicating their origin, "
                 "specify as the 'forward|reverse' tag, e.g. 'TGCCT|CCATA' ",
                    default=None)
    p.add_option("--format", dest="format", help=\
            "format of the reads. one of raw,fastq", default="raw")

    opts, args = p.parse_args()
    fmt = opts.format.strip()
    if len(args) != 1 or not fmt in ("raw", "fastq"):
        sys.exit(not p.print_help())
    tags = opts.tags and opts.tags.strip().split("|")
    reads = args[0]

    main(reads, tags, fmt)
