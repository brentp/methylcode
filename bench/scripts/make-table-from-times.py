"""
    %prog [options] path/to/directory-containing-time-files
"""
import optparse
import os.path as op
import sys
import glob
import itertools


def tail(f):
    for line in open(f): pass
    mem_kb, seconds = line.rstrip().split()
    return int(mem_kb), float(seconds)

def get_name(f):
    # remove .time
    return op.basename(f)[:-5]

def get_count(group, directory):
    cf = glob.glob(op.join(directory, "%s*.count" % group))
    assert len(cf) == 1, (cf, group, directory)
    return int(open(cf[0]).read())


def make_table(time_files, directory):
    tails = [(get_name(f), tail(f)) for f in time_files]
    longest = sorted(tails, key=lambda a: len(a[0]))[-1][0]
    llen = len(longest)
    sep = "=" * 19
    row_format = "%21s %" + str(llen) + "s %19s %19s %19s"
    header_string = "%s %s %s %s %s" % ("=" * 21, "=" * (llen + 3), sep, sep, sep)
    print header_string
    print row_format % ("program", "process", "memory (KB)", "time (seconds)",
                                                                "pairs-mapped")
    print header_string
    for group, li in itertools.groupby(sorted(tails), lambda a: a[0].split(".")[0]):
        group = "-".join(group.split("-")[:2])
        li = list(li)
        if len(li) == 1: group = "*%s*" % group
        count = get_count(group, directory) if len(li) == 1 else ""
        for name, (mem, seconds) in li:
            name = "-".join(name.split("-")[-2:])
            print row_format % (group, name, "%s" % mem, "%.1f" % seconds,
                    str(count))

        count = get_count(group, directory)
        if len(li) > 1:
            max_mem = max(ms[1][0] for ms in li)
            total_time = sum(ms[1][1] for ms in li)

            print row_format % ("*%s*"%  group, "total", "%s" % max_mem, "%.1f"
                    % total_time, str(count))
    print header_string




def main():
    p = optparse.OptionParser(__doc__)
    opts, args = p.parse_args()
    if (len(args) == 0):
        sys.exit(not p.print_help())

    directory = args[0]
    time_files = glob.glob(op.join(directory, "*.time"))
    make_table(time_files, directory)


if __name__ == "__main__":
    import doctest
    if doctest.testmod(optionflags=doctest.ELLIPSIS |\
                                   doctest.NORMALIZE_WHITESPACE).failed == 0:
        main()