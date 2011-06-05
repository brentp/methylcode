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
    try:
        mem_kb, seconds = line.rstrip().split()
    except:
        print >>sys.stderr, "skipping ", f, "bad memory/time info"
        return None, None
    return int(mem_kb) / 1024., float(seconds) / 60

def get_name(f):
    # remove .time
    return op.basename(f)[:-5]

def get_count(group, directory):
    cf = glob.glob(op.join(directory, "%s*.count" % group))
    assert len(cf) == 1, (cf, group, directory)
    return int(float(open(cf[0]).read()))


def make_table(time_files, directory):
    tails = [(get_name(f), tail(f)) for f in time_files]
    tails = [t for t in tails if not None in t[1]]
    longest = sorted(tails, key=lambda a: len(a[0]))[-1][0]
    llen = len(longest)
    sep = "=" * 19
    row_format = "%22s %" + str(llen) + "s %19s %19s %19s"
    header_string = "%s %s %s %s %s" % ("=" * 22, "=" * (llen + 3), sep, sep, sep)
    print header_string
    print row_format % ("program", "process", "memory (MB)", "time (minutes)",
                                                                "pairs-mapped")
    print header_string
    info = {}
    for full_group, li in itertools.groupby(sorted(tails), lambda a: a[0].split(".")[0]):
        group = "-".join(full_group.split("-")[:2])
        li = list(li)
        if len(li) == 1: group = "**%s**" % group
        count = get_count(group, directory) if len(li) == 1 else ""
        for name, (mem, seconds) in li:
            name = "-".join(name.split("-")[-2:])
            print row_format % (group, name, "%1.f" % mem, "%.1f" % seconds,
                    str(count))

        count = get_count(group, directory)
        max_mem = max(ms[1][0] for ms in li)
        total_time = sum(ms[1][1] for ms in li)
        info[group] = {'count': count, 'mem': max_mem, 'time': total_time }

        if len(li) > 1:
            print row_format % ("**%s**"%  group, "total", "%1.f" % max_mem, "%.1f"
                    % total_time, str(count))
    print header_string
    return info

GCE = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789-.'

def encode(numbers, GCE=GCE):
    """
    do extended encoding on a list of numbers for the google chart api
    >>> encode([1690, 90,1000])
    'chd=e:aaBaPo'
    """
    encoded = []
    for number in numbers:
        if number > 4095: raise ValueError('too large')
        first, second = divmod(number, len(GCE))
        encoded.append("%s%s" % (GCE[first], GCE[second]))
    return "chd=e:%s" % ''.join(encoded)


def make_chart_urls(info):
    url= 'https://chart.googleapis.com/chart?cht=bvs&chs=200x125&'





def main():
    p = optparse.OptionParser(__doc__)
    opts, args = p.parse_args()
    if (len(args) == 0):
        sys.exit(not p.print_help())

    directory = args[0]
    time_files = glob.glob(op.join(directory, "*.time"))
    info = make_table(time_files, directory)
    make_chart_urls(info)


if __name__ == "__main__":
    import doctest
    if doctest.testmod(optionflags=doctest.ELLIPSIS |\
                                   doctest.NORMALIZE_WHITESPACE).failed == 0:
        main()
