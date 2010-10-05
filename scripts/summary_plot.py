"""
plot a histogram of the methylation data per chromosome.
usage:
    %prog [options] dir1/ dir2/ ...

where dir1/ and dir2/ contain the .bin files from a run of methylcoder.
as many directories as you want.
"""
import matplotlib
matplotlib.rc('axes', edgecolor='#aaaaaa', linewidth=0.9)
matplotlib.rc('font', family='sans-serif')
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
import glob
import os.path as op
import numpy as np

def avg(c, t, m):
    ac = np.fromfile(c, dtype=np.uint32)
    at = np.fromfile(t, dtype=np.uint32)
    am = np.fromfile(m, dtype=np.uint8)
    mes = []
    for mt in range(1, 4):
        ctx = (am == mt) | (am == mt + 3)
        ctx &= (at + ac) > 0
        me = ac[ctx].sum() / float(ac[ctx].sum() + at[ctx].sum())
        mes.append(me)
    return mes

import collections
def get_data(ctms, labels):
    names = collections.defaultdict(list)
    for label, ctm in zip(labels, ctms):
        avgs = []
        for c, t, m in ctm:
            mes = avg(c, t, m)
            name = op.basename(c).replace('.c.bin', '')
            pidx = name.find(".") + 1
            name = name[pidx:]
            avgs.append((name, mes))
        avgs.sort()
        for name, (cgs, chgs, chhs) in avgs:
            names[name].append((cgs, chgs, chhs))
    return dict(names)

def main(ctms, labels, opts):

    names = get_data(ctms, labels)
    plot(names, labels, opts)

def add_label(rect, label):
    y = rect.get_height() + rect.get_y() + 0.01
    x = rect.get_x() + rect.get_width() * 0.15
    plt.text(x, y, label, horizontalalignment='left',
                    verticalalignment='bottom',
                    rotation='vertical'
                )


def plot(names, labels, opts):
    xlabels = sorted(names.keys())
    n = len(names[xlabels[0]])


    fig = plt.figure()

    #ax = fig.add_axes([0, 0, 1, 1])
    spacing = 0.04
    width = ((1 - spacing * (n - 1)) / (n + 1))
    xs = 2.0 * spacing + np.arange(len(xlabels))

    for i in range(n):
        cgs = np.array([names[schr][i][0] for schr in xlabels])
        chgs = np.array([names[schr][i][1] for schr in xlabels])
        chhs = np.array([names[schr][i][2] for schr in xlabels])

        p1 = plt.bar(xs, cgs, width, color='r', linewidth=0)
        p2 = plt.bar(xs, chgs, width, color='g', bottom=cgs, linewidth=0)
        p3 = plt.bar(xs, chhs, width, color='b', bottom=cgs + chgs, linewidth=0)
        xs += width + spacing
        if n > 1:
            add_label(p3[0], labels[i])

    plt.xticks(xs - 0.66 * width  - (width + spacing) * i/2., xlabels)
    plt.legend((p1[0], p2[0], p3[0]), ('CG', 'CHG', 'CHH'))
    plt.xlabel("Sequence Id")
    plt.ylabel("Methylation  c/(c + t)")
    if opts.title: plt.title(opts.title)

    print >>sys.stderr, "saving to %s" % opts.out
    plt.savefig(opts.out)


if __name__ == "__main__":
    import optparse
    p = optparse.OptionParser(__doc__)
    p.add_option("--out", dest="out", default="/var/www/t/met.png", help=\
                 "path to save image. extension (.png/.pdf) will determine file-type")
    p.add_option("--title", dest="title", default="", help="image title")
    p.add_option("--exclude", dest="exclude", default=None,
                help="exclude these chromosomes. e.g. --exclude 'm|c'")
    p.add_option("--labels", dest="labels", default=None, help="by default, "
                 "the directories specified in args are used to label the series "
                 "when plotting more than one directory. if --labels are specified "
                 "these will be used instead. must be pipe '|' delimited: e.g. : "
                 "--labels 'series a|series b|other series' the number specified"
                 " must match the number of directories sent on args")
    opts, args = p.parse_args()

    if len(args) == 0:
        sys.exit(p.print_help())

    ctms = []
    exclude = [] if opts.exclude is None else opts.exclude.split("|")
    for bin_dir in args:
        cs = glob.glob("%s/*.c.bin" % bin_dir)
        cs = [c for c in cs if not any(c.endswith("%s.c.bin" % e) for e in exclude)]
        assert cs
        ts = [c.replace('.c.bin', '.t.bin') for c in cs]
        mt = [c.replace('.c.bin', '.methyltype.bin') for c in cs]

        ctms.append(zip(cs, ts, mt))

    labels = opts.labels.split("|") if opts.labels \
                    else [x.rstrip("/") for x in args]
    main(ctms, labels, opts)
