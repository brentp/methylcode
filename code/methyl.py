"""
convenience class to handle the files and common tasks in the output files.
currently loads the data into memory--but only per-chromosome. at some point,
can use numpy.memmap or pytables ...
or use numexpr for better efficiency
example use. let's go by chromosome and calculate the proportion of sites that are methylated
at > 0.1 ratio the .bin files in a directory: out1234n'

    >>> from methyl import MethylGroup
    >>> mg = MethylGroup('/home/brentp/work/out1234n')
    >>> mg.seqids
    ['1', '2', '3', '4', '5', 'c', 'm']

    >>> mg.prefix, mg.dir, mg.pattern
    ('thaliana_v9', '/home/brentp/work/out1234n/', 'thaliana_v9.%s.%s.bin')

    # iteritems return seqid, object.
    >>> for seqid, meth in mg.iteritems():
    ...     cg_cs, cg_ts, cg_mask = meth.as_context('CG')
    ...     methylation = cg_cs.astype('f') / (cg_ts + cg_cs)
    ...     n_methylated = (methylation > 0.1).sum()
    ...     possible_methylated = cg_mask.sum()
    ...     proportion_methylated = float(n_methylated) / possible_methylated
    ...     seqid, proportion_methylated, cg_mask.sum()
    ('1', 0.23593716391585529, 1394740)
    ('2', 0.29339426363501264, 915144)
    ('3', 0.26254715749216057, 1118062)
    ('4', 0.28750867295289878, 879170)
    ('5', 0.25998534029087783, 1260598)
    ('c', 0.0014011640439749945, 9278)
    ('m', 0.015039789734978463, 27394)


this means that 23 to 30% of C's that can be methylated in the CG context are.
"""
import numpy as np
import os.path as op
import glob
from pyfasta import Fasta

class MethylGroup(object):
    __slots__ = ('dir', # directory containing the .bin files.
              'prefix', # e.g.: thaliana_v9 (then files are thaliana_v9.1.total.bin)
              'seqids', # list of chromosomes
              'pattern', # e.g.: thaliana_v9.%s.%s.bin (no directory)
              'fasta_path', # path to fasta file.
              'fasta'#: fasta object.
              )
    def __init__(self, prefix, fasta=None, in_memory=False):
        if fasta is not None:
            self.fasta = Fasta(fasta)
            self.fasta_path = fasta
        else:
            self.fasta = self.fasta_path = None
        self._setup_paths(prefix)


    def _setup_paths(self, prefix):
        if op.isdir(prefix): # prefix is a dir
            self.dir = prefix.rstrip("/") + "/"
            files = glob.glob("%s*methyltype.bin" % self.dir)
            prefix = op.basename(files[0].replace(".methyltype.bin", ""))
            self.prefix = ".".join(prefix.split(".")[:-1])

        else: # prefix is a dir/start_of_name
            prefix = prefix.rstrip(".")
            self.dir = op.dirname(prefix) + "/"
            self.prefix = op.basename(prefix)
            files = glob.glob(prefix + "*.methyltype.bin")

        files = (op.basename(f).replace('.methyltype.bin', '') for f in files)
        self.seqids = sorted(f.replace(self.prefix + ".", "") for f in files)
        self.pattern = self.prefix + ".%s.%s.bin"

    def keys(self):
        return self.seqids

    def __getitem__(self, seqid):
        return Methyl(self, seqid, fasta=self.fasta)

    def iteritems(self):
        for seqid in self.seqids:
            yield seqid, Methyl(self, seqid, fasta=self.fasta)

class Methyl(object):
    __slots__ = ('mg', 'seqid', 'ts', 'cs', 'methyltype', 'mt', 'sequence')
    contexts = {'CG': (1, 4), 'CHG': (2, 5), 'CHH': (3, 6) }
    def __init__(self, methylgroup, seqid, fasta=None):
        self.mg = methylgroup
        self.seqid = seqid

        self.ts = np.fromfile(self.mg.dir + self.mg.pattern % \
                        (self.seqid, 't'), dtype=np.uint32)

        self.cs = np.fromfile(self.mg.dir + self.mg.pattern % \
                        (self.seqid, 'c'), dtype=np.uint32)

        self.mt = self.methyltype = np.fromfile(self.mg.dir + self.mg.pattern % \
                (self.seqid, 'methyltype'), dtype=np.uint8)

        if not fasta is None:
            self.sequence = fasta[seqid]


    def as_context(self, context):
        """context is either an integer from 1 to 6, or a string of:
        CG, CHG or CGG
        returns cs, ts, mask
        where cs and ts are masked (to zero) except for the requested context. 
        and mask is True at sites in the requested context and False in all other 
        positions.
        """
        if isinstance(context, basestring):
            context = Methyl.contexts[context.upper()]
        elif isinstance(context, (int, long)):
            context = (context, )
        for c in context:
            assert 0 < c < 7, c

        no_mask = self.mt == context[0]
        for c in context[1:]:
            no_mask |= (self.mt == c)
        cs = self.cs.copy()
        ts = self.ts.copy()

        cs[~no_mask] = 0
        ts[~no_mask] = 0
        return cs, ts, no_mask
        

    def __repr__(self):
        return "Methyl(%s, %ibp)" % (self.seqid, len(self.ts))

if __name__ == "__main__":
    m = MethylGroup('out1234n/thaliana_v9', fasta="/labdata/thaliana_v9/thaliana_v9.fasta")
    mcs, mts, mmask = m['1'].as_context('CG')

    print mcs[100:150]
    print mts[100:150]
    print mmask[100:150]
    print m['1'].sequence[100:150]
