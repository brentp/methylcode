import numpy as np
import os
from run_bowtie import is_up_to_date_b
from calculate_methylation_points import calc_methylation
from pyfasta import Fasta

window = 50

fasta = Fasta("thaliana_v8.fasta")
methyl = "thaliana_v8.fasta.%i.methyl.bin"
mtype  = "thaliana_v8.fasta.%i.methyltype.bin"
converted = "thaliana_v8.fasta.%i.converted.bin"
total = "thaliana_v8.fasta.%i.total.bin"

for ichr in range(1, 6):
    i = ichr
    if not is_up_to_date_b(converted % ichr, methyl % ichr):
        tot = np.fromfile(total % i, dtype=np.uint8)
        con = np.fromfile(converted % i, dtype=np.uint8)
        meth = (1.0 - (con/tot.astype(np.float32))).astype(np.float32)
        meth[np.where(np.isnan(meth))] = 0
        meth.tofile(methyl % ichr)
    else:
        print methyl % ichr, "up to date"
        meth = np.fromfile(methyl % ichr, dtype=np.float32)
    
    if not is_up_to_date_b(total % i, mtype % i):
        seq = str(fasta[str(ichr)])
        assert seq , (fasta.keys(), fasta[str(ichr)][:10])
        print len(seq)
        print mtype % ichr, "calculating methylation type"
        mt = calc_methylation(str(seq))
        print mt.shape
        mt.tofile(total.replace('.total.bin', '.methyltype.bin') % ichr)
    else:     
        print mtype % ichr, "up to date"
        mt = np.fromfile(mtype  % ichr, dtype=np.uint8)

    assert mt.shape == meth.shape
    print mt.shape, meth.shape        
    # only want cg methy
    meth[(mt != 1) & (mt != 4)] = 0
    avg = np.convolve(meth, np.ones((window,)) / float(window),
                      mode='same').astype(np.float32)
    assert avg.shape == mt.shape
    avg.tofile("at.cg.%i.%i.bin" % (window, ichr))
    del meth, mt
