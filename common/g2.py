import numpy
from itertools import product

def gen_sparse_smpl_freqs_ph():
    freqs_f = numpy.array([-101, -11, 0, 10, 100])
    freqs_b = numpy.array([-101, -11, 0, 10, 100])
    nw = (len(freqs_f)**2) * len(freqs_b)
    smpl_freqs = numpy.zeros((nw, 3), dtype=numpy.int32)
    for idx_w, (v1, v2, w) in enumerate(product(freqs_f, freqs_f, freqs_b)):
        smpl_freqs[idx_w, 0] = v1
        smpl_freqs[idx_w, 1] = v2
        smpl_freqs[idx_w, 2] = w
    return smpl_freqs

def gen_box_smpl_freqs_ph(nwb, nwf):
    """Generate (non-negative) sampling frequencies in a box
    """
    # -nwb+1, -nwb+1, ..., nwb-2, nwb-1
    freqs_b = range(-nwb+1, nwb)
    # -nwf, -nwf+1, ..., nwf-1
    freqs_f = range(-nwf, nwf)

    nw = (len(freqs_f)**2) * len(freqs_b)
    smpl_freqs = numpy.zeros((nw, 3), dtype=numpy.int32)
    for idx_w, (w, v1, v2) in enumerate(product(freqs_b, freqs_f, freqs_f)):
        smpl_freqs[idx_w, 0] = w
        smpl_freqs[idx_w, 1] = v1
        smpl_freqs[idx_w, 2] = v2
    return smpl_freqs
