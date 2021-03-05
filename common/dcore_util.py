#!/usr/bin/env python


from dcore.tools import make_block_gf, raise_if_mpi_imported

import sys, os, time
import numpy
from itertools import product

import triqs
from triqs.gf import *
from triqs.operators.util.extractors import extract_U_dict4

raise_if_mpi_imported()

def convert_to_dcore_format(gf_struct, h_int, G0_iw, beta, n_iw):
     """
     Convert input data to DCore format:
         * gf_struct is a dict.
         * Allowed block names are "ud", "up", "down".
     """
     # --------- Convert gf_struct ----------
     num_blocks = len(gf_struct)
     for b in gf_struct:
         assert b[0] in ['up', 'dn', 'bl']
     bname_tr = {'up':'up', 'dn':'down', 'bl':'ud'}
     gf_struct_dcore = {bname_tr[b[0]] : b[1] for b in gf_struct}
     gf_struct_dcore_list = [(bname_tr[b[0]], b[1]) for b in gf_struct]
     norb = len(gf_struct[0][1])//2 if num_blocks == 1 else len(gf_struct[0][1])

     if num_blocks == 1:
         idx_tr = {('bl', i): i for i in range(2*norb)}
         G0_iw_dcore = make_block_gf(GfImFreq, gf_struct_dcore, beta, n_iw)
         G0_iw_dcore['ud'] = G0_iw['bl']
     else:
         idx_tr = {(b, i): i+ispin*norb for ispin, b in enumerate(['up', 'dn']) for i in range(norb)}
         G0_iw_dcore = make_block_gf(GfImFreq, gf_struct_dcore, beta, n_iw)
         G0_iw_dcore['up'] = G0_iw['up']
         G0_iw_dcore['down'] = G0_iw['dn']

     # --------- Construct Coulomb tensor from h_int ----------
     U_dict = extract_U_dict4(h_int)
     u_mat = numpy.zeros((2*norb,)*4, dtype=complex)
     for idx4, v in list(U_dict.items()):
         #print("U_dict: ", idx4, v)
         idx4_ = [idx_tr[idx] for idx in idx4]
         u_mat[idx4_[0], idx4_[1], idx4_[2], idx4_[3]] += v
         #print(idx4, idx4_, u_mat[idx4_[0], idx4_[1], idx4_[2], idx4_[3]])

     return norb, gf_struct_dcore, u_mat, G0_iw_dcore

def convert_to_triqs_bname(G, gf_struct, beta, n_iw):
    gf_struct_dict = {b[0] : b[1] for b in gf_struct}
    G_copy = make_block_gf(GfImFreq, gf_struct_dict, beta, n_iw)
    bname_trans = {'up':'up', 'down':'dn', 'ud':'bl'}
    for bname, g in G:
        G_copy[bname_trans[bname]] = g
    return G_copy


def flatten_so_idx_G2loc(G2loc, nflavors):
    """Flatten the spin-orbital indices of G2loc
    Return: ndarray
      (nflavors, nflavors, nflavors, nflavors, nfreqs)
    """

    if isinstance(G2loc, numpy.ndarray):
        return G2loc

    nfreqs = next(iter(G2loc.values())).shape[-1]
    G2loc_flatten = numpy.zeros((nflavors,)*4 + (nfreqs,), dtype=numpy.complex128)
    for k, v in G2loc.items():
        G2loc_flatten[k[0], k[1], k[2], k[3]] = v
    return G2loc_flatten


def make_triqs_g2(beta, nwb, nwf, data, gf_struct):
    """ Convert an array of G2 to TRIQS gf object

    data: ndarray of shape (nf, nf, nf, nf, nfreqs)
        nf is the number of flavors, nfreqs is number of frequencies.
        nfreqs must be (2*nwb-1)*(2*nwf)*(2*nwf)

    gf_struct: 

    return: TRIQS two-particle gf object
      The order of blocks is AABB.
    """

    nblocks = len(gf_struct)
    nf = data.shape[0]
    norb = nf//2
    nfreqs = data.shape[-1]
    block_size = nf//nblocks
    block_shape = (block_size,)*4
    bnames = [x[0] for x in gf_struct]
    assert bnames == ['up', 'dn'] or bnames == ['bl']

    # Change the order of spin and orbital
    # For each flavor index,
    #  from  (up, dn, ..., up, dn) to (up, up, ..., dn, dn, ...)
    data = data.reshape((norb,2)*4 + (nfreqs,))
    data = data.transpose((1, 0, 3, 2, 5, 4, 7, 6, 8))
    data = data.reshape((2*norb,)*4 + (nfreqs,))

    iwmesh = MeshProduct(
                MeshImFreq(beta=beta, S='Boson', n_max=nwb),
                MeshImFreq(beta=beta, S='Fermion', n_max=nwf),
                MeshImFreq(beta=beta, S='Fermion', n_max=nwf))
    freq_shape = (2*nwb-1, 2*nwf, 2*nwf)

    G4iw_blocks = []
    offset1 = 0
    for name1, basis1 in gf_struct:
        offset2 = 0
        subblocks = []
        for name2, basis2 in gf_struct:
            G4iw_block = Gf(mesh=iwmesh, target_shape=block_shape)
            s1, e1 = offset1, offset1 + block_size
            s2, e2 = offset2, offset2 + block_size
            G4iw_block.data[...] = \
                    data[s1:e1, s1:e1, s2:e2, s2:e2, :].reshape(block_shape + freq_shape).transpose((4, 5, 6, 0, 1, 2, 3))
            offset2 += len(basis2)
            subblocks.append(G4iw_block)
        G4iw_blocks.append(subblocks)
        offset1 += len(basis1)

    return Block2Gf(bnames, bnames, G4iw_blocks)
