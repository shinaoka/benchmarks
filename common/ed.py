import sys, os
sys.path.append(os.getcwd() + "/..")
sys.path.append(os.getcwd() + "/../../common")
from model import *

from pytriqs.archive import HDFArchive
from pytriqs.utility import mpi
from pyed.TriqsExactDiagonalization import TriqsExactDiagonalization

# --------- Construct the ED solver ----------
ed = TriqsExactDiagonalization(h_imp, fundamental_operators, beta)

# --------- Calculate the single-particle Green function ----------
G_iw = G0_iw.copy() 
for bl, idx_lst in gf_struct:
    for i, j in product(idx_lst, idx_lst):
        ed.set_g2_iwn(G_iw[bl][i,j], c(bl,i), c_dag(bl,j))

# --------- Calculate chi3ph ----------
n_tau = 10 * (2 * n_iw + 1)
tau_mesh = MeshImTime(beta, 'Fermion', n_tau)
chi3pp_tau = Block2Gf_from_struct(mesh=MeshProduct(tau_mesh, tau_mesh), struct=gf_struct)
chi3ph_tau = chi3pp_tau.copy()

for (bl1, idx1_lst), (bl2, idx2_lst) in product(gf_struct, gf_struct):
    # FIXME chi3ph_tau[(bl1, bl2)][i,j,k,l] not assignable
    # for i, j, k, l in product(idx_lst, idx_lst, idx_lst, idx_lst): 
        # ed.set_g3_tau(  chi3ph_tau[(bl1, bl2)][i,j,k,l], \
                        # c_dag(bl1,i), c(bl1,j), c_dag(bl2,k) * c(bl2,l) )
        ed.set_g3_tau(  chi3pp_tau[(bl1, bl2)], \
                        c_dag(bl1,0), c_dag(bl2,0), c(bl1,0) * c(bl2,0) )
        chi3pp_tau[(bl1,bl2)] *= -1.0  # Resort operators, since we want <cdag c(0^{+}) cdag c(0)>
        ed.set_g3_tau(  chi3ph_tau[(bl1, bl2)], \
                        c_dag(bl1,0), c(bl1,0), c_dag(bl2,0) * c(bl2,0) )

#chi3pp_iw = Block2Gf_from_fourier2D(chi3pp_tau, n_iw, gf_struct)
#chi3ph_iw = Block2Gf_from_fourier2D(chi3ph_tau, n_iw, gf_struct)

# -------- Save in archive ---------
with HDFArchive("../results/ed.h5",'w') as res:
    res["G"] = G_iw
    res["chi3pp_tau"] = chi3pp_tau
    res["chi3ph_tau"] = chi3ph_tau
