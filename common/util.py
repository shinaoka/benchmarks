from pytriqs.gf import *

def BlockGf_from_struct(mesh, struct):
    # Without block structure
    if not isinstance(struct[0], list):
        return Gf(mesh=mesh, indices=struct)

    # With block structure
    G_lst = []
    for bl, idx_lst in struct:
        G_lst.append(Gf(mesh=mesh, indices=idx_lst))
    return BlockGf(name_list=[bl[0] for bl in struct], block_list=G_lst)

def Block2Gf_from_struct(mesh, struct):
    
    # Without block structure
    if not isinstance(struct[0], list):
        return Gf(mesh=mesh, indices=struct)

    bl_lst = [bl[0] for bl in struct]

    # With block structure
    G_lst = []
    for bl1, idx1_lst in struct:
        lst = []
        for bl2, idx2_lst in struct:
            lst.append(Gf(mesh=mesh, indices=[idx1_lst, idx1_lst, idx2_lst, idx2_lst]))
        G_lst.append(lst)

    return Block2Gf(name_list1=bl_lst, name_list2=bl_lst, block_list=G_lst)

# Not Working, tail problems...
def Block2Gf_from_fourier2D(G_tau, n_iw, struct):
    
    bl1, G1 = next(iter(G_tau))

    tau_meshes = [G1.mesh[0], G1.mesh[1]]
    beta = tau_meshes[0].beta
    iw_meshes = [MeshImFreq(beta, tau_mesh.statistic, n_iw) for tau_mesh in tau_meshes]

    #FIXME struct should be member of Block2Gf
    G_iw = Block2Gf_from_struct(mesh=MeshProduct(*iw_meshes), struct=struct)
    temp = Block2Gf_from_struct(mesh=MeshProduct(iw_meshes[0], tau_meshes[1]), struct=struct)

    for bl, G_tau_bl in G_tau:
        for tau in tau_meshes[1]:
            temp[bl][:,tau] << Fourier(G_tau_bl[:,tau])
        for iw in iw_meshes[0]:
            G_iw[bl][iw,:] << Fourier(temp[iw,:])

    return G_iw
    

def kronecker(iw, iwp):
    return iw == iwp
