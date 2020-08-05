import os
import subprocess
import tempfile
import numpy as np
import pandas as pd
import biotite.structure.io.pdb as PDB
import biotite.structure.io.xtc as XTC
import biotite.structure as bstr

def plot_cylinder(pt1,pt2,radius):
    out = """from pymol import cmd
    from pymol.cgo import *
    cmd.load_cgo([ALPHA,0.25,CYLINDER,{0},{1},{2},{3},{4},{5},1,0,0,1,0,0],"Cylinder")
    """.format(pt1[0],pt1[1],pt1[2],pt2[0],pt2[1],pt2[2])
    with open('test.py',"w") as outfile:
        outfile.write(out)

def temp_pdb(structure,tmp_dir):
    fd,temp_path = tempfile.mkstemp(dir=tmp_dir.name,text=True)
    pdb_out = PDB.PDBFile()
    PDB.set_structure(pdb_out,structure)
    pdb_out.write(temp_path)
    os.close(fd)
    return temp_path

def in_cylinder(pt1,pt2,radius,q):
    dif_vec = pt1 - pt2 
    dif_vec2 = q - np.transpose(np.repeat(np.reshape(pt1,(3,1)),q.shape[0],axis=1))
    rho = np.square(dif_vec2[:,0]) + np.square(dif_vec2[:,1])
    bnd1 = np.linalg.norm(dif_vec) > np.linalg.norm(dif_vec2,axis = 1)
    bnd2 = radius**2 > rho
    return np.logical_and(bnd1,bnd2)    

def calculate_densities(struct,traj,chain,radius = [6,8,9,10,11,12,13,14,15,16]):
    inds_chain = struct.chain_id == chain    
    mask_active_site = np.logical_and(inds_chain,np.in1d(struct.res_id,np.arange(3919,3930)))
    mask_exit = np.logical_and(inds_chain,np.in1d(struct.res_id, np.array([3586,3661,3668,3832,3840,4024])))
    active_sites = traj[:,mask_active_site,:]
    active_sites = np.mean(active_sites,axis = 1)
    coords_exit = traj[:,mask_exit,:]
    coords_exit = np.mean(coords_exit,axis = 1)
    density = np.empty((traj.shape[0],len(radius)))
    for j,r in enumerate(radius):
        for i in range(traj.shape[0]):
            n_atoms = np.where(in_cylinder(coords_exit[i],active_sites[i],r,traj[i]))[0].shape
            volume = np.linalg.norm(coords_exit[i] - active_sites[i]) * np.pi * r**2
            density[i,j] = n_atoms / volume
    return density

radius = [6,8,9,10,11,12,13,14,15,16]

struct = PDB.PDBFile()
struct.read("kinase_dimer.pdb")
struct = struct.get_structure()[0]

traj = XTC.XTCFile()
traj.read("kinase_dimer_nopbc_cluster_fit.xtc")
traj = traj.get_coord()

density_a = calculate_densities(struct,traj,"A",radius)
density_data_a = pd.DataFrame(data=density_a,columns = [str(i) for i in radius])
density_data_a.to_csv("density_data_kinase_a.csv")

density_d = calculate_densities(struct,traj,"D",radius)
density_data_d = pd.DataFrame(data=density_d,columns = [str(i) for i in radius])
density_data_d.to_csv("density_data_kinase_d.csv")

