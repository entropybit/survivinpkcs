
import biotite
import biotite.structure as struc
import biotite.structure.io as strucio
import biotite.structure.io.xtc as xtc
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
font = {'size'   : 16}
mpl.rc('font', **font)
import matplotlib.pyplot as plt
from biotite.structure.io import load_structure, save_structure
import os

def make_dir(directory):
    if not os.path.exists(directory):
            os.makedirs(directory)

# Put here the path of the downloaded files
templ_file_path = "output/kinase_dimer.pdb"
template_model = strucio.load_structure(templ_file_path)
#templ_file_path = "output/npt.gro"
#traj_file_path  = "output/kinase_dimer_md_center.xtc"
#traj_file_path  = "output/kinase_dimer_nowater_fit.xtc"
#traj_file_path  = "output/kinase_dimer_md.xtc"
traj_file_path = "output/kinase_dimer_nopbc_cluster_fit.xtc"

# Gromacs does not set the element symbol in its PDB files,
# but Biotite guesses the element names from the atom names,
# emitting a warning
protein = strucio.load_structure(templ_file_path)
# The structure still has water and ions, that are not needed for our
# calculations, we are only interested in the protein itself
# These are removed for the sake of computational speed using a boolean
# mask
protein_mask = struc.filter_amino_acids(protein)
template = protein[protein_mask]
# We could have loaded the trajectory also with
# 'strucio.load_structure()', but in this case we only want to load
# those coordinates that belong to the already selected atoms of the
# template structure.
# Hence, we use the 'XTCFile' class directly to load the trajectory
# This gives us the additional option that allows us to select the
# coordinates belonging to the amino acids.



print(" .. loading trajectory ...")
xtc_file = xtc.XTCFile()
#xtc_file.read(traj_file_path, 0, 10, atom_i=np.where(protein_mask)[0])
xtc_file.read(traj_file_path, atom_i=np.where(protein_mask)[0])
#xtc_file.read(traj_file_path)

trajectory = xtc_file.get_structure(template)
print(" ... done ... ")


#trajectory = struc.remove_pbc(trajectory)

print(trajectory[0].coord)

trajectory_kinase_left  = trajectory[:, trajectory.chain_id == "D"]
trajectory_kinase_right  = trajectory[:, trajectory.chain_id == "G"]

print(trajectory_kinase_left.coord.shape[0])
assert(trajectory_kinase_left.coord.shape[0] >0)

# Get simulation time for plotting purposes
time = xtc_file.get_time()
time = 10**(-3)*time

#ind_stable = np.where(time == 75)[0][0]
#print(" stable after 75ns which is frame[" +str(ind_stable) + "] ")

print("start time is ::"+ str(time[0]))
print("end time is   :: "+ str(time[-1])) 


print(" ... writing start frame ...")
frame_start = template_model.copy()
frame_start.coord = trajectory[0].coord
save_structure("frame_start_coord.pdb", frame_start)
save_structure("frame_start.pdb", trajectory[0])
print(" ... done ... ")

print(" ... writing frame[1] ... ")
frame_1 = template_model.copy()
frame_1.coord = trajectory[1].coord
save_structure("frame_1_coord.pdb", frame_1)
save_structure("frame_1.pdb", trajectory[1])
print(" ... done ... ")

print(" ... writing end frame ...")
frame_end = template_model.copy()
frame_end.coord = trajectory[-1].coord
save_structure("frame_end_coord.pdb", frame_end)
save_structure("frame_end.pdb", trajectory[-1])
print(" ... done ... ")


rmsd_overall = struc.rmsd(
    trajectory[0], trajectory
)
radius_overall = struc.gyration_radius(trajectory)



# kinase left
trajectory_kinase_left, transform = struc.superimpose(
    trajectory_kinase_left[0], trajectory_kinase_left
)
rmsd_kinase_left = struc.rmsd(
    trajectory_kinase_left[0], trajectory_kinase_left)
radius_kinase_left = struc.gyration_radius(trajectory_kinase_left)

# kinase right
trajectory_kinase_right, transform = struc.superimpose(
    trajectory_kinase_right[0], trajectory_kinase_right
)
rmsd_kinase_right = struc.rmsd(
    trajectory_kinase_right[0], trajectory_kinase_right
)
radius_kinase_right = struc.gyration_radius(trajectory_kinase_right)


figure, (ax1, ax2) = plt.subplots(2,1)

ax1.plot(time, rmsd_kinase_left, color=biotite.colors["dimorange"])
ax1.set_xlim(time[0], time[-1])
ax1.set_ylim(0, rmsd_kinase_left.max()*1.1)
ax1.set_title("Left Kinase")
ax1.set_xlabel("Time (ns)")
ax1.set_ylabel("RMSD (Å)")
ax1.set_xticks(np.arange(min(time), max(time)+1, 25))


ax12 = ax1.twinx()
ax12.plot(time, radius_kinase_left, color="gray")
ax12.set_ylim(radius_kinase_left.min()*0.8, radius_kinase_left.max()*1.1)
ax12.set_title("Left Kinase")
ax12.set_xlabel("Time (ns)")
ax12.set_ylabel("Radius (Å)")


ax2.plot(time, rmsd_kinase_right, color=biotite.colors["dimorange"])
ax2.set_xlim(time[0], time[-1])
ax2.set_ylim(0, rmsd_kinase_right.max()*1.1)
ax2.set_title("Right Kinase")
ax2.set_xlabel("Time (ns)")
ax2.set_ylabel("RMSD (Å)")
ax2.set_xticks(np.arange(min(time), max(time)+1, 25))

ax21 = ax2.twinx()
ax21.plot(time, radius_kinase_right, color="gray")
ax21.set_ylim(radius_kinase_right.min()*0.8, radius_kinase_right.max()*1.1)
ax21.set_title("Right Kinase")
ax21.set_xlabel("Time (ns)")
ax21.set_ylabel("Radius (Å)")

plt.tight_layout()
plt.savefig("rmsd_radius.png", dpi=600)
plt.clf()

#np.savetxt("rmsd_surivin_right", rmsd_survivin)
df_rmsd_kinase_right = pd.DataFrame(rmsd_kinase_right)
df_rmsd_kinase_left = pd.DataFrame(rmsd_kinase_left)
df_radius_kinase_right = pd.DataFrame(radius_kinase_right)
df_radius_kinase_left = pd.DataFrame(radius_kinase_left)

df_rmsd_kinase_right.to_csv("rmsd_kinase_right.csv")
df_rmsd_kinase_left.to_csv("rmsd_kinase_left.csv")
df_radius_kinase_right.to_csv("radius_kinase_right.csv")
df_radius_kinase_left.to_csv("radius_kinase_left.csv")

df_rmsd = pd.DataFrame(rmsd_overall)
df_radius = pd.DataFrame(radius_overall)

df_rmsd.to_csv("rmsd.csv")
df_radius.to_csv("radius.csv")





        
        
figure, ax1 = plt.subplots(1,1)

ax1.plot(time, rmsd_kinase_left, color=biotite.colors["dimorange"], label="left")
ax1.plot(time, rmsd_kinase_right, color=biotite.colors["darkgreen"], label="right")
ax1.set_xlim(time[0], time[-1])
ax1.set_ylim(0, np.max([rmsd_kinase_left.max(), rmsd_kinase_right.max()])*1.1)
ax1.set_title("FRB-FRB Dimer")
ax1.set_xlabel("Time (ns)")
ax1.set_ylabel("RMSD (Å)")
ax1.set_xticks(np.arange(min(time), max(time)+1, 25))

legend_x = 1
legend_y = 0.5
plt.legend(loc='center left', bbox_to_anchor=(legend_x, legend_y))
plt.tight_layout()
plt.savefig("rmsd_two.png", dpi=600)
plt.clf()


r_min = np.min(
    np.array(
        [
            radius_kinase_left.min(),
            radius_kinase_right.min()                        
        ]
    )
)

r_max = np.max(
    np.array(
        [
            radius_kinase_left.max(),
            radius_kinase_right.max()                        
        
        ]
    )
)


figure, ax1 = plt.subplots(1,1)


ax1.plot(time, radius_kinase_left,  color=biotite.colors["dimorange"], label="left")
ax1.plot(time, radius_kinase_right, color=biotite.colors["darkgreen"], label="right")
ax1.set_ylim(r_min*0.8, r_max*1.1)
ax1.set_title("FRB-FRB Dimer")
ax1.set_xlabel("Time (ns)")
ax1.set_ylabel("Radius (Å)")

legend_x = 1
legend_y = 0.5
plt.legend(loc='center left', bbox_to_anchor=(legend_x, legend_y))
#plt.legend(loc="best")
plt.tight_layout()
plt.savefig("radius_two.png", dpi=600)
plt.clf()


ca_trajectory_kinase_left = trajectory_kinase_left[
    :, 
    trajectory_kinase_left.atom_name == "CA"
]
ca_trajectory_kinase_right = trajectory_kinase_right[
    :, 
    trajectory_kinase_right.atom_name == "CA"
]


rmsf_kinase_left = struc.rmsf(
    struc.average(ca_trajectory_kinase_left),
    ca_trajectory_kinase_left
)
rmsf_upper_kinase_left = rmsf_kinase_left.max()*1.1


rmsf_kinase_right = struc.rmsf(
    struc.average(ca_trajectory_kinase_right),
    ca_trajectory_kinase_right
)
rmsf_upper_kinase_right = rmsf_kinase_right.max()*1.1


fig, (ax1,ax2) = plt.subplots(2,1)

res_count = struc.get_residue_count(trajectory_kinase_left)
ax1.plot(
    np.arange(1, res_count+1) + 2801, 
    rmsf_kinase_left,
    color=biotite.colors["dimorange"]
)
ax1.set_title("Kinase Left")
#ax1.axvline(3828, ls="--", color="k")
#ax1.axvline(3838, ls="--", color="k")
ax1.set_xlim(2801+1, 2801+res_count)
ax1.set_ylim(0, rmsf_upper_kinase_left)
ax1.set_xlabel("Residue")
ax1.set_ylabel("RMSF (Å)")

res_count = struc.get_residue_count(trajectory_kinase_right)
ax2.plot(
    np.arange(1, res_count+1) + 2801, 
    rmsf_kinase_right,
    color=biotite.colors["dimorange"]
)
ax2.set_title("Kinase Right")
ax2.set_xlim(2801+1, 2801+res_count)
ax2.set_ylim(0, rmsf_upper_kinase_right)
#ax2.axvline(3828, ls="--", color="k")
#ax2.axvline(3838, ls="--", color="k")
ax2.set_xlabel("Residue")
ax2.set_ylabel("RMSF (Å)")

plt.tight_layout()
plt.savefig("rmsf.png")
plt.clf()

exit()

# averaging
trajectory_stable=trajectory[ind_stable:]


time_stable = time[ind_stable:]
inds = np.arange(0, time_stable.shape[0], time_stable.shape[0]/100).astype(int)
times_movie = time_stable[inds]
np.save("times_movie", times_movie)
        
make_dir("trajectory_movie_pdbs")
for i in inds:
    strucio.save_structure(
        "trajectory_movie_pdbs/frame_" + str(i) + ".pdb", 
         trajectory[i, trajectory.atom_name == 'CA']
    )

avg_pdb = struc.average(trajectory_stable)
strucio.save_structure("average.pdb", avg_pdb)
ind_closest = np.argmin(struc.rmsd(avg_pdb.coord, trajectory_stable.coord))
strucio.save_structure("frame_closest_avg.pdb", trajectory[ind_closest])

traj_avg = struc.average(trajectory_stable)
strucio.save_structure("average_data.pdb", traj_avg)










