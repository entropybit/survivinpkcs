
import biotite
import biotite.structure as struc
import biotite.structure.io as strucio
import biotite.structure.io.xtc as xtc
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
font = {'size'   : 16}
mpl.rc('font', **font)
import matplotlib.pyplot as plt
import pandas as pd
from biotite.structure.io import load_structure, save_structure
import os

def make_dir(directory):
    if not os.path.exists(directory):
            os.makedirs(directory)

# Put here the path of the downloaded files
templ_file_path = "heterodimer.pdb"
template_model = strucio.load_structure(templ_file_path)
#templ_file_path = "output/npt.gro"
#traj_file_path  = "output/kinase_dimer_md_center.xtc"
#traj_file_path  = "output/kinase_dimer_nowater_fit.xtc"
#traj_file_path  = "output/kinase_dimer_md.xtc"
traj_file_path = "output/heterodimer_nopbc_cluster_fit.xtc"

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

trajectory_kinase_left  = trajectory[:, trajectory.chain_id == "A"]
trajectory_kinase_right  = trajectory[:, trajectory.chain_id == "D"]
trajectory_survivin_left  = trajectory[:, trajectory.chain_id == "B"]
trajectory_survivin_right  = trajectory[:, trajectory.chain_id == "C"]


# Get simulation time for plotting purposes
time = xtc_file.get_time()
time = 10**(-3)*time

ind_stable = np.where(time == 20)[0][0]
#print(" stable after 75ns which is frame[" +str(ind_stable) + "] ")

print("start time is ::"+ str(time[0]))
print("end time is   :: "+ str(time[-1])) 




#trajectory_kinase_left = struc.remove_pbc(trajectory_kinase_left)
#trajectory_kinase_right = struc.remove_pbc(trajectory_kinase_right)
#trajectory_survivin_left = struc.remove_pbc(trajectory_survivin_left)
#trajectory_survivin_right = struc.remove_pbc(trajectory_survivin_right)

print(" ... writing start frame ...")


#kinase_left_start = trajectory_kinase_left[0]
#kinase_left_start.chain_id = np.repeat('A', len(kinase_left_start))

#survivin_left_start = trajectory_survivin_left[0]
#survivin_left_start.chain_id = np.repeat('B', len(survivin_left_start))

#survivin_right_start = trajectory_survivin_right[0]
#survivin_right_start.chain_id = np.repeat('C', len(survivin_right_start))
#
#kinase_right_start = trajectory_kinase_right[0]
#kinase_right_start.chain_id = np.repeat('D', len(kinase_right_start))


#frame_start = kinase_left_start + survivin_left_start
#frame_start += survivin_right_start + kinase_right_start



frame_start = template_model.copy()
#frame_start, _ = struc.superimpose(trajectory[0], frame_start)
frame_start.coord = trajectory[0].coord
save_structure("frame_start_coord.pdb", frame_start)
save_structure("frame_start.pdb", trajectory[0])
print(" ... done ... ")

print(" ... writing frame[1] ... ")

frame_1 = template_model.copy()
#frame_start, _ = struc.superimpose(trajectory[0], frame_start)
frame_1.coord = trajectory[1].coord
save_structure("frame_1_coord.pdb", frame_1)
save_structure("frame_1.pdb", trajectory[1])


print(" ... done ... ")

print(" ... writing end frame ...")
#kinase_left_end = trajectory_kinase_left[-1]
#kinase_left_end.chain_id = np.repeat('A', len(kinase_left_end))

#survivin_left_end = trajectory_survivin_left[-1]
#survivin_left_end.chain_id = np.repeat('B', len(survivin_left_end))

#survivin_right_end = trajectory_survivin_right[-1]
#survivin_right_end.chain_id = np.repeat('C', len(survivin_right_end))

#kinase_right_end = trajectory_kinase_right[-1]
#kinase_right_end.chain_id = np.repeat('D', len(kinase_right_end))
#frame_end += survivin_right_end + kinase_right_end

frame_end = template_model.copy()
#frame_end, _ = struc.superimpose(trajectory[-1], frame_end)
frame_end.coord = trajectory[-1].coord
save_structure("frame_end_coord.pdb", frame_end)
save_structure("frame_end.pdb", trajectory[-1])
print(" ... done ... ")


rmsd_overall = struc.rmsd(
    trajectory[0], trajectory
)
radius_overall = struc.gyration_radius(trajectory)





figure, (ax1, ax2) = plt.subplots(2,1)

ax1.plot(time, rmsd_overall, color=biotite.colors["dimorange"])
ax1.set_xlim(time[0], time[-1])
ax1.set_ylim(0, rmsd_overall.max()*1.1)
ax1.set_xlabel("Time (ns)")
ax1.set_ylabel("RMSD (Å)")
ax1.set_xticks(np.arange(min(time), max(time)+1, 25))

ax2.plot(time, radius_overall, color=biotite.colors["dimorange"])
ax2.set_xlim(time[0], time[-1])
ax2.set_ylim(0, radius_overall.max()*1.1)
ax2.set_xlabel("Time (ns)")
ax2.set_ylabel("Radius (Å)")
ax2.set_xticks(np.arange(min(time), max(time)+1, 25))

plt.tight_layout()
plt.savefig("overall.png", dpi=600)
plt.clf()


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

# survivin left
trajectory_survivin_left, transform = struc.superimpose(
    trajectory_survivin_left[0], trajectory_survivin_left
)
rmsd_survivin_left = struc.rmsd(
    trajectory_survivin_left[0], trajectory_survivin_left
)
radius_survivin_left = struc.gyration_radius(trajectory_survivin_left)

# survivin right
trajectory_survivin_right, transform = struc.superimpose(
    trajectory_survivin_right[0], trajectory_survivin_right
)
rmsd_survivin_right = struc.rmsd(
    trajectory_survivin_right[0], trajectory_survivin_right
)
radius_survivin_right = struc.gyration_radius(trajectory_survivin_right)



figure, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2,2)

ax1.plot(time, rmsd_kinase_left, color=biotite.colors["dimorange"])
ax1.set_xlim(time[0], time[-1])
ax1.set_ylim(0, rmsd_kinase_left.max()*1.1)
ax1.set_title("Left Kinase")
ax1.set_xlabel("Time (ns)")
ax1.set_ylabel("RMSD (Å)")
ax1.set_xticks(np.arange(min(time), max(time)+1, 25))
ax11 = ax1.twinx()
ax11.plot(time, radius_kinase_left, color="gray")
ax11.set_ylim(radius_kinase_left.min()*0.8, radius_kinase_left.max()*1.1)
ax11.set_ylabel("Radius (Å)")


ax3.plot(time, rmsd_survivin_left, color=biotite.colors["dimorange"])
ax3.set_xlim(time[0], time[-1])
ax3.set_ylim(0, rmsd_survivin_left.max()*1.1)
ax3.set_title("Left Survivin")
ax3.set_xlabel("Time (ns)")
ax3.set_ylabel("RMSD (Å)")
ax3.set_xticks(np.arange(min(time), max(time)+1, 25))
ax31 = ax3.twinx()
ax31.plot(time, radius_survivin_left, color="gray")
ax31.set_ylim(radius_survivin_left.min()*0.5,  radius_survivin_left.max()*1.1)
ax31.set_ylabel("Radius (Å)")

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
ax21.set_ylabel("Radius (Å)")


ax4.plot(time, rmsd_survivin_right, color=biotite.colors["dimorange"])
ax4.set_xlim(time[0], time[-1])
ax4.set_ylim(0, rmsd_survivin_right.max()*1.1)
ax4.set_title("Right Survivin")
ax4.set_xlabel("Time (ns)")
ax4.set_ylabel("RMSD (Å)")
ax4.set_xticks(np.arange(min(time), max(time)+1, 25))
ax41 = ax4.twinx()
ax41.plot(time, radius_survivin_right, color="gray")
ax41.set_ylim(radius_survivin_right.min()*0.8, radius_survivin_right.max()*1.1)
ax41.set_ylabel("Radius (Å)")

plt.tight_layout()
plt.savefig("rmsd_radius.png", dpi=600)
plt.clf()


rmsd_min = np.min(
    np.array([
        rmsd_survivin_left.min(),
        rmsd_survivin_right.min(),
    ])        
)
rmsd_max = np.max(
    np.array([
        rmsd_survivin_left.max(),
        rmsd_survivin_right.max(),
    ])        
    
)


r_surv_min = np.min(
    np.array([
        radius_survivin_left.min(),
        radius_survivin_right.min(),
    ])        
)
r_surv_max = np.max(
    np.array([
        radius_survivin_left.max(),
        radius_survivin_right.max(),
    ])        
    
)
r_kin_min = np.min(
    np.array([
        radius_kinase_left.min(),
        radius_kinase_right.min(),
    ])        
)
r_kin_max = np.max(
    np.array([
        radius_kinase_left.max(),
        radius_kinase_right.max(),
    ])        
    
)

figure, ((ax1, ax2)) = plt.subplots(2, 1)

ax1.plot(time, rmsd_kinase_left, color=biotite.colors["dimorange"], label="left")
ax1.plot(time, rmsd_kinase_right, color=biotite.colors["darkgreen"], label="righ")
ax1.set_xlim(time[0], time[-1])
ax1.set_ylim(0, rmsd_kinase_left.max()*1.1)
ax1.set_title("Kinase")
ax1.set_xlabel("Time (ns)")
ax1.set_ylabel("RMSD (Å)")
ax1.set_xticks(np.arange(min(time), max(time)+1, 25))


ax2.plot(time, rmsd_survivin_left, color=biotite.colors["dimorange"], label="left")
ax2.plot(time, rmsd_survivin_right, color=biotite.colors["darkgreen"], label="right")
ax2.set_xlim(time[0], time[-1])
ax2.set_ylim(0, rmsd_max*1.1)
ax2.set_title("Survivin")
ax2.set_xlabel("Time (ns)")
ax2.set_ylabel("RMSD (Å)")
ax2.set_xticks(np.arange(min(time), max(time)+1, 25))

legend_x = 1
legend_y = 0.5
plt.legend(loc='center left', bbox_to_anchor=(legend_x, legend_y))
plt.tight_layout()
plt.savefig("rmsd_two.png", dpi=600)
plt.clf()


# store data kinase
df_rmsd_kinase_right = pd.DataFrame(rmsd_kinase_right)
df_rmsd_kinase_left = pd.DataFrame(rmsd_kinase_left)
df_radius_kinase_right = pd.DataFrame(radius_kinase_right)
df_radius_kinase_left = pd.DataFrame(radius_kinase_left)

print(" ... writing dataframes ... ")

df_rmsd_kinase_right.to_csv("rmsd_kinase_right.csv")
df_rmsd_kinase_left.to_csv("rmsd_kinase_left.csv")
df_radius_kinase_right.to_csv("radius_kinase_right.csv")
df_radius_kinase_left.to_csv("radius_kinase_left.csv")

print(df_rmsd_kinase_right)
print(df_rmsd_kinase_left)
print(df_radius_kinase_right)
print(df_radius_kinase_left)

# store data survivin
df_rmsd_survivin_right = pd.DataFrame(rmsd_survivin_right)
df_rmsd_survivin_left = pd.DataFrame(rmsd_survivin_left)
df_radius_survivin_right = pd.DataFrame(radius_survivin_right)
df_radius_survivin_left = pd.DataFrame(radius_survivin_left)

df_rmsd_survivin_right.to_csv("rmsd_survivin_right.csv")
df_rmsd_survivin_left.to_csv("rmsd_survivin_left.csv")
df_radius_survivin_right.to_csv("radius_survivin_right.csv")
df_radius_survivin_left.to_csv("radius_survivin_left.csv")

df_rmsd = pd.DataFrame(rmsd_overall)
df_radius = pd.DataFrame(radius_overall)

df_rmsd.to_csv("rmsd.csv")
df_radius.to_csv("radius.csv")

print(" ... done ... ")



figure, ((ax1, ax2)) = plt.subplots(2, 1)

ax1.plot(time, radius_kinase_left, color=biotite.colors["dimorange"])
ax1.plot(time, radius_kinase_right, color=biotite.colors["darkgreen"])
ax1.set_xlim(time[0], time[-1])
ax1.set_ylim(0.8*r_kin_min, r_kin_max*1.1)
ax1.set_title("Kinase")
ax1.set_xlabel("Time (ns)")
ax1.set_ylabel("Radius (Å)")
ax1.set_xticks(np.arange(min(time), max(time)+1, 25))


ax2.plot(time, radius_survivin_left, color=biotite.colors["dimorange"], label="left")
ax2.plot(time, radius_survivin_right, color=biotite.colors["darkgreen"], label="right")
ax2.set_xlim(time[0], time[-1])
ax2.set_ylim(0.8*r_surv_min, r_surv_max*1.1)
ax2.set_title("Survivin")
ax2.set_xlabel("Time (ns)")
ax2.set_ylabel("Radius (Å)")
ax2.set_xticks(np.arange(min(time), max(time)+1, 25))

legend_x = 1
legend_y = 0.5
plt.legend(loc='center left', bbox_to_anchor=(legend_x, legend_y))
plt.tight_layout()
plt.savefig("radius_two.png", dpi=600)
plt.clf()



#figure, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2,2)

#ax1.plot(time, radius_kinase_left, color=biotite.colors["dimorange"])
#ax1.set_xlim(time[0], time[-1])
#ax1.set_ylim(0, radius_kinase_left.max()*1.1)
#ax1.set_title("Left Kinase")
#ax1.set_xlabel("Time (ns)")
#ax1.set_ylabel("Radius (Å)")

#ax3.plot(time, radius_survivin_left, color=biotite.colors["dimorange"])
#ax3.set_xlim(time[0], time[-1])
#ax3.set_ylim(0, radius_survivin_left.max()*1.1)
#ax3.set_title("Left Survivin")
#ax3.set_xlabel("Time (ns)")
#ax3.set_ylabel("Radius (Å)")

#ax2.plot(time, radius_kinase_right, color=biotite.colors["dimorange"])
#ax2.set_xlim(time[0], time[-1])
#ax2.set_ylim(0, radius_kinase_right.max()*1.1)
#ax2.set_title("Right Kinase")
#ax2.set_xlabel("Time (ns)")
#ax2.set_ylabel("Radius (Å)")

#ax4.plot(time, radius_survivin_right, color=biotite.colors["dimorange"])
#ax4.set_xlim(time[0], time[-1])
#ax4.set_ylim(0, radius_survivin_right.max()*1.1)
#ax4.set_title("Right Survivin")
#ax4.set_xlabel("Time (ns)")
#ax4.set_ylabel("Radius (Å)")

#plt.tight_layout()
#plt.savefig("radius.png", dpi=600)
#plt.clf()



ca_trajectory_kinase_left = trajectory_kinase_left[
    :, 
    trajectory_kinase_left.atom_name == "CA"
]
ca_trajectory_survivin_left = trajectory_survivin_left[
    :, 
    trajectory_survivin_left.atom_name == "CA"
]
ca_trajectory_kinase_right = trajectory_kinase_right[
    :, 
    trajectory_kinase_right.atom_name == "CA"
]
ca_trajectory_survivin_right = trajectory_survivin_right[
    :, 
    trajectory_survivin_right.atom_name == "CA"
]





rmsf_kinase_left = struc.rmsf(
    struc.average(ca_trajectory_kinase_left),
    ca_trajectory_kinase_left
)
rmsf_upper_kinase_left = rmsf_kinase_left.max()*1.1

rmsf_survivin_left = struc.rmsf(
    struc.average(ca_trajectory_survivin_left),
    ca_trajectory_survivin_left
)
rmsf_upper_survivin_left = rmsf_survivin_left.max()*1.1

rmsf_kinase_right = struc.rmsf(
    struc.average(ca_trajectory_kinase_right),
    ca_trajectory_kinase_right
)
rmsf_upper_kinase_right = rmsf_kinase_right.max()*1.1

rmsf_survivin_right = struc.rmsf(
    struc.average(ca_trajectory_survivin_right),
    ca_trajectory_survivin_right
)
rmsf_upper_survivin_right = rmsf_survivin_right.max()*1.1
#rmsf_to_cm = struc.rmsf(cm_candidate_ca, ca_trajectory)
#rmsf_upper = max([rmsf.max(), rmsf_to_cm.max()])*1.1


#fig, (ax1, ax2) = plt.subplots(2,1)

fig, ((ax1,ax2), (ax3,ax4)) = plt.subplots(2,2)

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

res_count = struc.get_residue_count(trajectory_survivin_left)
ax3.plot(
    np.arange(1, res_count+1), 
    rmsf_survivin_left,
    color=biotite.colors["dimorange"]
)
ax3.set_title("Survivin Left")
ax3.set_xlim(1, res_count)
ax3.axvline(67, ls="--", color="k")
ax3.axvline(20, ls="--", color="k")
ax3.set_ylim(0, rmsf_upper_survivin_left)
ax3.set_xlabel("Residue")
ax3.set_ylabel("RMSF (Å)")

res_count = struc.get_residue_count(trajectory_survivin_right)
ax4.plot(
    np.arange(1, res_count+1), 
    rmsf_survivin_right,
    color=biotite.colors["dimorange"]
)
ax4.set_title("Survivin Right")
ax4.set_xlim(1, res_count)
ax4.axvline(67, ls="--", color="k")
ax4.axvline(20, ls="--", color="k")
ax4.set_ylim(0, rmsf_upper_survivin_right)
ax4.set_xlabel("Residue")
ax4.set_ylabel("RMSF (Å)")

#ax2.plot(np.arange(1, res_count+1), rmsf_to_cm, color=biotite.colors["dimorange"])
#ax2.set_title("Relative to S_14771")
#ax2.set_xlim(1, res_count)
#ax2.set_ylim(0, rmsf_upper)
#ax2.set_xlabel("Residue")
#ax2.set_ylabel("RMSF (Å)")

plt.tight_layout()
plt.savefig("rmsf.png")
plt.clf()


figure, (ax1, ax2) = plt.subplots(2,1)


res_count = struc.get_residue_count(trajectory_kinase_left)
ax1.plot(
    np.arange(1, res_count+1) + 2801, 
    rmsf_kinase_left,
    color=biotite.colors["dimorange"]
)
ax1.set_title("Kinase Left")
ax1.axvline(3828, ls="--", color="k")
ax1.axvline(3838, ls="--", color="k")
ax1.set_xlim(2500, 2801+res_count)
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
ax2.set_xlim(2500, 2801+res_count)
ax2.set_ylim(0, rmsf_upper_kinase_right)
ax2.axvline(3828, ls="--", color="k")
ax2.axvline(3838, ls="--", color="k")
ax2.set_xlabel("Residue")
ax2.set_ylabel("RMSF (Å)")

plt.tight_layout()
plt.savefig("kinase_survivin_interact.png")
plt.clf()

#exit()


## plot distances S20 and S67 to active site
#s20 = trajectory[:
#    (trajectory.chain_id == 'B') 
#    & (trajectory.res_id == 20)
#]
#s20_com = s20.mean(axis=1)

#pi3k_left = trajectory[:,
#    (trajectory.chain_id == 'A')
#    & (trajectory.res_id >= 3747)
#    & (trajectory.res_id <= 4015)
#]    
#as_left = trajectory[:,
#    (trajectory.chain_id == 'A')
#    & (trajectory.res_id >= 3747)
#    & (trajectory.res_id <= 4015)
#]    

#s67 = trajectory[:
#    (trajectory.chain_id == 'C') 
#    & (trajectory.res_id == 67)
#]
#s67_com = s67.mean(axis=1)

#pi3k_right = trajectory[:,
#    (trajectory.chain_id == 'D')
#    & (trajectory.res_id >= 3747)
#    & (trajectory.res_id <= 4015)
#]    


#as_right = trajectory[:,
#    (trajectory.chain_id == 'A')
#    & (trajectory.res_id >= 3747)
#    & (trajectory.res_id <= 4015)
#]    





# averaging
trajectory_stable=trajectory[ind_stable:]


time_stable = time[ind_stable:]
inds = np.arange(0, time_stable.shape[0], time_stable.shape[0]/100).astype(int)
times_movie = time_stable[inds]
np.save("times_movie", times_movie)

pdbs_list = [

]
        
make_dir("trajectory_movie_pdbs")
for i in inds:
    strucio.save_structure(
        "trajectory_movie_pdbs/frame_" + str(i) + ".pdb", 
         trajectory[i, trajectory.atom_name == 'CA']
    )
    pdbs_list.append("frame_" + str(i) + ".pdb")

avg_pdb = struc.average(trajectory_stable)
strucio.save_structure("average.pdb", avg_pdb)
ind_closest = np.argmin(struc.rmsd(avg_pdb.coord, trajectory_stable.coord))
strucio.save_structure("frame_closest_avg.pdb", trajectory[ind_closest])

traj_avg = struc.average(trajectory_stable)
strucio.save_structure("average_data.pdb", traj_avg)


with open("pdb_list.txt", "w") as f:
    for name in pdbs_list:
        f.write(name + "\n")











