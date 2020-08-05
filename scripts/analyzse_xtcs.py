import biotite
import biotite.structure as struc
import biotite.structure.io as strucio
import biotite.structure.io.xtc as xtc
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import scipy
from scipy.spatial import distance_matrix
from multiprocessing import Pool


def f_worker(line):

    X, Y = line

    d = distance_matrix(X, Y)



    return np.linalg.norm(d), np.min(d), np.max(d), np.mean(d)


def distance(traj1, traj2):
    """
        calculates distance between all particles of traj1 
        and traj2 in parallel using f_worker and parallel Pool
        returns an array of distance matrices                
    """

    assert(traj1.coord.shape[0] == traj2.coord.shape[0])
    Ts = range(traj1.shape[0])       
    lines = [[traj1[t].coord, traj2[t].coord] for t in Ts]
    distance_matrices = None
    


    with Pool() as p:
        distance_matrices = p.map(f_worker, lines)


    return np.array(distance_matrices)

    
if __name__=="__main__":

    

    # take first output structure as structure to align to    
    template_dimer = strucio.load_structure(
        "dimer_refined/pk_mono_sur_di_0001_000001_0001.pdb"
    )

    print(" ... loading XTC files ... ")
    xtc_dimer = xtc.XTCFile()
    xtc_dimer.read("dimers_ordered_by_cleaned.xtc")#, 1, 10)
    print(" ... done ... ")
    print("")
    print("")

    trajectory_dimer = xtc_dimer.get_structure(template_dimer)


    pkcs_start = trajectory_dimer[0][trajectory_dimer[0].chain_id == 'A']
    pkcs_start = pkcs_start[pkcs_start.atom_name == "CA"]
    trajectory_dimer, transform = struc.superimpose(
        trajectory_dimer[0], 
        trajectory_dimer,
        (trajectory_dimer[0].chain_id == 'A') 
        & (trajectory_dimer[0].atom_name == 'CA')
    )

    trajectory_dimer_ca = trajectory_dimer[:, 
        (trajectory_dimer.atom_name == "CA")
        & ((trajectory_dimer.res_id < 3206) | (trajectory_dimer.res_id > 3226))
    ]
    trajectory_dimer_activesite = trajectory_dimer[:,
        (trajectory_dimer.res_id >= 3747)
        & (trajectory_dimer.res_id <= 4015)

    ]
    trajectory_dimer_survivin_1 = trajectory_dimer[:,
        trajectory_dimer.chain_id == 'B'
    ]
    trajectory_dimer_survivin_2 = trajectory_dimer[:,
        trajectory_dimer.chain_id == 'C'
    ]

    # superimpose on kinase
    trajectory1, transform1 = struc.superimpose(
        surv_dnapk_dimer_ca1, 
        trajectory_dimer_ca,
        trajectory_dimer_ca.chain_id == 'A'
    )

    # superimpose on kinase
    trajectory2, transform2 = struc.superimpose(
        surv_dnapk_dimer_ca2,
        trajectory_dimer_ca,
        trajectory_dimer_ca.chain_id == 'A'    
    )

    # modify this accordingly or just create a subfolder contacts
    contacts_folder= "contacts/"
    residues_selection = list(np.arange(15,88))


    for resid in residues_selection:

        print(" ... distance calculations for S" +str(resid) + " ... ")

        d_resid_site1 = distance(
            trajectory_dimer_activesite,
            trajectory_dimer_survivin_1[:,    
                (trajectory_dimer_survivin_1.res_id == resid)
                & ~np.isin(
                    trajectory_dimer_survivin_1.atom_name,
                    (
                    "CA", "N", "C", "O", "1H", "H",
                    "2H", "3H", "OXT", "HA",
                    )
                )
            ]    
        )

        np.save(
            contacts_folder+"d_activesite_S"+str(resid)+ "_dimer_1", 
            d_resid_site1
        )

        d_resid_site2 = distance(
            trajectory_dimer_activesite,
            trajectory_dimer_survivin_2[:,
                (trajectory_dimer_survivin_2.res_id == resid)
                & ~np.isin(
                    trajectory_dimer_survivin_2.atom_name,
                    (
                    "CA", "N", "C", "O", "1H", "H",
                    "2H", "3H", "OXT", "HA",
                    )
                )            
            ]
        )
        np.save(
            contacts_folder + "d_activesite_S"+str(resid)+"_dimer_2", 
            d_resid_site2
        )


   

