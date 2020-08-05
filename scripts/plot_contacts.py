import numpy as np
import matplotlib as mpl
import itertools
mpl.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import os
from matplotlib.colors import LogNorm


folder = "../dimer_refined"


def create_dir(path):
    if not os.path.exists(path):
        os.makedirs(path)
        
        
def f_single(x):
    return 1/np.pi * np.arctan(x) + 0.5

def ranking_score(dist, I_sc, d_cutoff = 3):    
    return f_single(-I_sc)*np.exp(-1.0/200*dist)

    

#    return np.exp(-1.0*dist)*np.exp(-1.0*I_sc)
#    return np.exp(-0.1*(0.05*dist + I_sc))
    

def relative_pdb_path(s):
    return "../" + s.split("/")[-2] + "/" + s.split("/")[-1] + ".pdb"


def copy_bir_and_nonbir(scores):


    #printing and copying selection of non bir <-> activesite dockings
    cut_nonbir = scores[
        (scores['I_sc'] >= -17)
        & (scores['I_sc'] <= -16
        )
        & (scores['d_bir_activesite'] >= 19)    
        & (scores['d_bir_activesite'] <= 21)
    ]
    create_dir("nonbir_cluster_pdbs")
    print("best non bir binding models")
    print(cut_nonbir['description'])
    print("")

    i = 0
    for desc in cut_nonbir['description']:
        pdb_path = relative_pdb_path(desc)
        print("cp " + str(pdb_path) + " nonbir_cluster_pdbs/" + str(i) + ".pdb")
        os.system("cp " + str(pdb_path) + " nonbir_cluster_pdbs/" + str(i) + ".pdb")
        i = i+1
        
        
    print("")
    print("")    
    print("")
    print("")


    # printing and copying best BIR binding models
    cut_bestbir = scores[
        (scores['I_sc'] <= -27)
        & (scores['d_bir_activesite'] <= 4)
    ]

    indices = np.argsort(cut_bestbir['I_sc'])

    create_dir("bir_cluster_pdbs")
    print("best bir binding models ::")
    print(cut_bestbir['description'])
    print("")
    for index in indices:
        row = cut_bestbir.iloc[index]    
        desc = row['description']        
        pdb_path = relative_pdb_path(desc)
        print("cp " + str(pdb_path) + " bir_cluster_pdbs/" + str(np.abs(row['I_sc'])) + ".pdb")
        os.system("cp " + str(pdb_path) + " bir_cluster_pdbs/" + str(np.abs(row['I_sc'])) + ".pdb")
        
        
    print("")
    print("")    
    print("")
    print("")

def load_scores(cut_selector, folder="contacts/"):
    """
        loads the scores.csv and loads distance .npys for all residues in BIR
        domain.
        Selects the minimum distance from the two survivin columns and stores
        it in new column.
        Generates histograms for all pair scores on (I_sc, distance) per residue        
    
        todo: add parameters in function
    """

    rankings = []
    scores_cut = None
    d_cutoff = np.inf


    scores = pd.read_csv("../dnapkcs_dimersscores.csv")
    create_dir("plots") 

    create_dir("pdbs_d" +str(d_cutoff))


    for rid in residues_selection:
        print(" ... processing S" + str(rid) + " ... ")

        distances1 = np.load(folder+"d_activesite_S" + str(rid) + "_dimer_1.npy")
        distances1 = distances1[:, 1]



        distances2 = np.load(folder+"d_activesite_S" + str(rid) + "_dimer_2.npy")
        distances2 = distances2[:, 1]

        scores["S" + str(rid) + "_1"] = distances1
        scores["S" + str(rid) + "_2"] = distances2

        scores_cut = scores[
              (scores["S" + str(rid) + "_1"] <= d_cutoff)
            | (scores["S" + str(rid) + "_2"] <= d_cutoff)
        ]

        distances1 = scores_cut["S" + str(rid) + "_1"]
        distances2 = scores_cut["S" + str(rid) + "_2"]
        distances = scores_cut[["S" + str(rid) + "_1", "S" + str(rid) + "_2"]].min(axis=1)
        
        
        scores["S" + str(rid)] = scores[["S" + str(rid) + "_1", "S" + str(rid) + "_2"]].min(axis=1)
        
        
        
        I_sc = scores_cut["I_sc"]

    #    plt.hist2d(distances1, distances2, 50)
    #    plt.hist2d(distances1, distances2, bins=(100,100), norm=LogNorm())
    #    plt.colorbar()
    #    plt.xlabel("d(dimer_1, PI3K)")
    #    plt.ylabel("d(dimer_2, PI3K)")
    #    plt.title("S"+str(rid)+ " to PI3K active site")
    #    plt.savefig("plots/dimer_vs_S"+str(rid)+".png", dpi=600)   
    #    plt.clf()    

    #    plt.hist2d(distances1, I_sc, bins=(100,100))
    #    plt.colorbar()
    #    plt.xlabel("d(dimer_1, PI3K)")
    #    plt.ylabel("I_sc")
    #    plt.title("S"+str(rid)+ " to PI3K active site")
    #    plt.savefig("plots/score_dimer_1_S"+str(rid)+".png", dpi=600)  
    #    plt.clf()

        plt.hist2d(distances, I_sc, 50)
    #    plt.hist2d(distances, I_sc, bins=(100,100), norm=LogNorm())
        plt.colorbar()
        plt.xlabel("d(dimer, PI3K)")
        plt.ylabel("I_sc")
        plt.title("S"+str(rid)+ " to PI3K active site")
        plt.axvline(6.0, ls="--", color="r")
        plt.savefig("plots/score_dimer_S"+str(rid)+".png", dpi=600)  
        plt.clf()

        ranking = ranking_score(distances, I_sc)
        rankings.append(ranking)

        scores_cut["S" + str(rid) + "_ranking"] = np.array(ranking)
        scores["S" + str(rid) + "_ranking"] = np.array(ranking)


        indmax = scores_cut["S" + str(rid) + "_ranking"].argmax()

        print("       pdb for max ranked S" +str(rid) + " ::")
        print("           " + str(scores_cut['description'][indmax]))
        print(" min(I_sc) :: " + str(scores_cut['I_sc'][indmax]))
        print(" min(d) :: " + str(distances[indmax]))    
        
        pdb_path = relative_pdb_path(str(scores_cut['description'][indmax]))
        

        os.system("cp " + str(pdb_path) + " pdbs_d" +str(d_cutoff) + "/S" + str(rid) + ".pdb")

        print("       --> done")

    distances = scores[cut_selector].min(axis=1)
    scores["d_bir_activesite"] = distances
    scores.to_csv("scores_extended.csv")


    scores_cut = scores[scores["d_bir_activesite"] <= 20]
    scores_cut_large = scores[scores["d_bir_activesite"] <= 100]


    plt.hist2d(
        scores['I_sc'], scores["d_bir_activesite"], 
        bins=100, norm=LogNorm()
    )
    plt.colorbar()
    plt.xlabel("$I_{sc}$")
    plt.ylabel("$d_{BIR} [\AA]$")
    plt.yscale('log')
    #plt.xscale('log')
    plt.tight_layout()
    plt.savefig("hist2d_bir_activesite.png", dpi=600)
    plt.clf()

    plt.hist2d(
        scores_cut_large['I_sc'], scores_cut_large["d_bir_activesite"], 
        bins=100, norm=LogNorm()
    )
    plt.colorbar()
    plt.xlabel("$I_{sc}$")
    plt.ylabel("$d_{BIR} [\AA]$")
    plt.tight_layout()
    plt.savefig("hist2d_bir_activesite_large.png", dpi=600)
    plt.clf()


    plt.hist2d(
        scores_cut['I_sc'], scores_cut["d_bir_activesite"], 
        bins=100, norm=LogNorm()
    )
    plt.colorbar()
    plt.xlabel("$I_{sc}$")
    plt.ylabel("$d_{BIR} [\AA]$")
    plt.tight_layout()
    plt.savefig("hist2d_bir_activesite_cut.png", dpi=600)
    plt.clf()

    return scores

     
#residues_good_pairs = list(itertools(residues_good, 2))
#residues_good_triples = list(itertools(residues_good, 3))






#residues_selection = [
#    20,
#    25,
#    27,
#    28,
#    29,
#    31,
#    34,
#    53,
#    67,
#    76,
#    79,
#    117            
#]

residues_selection = list(np.arange(15,88))

residues_good = [
    67, 20, 27, 
    34, 53, 31
] 
residues_good_pairs = list(itertools.combinations(residues_good, 2))
residues_good_triples = list(itertools.combinations(residues_good, 3))



#weights_good = [
#    10, 10, 8,
#    7,  6,  5
#]

cut_selector = ["S" + str(rid) for rid in residues_selection]
cut_selector_good = ["S" + str(rid) for rid in residues_good]
cut_selector_good_pairs = [
    ["S" + str(rid) for rid in pair] 
    for pair in residues_good_pairs
]    
cut_selector_good_triples = [
    ["S" + str(rid) for rid in triple] 
    for triple in residues_good_triples
]    

scores = None
if "scores_extended.csv" in os.listdir("."):
    scores = pd.read_csv("scores_extended.csv")
else:
    scores = load_scores(cut_selector)

d_cutoff = 4





cut_selection = scores[scores['d_bir_activesite'] <=4]
cut_selection_numeric = cut_selection[['I_sc']+cut_selector_good]
cut_min = cut_selection_numeric.min(axis=0)
cut_max = cut_selection_numeric.max(axis=0)



cut_normed = (cut_selection_numeric - cut_min)/(cut_max - cut_min)




ranking_mean_all = cut_normed.mean(axis=1)
ranking_prod_all = cut_normed.prod(axis=1)/(len(cut_selector)+1.0)

ranking_mean_pairs = [
    cut_normed[['I_sc'] + selection].mean(axis=1) 
    for selection in cut_selector_good_pairs
]
ranking_prod_pairs = [
    cut_normed[['I_sc'] + selection].prod(axis=1)/(len(selection)+1.0)
    for selection in cut_selector_good_pairs
]

ranking_mean_triples = [
    cut_normed[['I_sc'] + selection].mean(axis=1) 
    for selection in cut_selector_good_triples
]
ranking_prod_triples = [
    cut_normed[['I_sc'] + selection].prod(axis=1)/(len(selection)+1.0)
    for selection in cut_selector_good_triples
]

create_dir("plot_ranks_pairs")
create_dir("pdbs_ranks_pairs")
#print(ranking_mean_pairs)


best_pairs = []

i = 0
for energy in ranking_mean_pairs:

#    print("energy: ")
#    print(energy)
#    print("")
#    print("")
#    print("i:")
#    print(i)
#    print("")
#    print("")
#    print("")
    
    x1, x2 = cut_selector_good_pairs[i]    
#    print(x1)
#    print(x2)
    ind = energy.idxmin()   
    pdb = relative_pdb_path(
        cut_selection.loc[ind,'description']
    )   

    cut_subselect = cut_selection[
        (cut_selection[str(x1)] <= 4.0)
        & (cut_selection[str(x2)] <= 4.0)
    ]

    if len(cut_subselect.index) > 0:

        indbest = cut_subselect['I_sc'].idxmin()
        pdb_best = relative_pdb_path(
            cut_subselect.loc[indbest,'description']
        )  

        best_pairs.append([
            cut_subselect.loc[indbest,'I_sc'], 
            float(x1[1:]), float(x2[1:])
        ])

        command = "cp " + str(pdb_best) + " pdbs_ranks_pairs/" + x1 + "_" + x2 
        command += ".pdb"    
        print(command)
        os.system(command) 

        print(" best for pair [" + str(x1) + ", " + str(x2) + " ] ::")
        print(cut_subselect.loc[indbest, ['I_sc', x1, x2, 'description']])
            
    command = "cp " + str(pdb) + " pdbs_ranks_pairs/" + x1 + "_" + x2 
    command += "_mean.pdb"    
    print(command)
    os.system(command)



    
    
    f, (ax1, ax2, ax3) = plt.subplots(3,1)
    ax1.hist2d(cut_selection['I_sc'], energy, bins=100, norm=LogNorm())
    ax1.set_xlabel("$I_{sc}$")
    ax1.set_ylabel("$\langle d_i \\rangle$")    
    ax2.hist2d(cut_selection[x1], energy, bins=100, norm=LogNorm())
    ax2.set_xlabel(x1)
    ax2.set_ylabel("$\langle d_i \\rangle$")
    ax3.hist2d(cut_selection[x2], energy, bins=100, norm=LogNorm())
    ax3.set_xlabel(x2)
    ax3.set_ylabel("$\langle d_i \\rangle$")
    plt.tight_layout()
    plt.savefig("plot_ranks_pairs/" + x1 + "_" + x2 + "_mean.png", dpi=600)
    plt.clf()      
    i = i+1 
       
#for energy, i in np.ndenumerate(ranking_prod_pairs):
i = 0
for energy in ranking_prod_pairs:

    
    x1, x2 = cut_selector_good_pairs[i]
    ind = energy.idxmin()   
    pdb = relative_pdb_path(
        cut_selection.loc[ind,'description']
    )           
    command = "cp " + str(pdb) + " pdbs_ranks_pairs/" + x1 + "_" + x2 
    command += "_prod.pdb"       
    print(command)
    os.system(command)


    
    
    f, (ax1, ax2, ax3) = plt.subplots(3,1)
    ax1.hist2d(cut_selection['I_sc'], energy, bins=100, norm=LogNorm())
    ax1.set_xlabel("$I_{sc}$")
    ax2.set_ylabel("$\langle d_i \\rangle$")    
    ax2.hist2d(cut_selection[x1], energy, bins=100, norm=LogNorm())
    ax2.set_xlabel(x1)
    ax2.set_ylabel("$\langle d_i \\rangle$")
    ax3.hist2d(cut_selection[x2], energy, bins=100, norm=LogNorm())
    ax3.set_xlabel(x2)
    ax3.set_ylabel("$\langle d_i \\rangle$")
    plt.tight_layout()
    plt.savefig("plot_ranks_pairs/" + x1 + "_" + x2 + "_prod.png", dpi=600)
    plt.clf()  
    i = i+1     
    
#    

best_triples = []
create_dir("plot_ranks_triples")
create_dir("pdbs_ranks_triples")
#for energy, i in np.ndenumerate(ranking_mean_triples):
i = 0
for energy in ranking_mean_triples:
    
    x1, x2, x3 = cut_selector_good_triples[i]    
    ind = energy.idxmin()   
    pdb = relative_pdb_path(
        cut_selection.loc[ind,'description']
    )           
    command  = "cp " + str(pdb) + " pdbs_ranks_triples/" 
    command += x1 + "_" + x2 + "_" + x3
    command += "_mean.pdb"    
    print(command)
    os.system(command)


    cut_subselect = cut_selection[
        (cut_selection[str(x1)] <= d_cutoff)
        & (cut_selection[str(x2)] <= d_cutoff)
        & (cut_selection[str(x3)] <= d_cutoff)
    ]

    if len(cut_subselect.index) > 0:

        indbest = cut_subselect['I_sc'].idxmin()
        pdb_best = relative_pdb_path(
            cut_subselect.loc[indbest,'description']
        )   

        best_triples.append([
            cut_subselect.loc[indbest,'I_sc'], float(x1[1:]), 
            float(x2[1:]), float(x3[1:])
        ])

        command  = "cp " + str(pdb_best) + " pdbs_ranks_triples/" 
        command += x1 + "_" + x2 + "_" + x3
        command += ".pdb"    
        print(command)
        os.system(command)
        print(" best for triple [" + str(x1) + ", " + str(x2) + ", " + str(x3) + " ] ::")
        print(cut_subselect.loc[indbest, ['I_sc', x1, x2, x3, 'description']])
    
    f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2,2)
    ax1.hist2d(cut_selection['I_sc'], energy, bins=100, norm=LogNorm())
    ax1.set_xlabel("$I_{sc}$")
    ax2.set_ylabel("$\langle d_i \\rangle$")    
    ax2.hist2d(cut_selection[x1], energy, bins=100, norm=LogNorm())
    ax2.set_xlabel(x1)
    ax2.set_ylabel("$\langle d_i \\rangle$")
    ax3.hist2d(cut_selection[x2], energy, bins=100, norm=LogNorm())
    ax3.set_xlabel(x2)
    ax3.set_ylabel("$\langle d_i \\rangle$")
    ax4.hist2d(cut_selection[x3], energy, bins=100, norm=LogNorm())
    ax4.set_xlabel(x3)
    ax4.set_ylabel("$\langle d_i \\rangle$")

    plt.tight_layout()

    plt.savefig(
        "plot_ranks_triples/" + x1 + "_" + x2 + "_" + x3 + "_mean.png", dpi=600
    )
    plt.clf()    
    i = i+1   

i = 0
for energy in ranking_prod_triples:
#for energy, i in np.ndenumerate(ranking_prod_pairs):
    
    x1, x2, x3 = cut_selector_good_triples[i]    
    ind = energy.idxmin()   
    pdb = relative_pdb_path(
        cut_selection.loc[ind,'description']
    )           
    command = "cp " + str(pdb) + " pdbs_ranks_triples/"
    command += x1 + "_" + x2 + "_" + x3    
    command += "_prod.pdb"       
    print(command)
    os.system(command)
    
    
    f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2,2)
    ax1.hist2d(cut_selection['I_sc'], energy, bins=100, norm=LogNorm())
    ax1.set_xlabel("$I_{sc}$")
    ax1.set_ylabel("$\langle d_i \\rangle$")    
    ax2.hist2d(cut_selection[x1], energy, bins=100, norm=LogNorm())
    ax2.set_xlabel(x1)
    ax2.set_ylabel("$\langle d_i \\rangle$")
    ax3.hist2d(cut_selection[x2], energy, bins=100, norm=LogNorm())
    ax3.set_xlabel(x2)
    ax3.set_ylabel("$\langle d_i \\rangle$")
    ax4.hist2d(cut_selection[x3], energy, bins=100, norm=LogNorm())
    ax4.set_xlabel(x3)
    ax4.set_ylabel("$\langle d_i \\rangle$")
    plt.tight_layout()
    plt.savefig(
        "plot_ranks_triples/" + x1 + "_" + x2 + "_" + x3 + "_prod.png", dpi=600
    )
    plt.clf()  
    i = i+1     


#print("ranking ::")
#print(ranking)
#print("")
#print("")

#plt.scatter(cut_selection['I_sc'], ranking)
plt.hist2d(cut_selection['I_sc'], ranking_mean_all, bins=100, norm=LogNorm())
plt.colorbar()
plt.xlabel('$I_{sc}$')
plt.ylabel("$\langle d_i \\rangle$")
plt.savefig("ranking_mean_all_Isc.png", dpi=600)
plt.clf()

f, (ax1, ax2) = plt.subplots(2,1)
ax1.hist2d(cut_selection['S67'], ranking_mean_all, bins=100, norm=LogNorm())
ax1.set_xlabel('$d_{67}$')
ax1.set_ylabel("$\langle d_i \\rangle$")
ax2.hist2d(cut_selection['S20'], ranking_mean_all, bins=100, norm=LogNorm())
ax2.set_xlabel('$d_{20}$')
ax2.set_ylabel("$\langle d_i \\rangle$")
plt.tight_layout()
plt.savefig("ranking_mean_all_dists.png", dpi=600)
plt.clf()

plt.hist2d(cut_selection['I_sc'], ranking_prod_all, bins=100, norm=LogNorm())
plt.colorbar()
plt.ylabel("$ \prod d_i / N $")
plt.xlabel('$I_{sc}$')
plt.savefig("ranking_prod_all.png", dpi=600)
plt.clf()


f, (ax1, ax2) = plt.subplots(2,1)
ax1.hist2d(cut_selection['S67'], ranking_prod_all, bins=100, norm=LogNorm())
ax1.set_xlabel('$d_{67}$')
ax1.set_ylabel("$ \prod d_i / N $")
ax2.hist2d(cut_selection['S20'], ranking_prod_all, bins=100, norm=LogNorm())
ax2.set_xlabel('$d_{20}$')
ax2.set_ylabel("$ \prod d_i / N $")
plt.tight_layout()
plt.savefig("ranking_prod_all_dists.png", dpi=600)
plt.clf()


best_pairs = np.array(best_pairs)
best_triples = np.array(best_triples)

pairs_ind_sorted = np.argsort(best_pairs[:,0])
triples_ind_sorted = np.argsort(best_triples[:,0])

print(" best pairs ::")
print(best_pairs[pairs_ind_sorted])
print("")

print(" best triples ::")
print(best_triples[triples_ind_sorted])
print("")


scores.to_csv("scores_extended_final.csv")



    


