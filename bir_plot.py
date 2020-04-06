import numpy as np
import matplotlib as mpl
import itertools
mpl.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import os
from matplotlib.colors import LogNorm

font = {'family' : 'normal',
#        'weight' : 'bold',
        'size'   : 21}

mpl.rc('font', **font)

def load_scores(scores_path, cut_selector, folder):
    """
        loads the scores.csv and loads distance .npys for all residues in BIR
        domain.
        Selects the minimum distance from the two survivin columns and stores
        it in new column.
           
    """

    rankings = []
    scores_cut = None
    d_cutoff = np.inf
    scores = pd.read_csv(scores_path)
    
    listing = os.listdir(".")
    
    if "scores_extended.csv" in listing:
        scores = pd.read_csv("scores_extended.csv")
    else:        
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

        
        distances = scores[cut_selector].min(axis=1)
        scores["d_bir_activesite"] = distances
        scores.to_csv("scores_extended.csv")                 

    return scores
    
    
# generate a selection of residue ids
residues_selection = list(np.arange(15,88))
residues_good = [
    67, 20, 27, 
    34, 53, 31
] 
residues_good_pairs = list(itertools.combinations(residues_good, 2))
residues_good_triples = list(itertools.combinations(residues_good, 3))
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
    
# load scores with residue distances 
scores = load_scores(
    "dnapkcs_dimersscores.csv",  cut_selector,
    "contacts/"
    
)



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
