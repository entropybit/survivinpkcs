

# Modeling of Survivin - PKCS - Survivin complexes.
This repository holds the scripts that were used for the modeling done 
for the publication *Spatial and functional interaction of a 
heterotetramer Survivin-DNA-PKcs complex in DNA double-strand break repair* 
(Submitted). 

## Dependencies.
Theses are necessary to use all the scripts in the complete workflow described in 
the next section. 
- [GNU parallel](https://www.gnu.org/software/parallel)
- [Rosetta Software Suite](https://www.rosettacommons.org/software)
- PyMol
- Gromacs 
- Python>=3.7
- [biotite](https://www.biotite-python.org/)

Although it will be explicitly mentioned, which script needs which software.

## Modeling step by step.
The basis for the here documented modeling steps were [1E31](https://www.rcsb.org/structure/1E31) for Survivin
and [5LUQ](https://www.rcsb.org/structure/5LUQ) for the DNA-PKcs.
After an initial loop modeling using [Modeller](https://salilab.org/modeller/) a closed structure was found for the
Kinase and FAT Domain of DNA-PKcs, which we called [head](https://github.com/entropybit/survivinpkcs/blob/master/pdbs/pkcs_head.pdb) structure in the paper. 
Furthermore, the here used [head](https://github.com/entropybit/survivinpkcs/blob/master/pdbs/pkcs_head.pdb) structure is the result of a 50ns productive run Molecular Dynamics Simulation in gromacs, after initial minimization an equilibration in NVT and NPT ensemble. The scripts for this can be found in the subfolder
[md_simulations/head_structur](https://github.com/entropybit/survivinpkcs/blob/master/md_simulations/head_structure/).

The justification for only using the head structure is as follows: Due to a large gap of missing structural information
that can not be closed with loop modeling, a complete model of the DNA-PKcs structure is not available. The head structure consists of
the Kinase and the FAT and FAT-C domain, which together make up the bulk of the protein surrounding the PI3K region.
As for this work only the kinase activity is of interest, it is argued that this region contains enough structural information
to sufficiently model the DNA-PKcs.

From a perliminary docking with the Schrödinger suite, the best two poses were available, which are stored as [schroedinger_pose1](https://github.com/entropybit/survivinpkcs/blob/master/pdbs/reworked/schroedinger_pose1.pdb) and [schroedinger_pose2](https://github.com/entropybit/survivinpkcs/blob/master/pdbs/reworked/schroedinger_pose2.pdb). After cutting away the cradle region of the
DNA-PKcs as well as an additional Protein docked to Survivin, which is not consieder in this study and should be thought of as artefact one get's
the according strcutures [schroedinger_pose1_aligned](https://github.com/entropybit/survivinpkcs/blob/master/pdbs/reworked/schroedinger_pose1_aligned.pdb) and 
[schroedinger_pose1_aligned](https://github.com/entropybit/survivinpkcs/blob/master/pdbs/reworked/schroedinger_pose1_aligned.pdb). Visualizing both, aligned on the head structure yields the following picture

<img src="https://github.com/entropybit/survivinpkcs/blob/master/pdbs/reworked/differing_states.png" width="800">

which is also stored as a [PyMol session](https://github.com/entropybit/survivinpkcs/blob/master/pdbs/reworked/differing_states.pse).
Here, we see very different positions of Survivin relative to its docking partner, placing two different Survivin residues, W67 and S20 in close proximity to
the PI3K region. These poses can not be both fulfilled at the same time.
However, mutational studies on the Survivin residues show that W67 and S20 are important residues for the Survivin - DNA-PKcs interaction. since mutation these
lead to the biggest change in Kinase activity.
Therefore the idea is, that both poses are fulfilled by Survivin binding in according poses to two
DNA-PKcs Proteins, on each of its monomers. 
An [initial structure](https://github.com/entropybit/survivinpkcs/blob/master/pdbs/kinase_heterodimer_aligned_0001.pdb) for this was generated by aligning both poses on Survivin accordingly. 

<img src="https://github.com/entropybit/survivinpkcs/blob/master/pdbs/aligned_heterotetramer.png" width="800">

To achieve this we first replaced the DNA-PKcs head region from the Schroedinger poses with the closed head structure, resulting
in the two structures [state1](https://github.com/entropybit/survivinpkcs/blob/master/pdbs/state1.pdb) and [state2](https://github.com/entropybit/survivinpkcs/blob/master/pdbs/state2.pdb). However, this thas the effect that proximity for S20 to the PI3k region
is not given anymore. A Molecular Dynamics simuation starting from this structure did not result in a stable simuation instead the two DNA-PKcs - Survivin Monomer complexes seperated. 

To overcome this and also to analyze the importance of the BIR region with a larger sample size, a global docking was performed with Rosetta.

### Protein-Protein Docking with Rosetta

First we perforemd a global docking, which means that a complete uninformed docking is run where several positions of the *ligand* protein on the surface of the *receptor* protein are tried. This was done using our global_docking.sh skript which can simply be executed 
after making it executable:
```
chmod u+x global_docking.sh
./global_docking.sh
```
Within this skript global docking of the survivin dimer with the head domain is executed. 
Before using it the Rosetta paths and target path have to be updated:
```
ROSETTA_PATH=/path/to/rosetta
FOLDER=/path/to/execution/folder
```

Furthermore, the hostfile needs to be adapted to ones cluster if necessary.
In our case the calculation has been evaluated on 16 cores per machine over 3 cluster nodes. 
In general the hostfile should look like this
```
192.168.0.11 slots=64
192.168.0.12 slots=64
...
```
with the ips of the according nodes. Hostnames can also be used as long as these can be resolved from each 
single node.

The global docking script generates a folder *sur_dimer* inside the working direcotry of the script.
Inside this folder all the resulting structures from the global docking are stored. As the global 
docking uses a rough backbone protocoll these structures do not have sidechains for the mobile structure.
To remedy that, and also to increase the resolution of the docking overall a local refinement docking
is done for every single one of these structures. Simliar to global docking there is a script for
this task in the repo. So one can simply do
```
chmod u+x refinement_docking.sh
./refinement_docking.sh
```
This is trivially parallelized over the files using gnu parallel and the results are stored in the new folder *dimer_refined*.
