#!/bin/bash

OUT_NAME=$1
JOBS=$2

gmx_mpi grompp -f prod.mdp -c npt.gro -t npt.cpt -p topol.top -o $OUT_NAME""_md.tpr
mpirun -np $JOBS gmx_mpi mdrun -deffnm $OUT_NAME""_md -v
