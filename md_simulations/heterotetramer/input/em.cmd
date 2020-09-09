#!/bin/bash

ENERGY_DATA=$1
OUT_NAME=$2
JOBS=$3


gmx_mpi grompp -f min.mdp -c $OUT_NAME""_ions.gro -p topol.top -o em.tpr
mpirun -np $JOBS --hostfile hostfile gmx_mpi mdrun -v -deffnm em

if [ $ENERGY_DATA == 1 ]
then
	echo "10\n0\n" | gmx energy -f em.edr -o potential.xvg
fi
