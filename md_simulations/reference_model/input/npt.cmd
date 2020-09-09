#!/bin/bash

ENERGY_DATA=$1
JOBS=$2

gmx_mpi grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr
mpirun --np $JOBS gmx_mpi mdrun -deffnm npt -v

if [ $ENERGY_DATA == 1 ]
then
	echo "18\n0\n" | gmx energy -f npt.edr -o pressure.xvg
	echo "24\n0\n" | gmx energy -f npt.edr -o density.xvg
fi
