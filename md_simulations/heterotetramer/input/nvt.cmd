#!/bin/bash

ENERGY_DATA=$1
JOBS=$2


gmx_mpi grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr
mpirun --np $JOBS gmx_mpi mdrun -deffnm nvt -v
if [ $ENERGY_DATA == 1 ]
then
	echo "16\n0\n" | gmx energy -f nvt.edr -o temperature.xvg
fi
