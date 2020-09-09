#!/bin/bash

PDB_NAME=$1
OUT_NAME=$2

gmx_mpi pdb2gmx -f $PDB_NAME -o $OUT_NAME"".gro -p topol.top -ff charmm36-mar2019 -water tip3p -ignh
gmx_mpi editconf -f $OUT_NAME"".gro -o $OUT_NAME""_boxed.gro -c -d 0.7 -bt dodecahedron
gmx_mpi solvate -cp $OUT_NAME""_boxed.gro -cs spc216.gro -o $OUT_NAME""_solv.gro -p topol.top
gmx_mpi grompp -f ions.mdp -c $OUT_NAME""_solv.gro -p topol.top -o ions.tpr -maxwarn 1
echo "13" | gmx_mpi genion -s ions.tpr -o $OUT_NAME""_ions.gro -p topol.top -pname NA -nname CL -neutral
