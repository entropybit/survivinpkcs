#!/bin/bash
ROSETTA_PATH=/opt/rosetta_src_2018.33.60351_bundle
FOLDER=$HOME/nfs/DNAPKcs
RS_EXEC=$ROSETTA_PATH/main/source/bin
RS_DATA=$ROSETTA_PATH/main/database
NCORES=48
N=100000
PREPACK=docking_prepack_protocol.mpi.linuxgccrelease
DOCK=docking_protocol.mpi.linuxgccrelease

DOCK_DIMER=true
DOCK_MONOMER=false

#### DOCKING ####################################

$RS_EXEC/$PREPACK -docking:partners BC_A -docking::sc_min -s $FOLDER/pk_mono_sur_di.pdb -database $RS_DATA
mkdir $FOLDER/sur_dimer
mv $FOLDER/pk_mono_sur_di_0001.pdb $FOLDER/sur_dimer/
cd $FOLDER/sur_dimer
mpiexec --hostfile $FOLDER/hostfile -n $NCORES $RS_EXEC/$DOCK -docking:partners BC_A -randomize1 -randomize2 -spin -dock_pert 3 8 -low_res_protocol_only -nstruct $N -native $FOLDER/pk_mono_sur_di.pdb -s $FOLDER/sur_dimer/pk_mono_sur_di_0001.pdb -database $RS_DATA -score:docking_interface_score 1
cd $FOLDER

################################################
