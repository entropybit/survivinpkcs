#!/bin/bash
ROSETTA_PATH=/opt/rosetta_src_2018.33.60351_bundle
FOLDER=$HOME/nfs/DNAPKcs
RS_EXEC=$ROSETTA_PATH/main/source/bin
RS_DATA=$ROSETTA_PATH/main/database
NCORES=80
N=1
#PREPACK=docking_prepack_protocol.mpi.linuxgccrelease
DOCK=docking_protocol.default.linuxgccrelease


mkdir $FOLDER/dimer_refined/
ls $FOLDER/sur_dimer | grep _0001_ | grep .pdb | parallel $RS_EXEC/$DOCK -docking:partners BC_A -docking_local_refine -nstruct 1 -database $RS_DATA -score:docking_interface_score 1 -native $FOLDER/pk_mono_sur_di.pdb --overwrite -out:prefix $FOLDER/dimer_refined/ -out:file:scorefile dimer_scores.sc -s $FOLDER/sur_dimer/{} ::: ${FILES[@]} 



