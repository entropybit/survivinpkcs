#!/bin/bash

OUT_NAME=$1

echo "1 \n 0 \n" | gmx trjconv -s $OUT_NAME""_md.tpr -f $OUT_NAME""_md.xtc -o $OUT_NAME""_md_noPBC.xtc -pbc mol -center
echo "4 \n 4 \n" | gmx rms -s $OUT_NAME""_md.tpr -f $OUT_NAME""_md_noPBC.xtc -o rmsd.xvg -tu ns
echo "1\n" | gmx gyrate -s $OUT_NAME""_md.tpr -f $OUT_NAME""_md_noPBC.xtc -o gyrate.xvg
