#!/bin/bash

names=`ls Iter_*.p* | sed "s/Iter_\(.*\)\.i\w*\.layer\w*\.p\w*/\1/g"`
prename="NONE"
for name in $names
do
    if [ "$name" = "$prename" ]
    then
        continue
    fi
    prename=$name
    homeDir=$PWD
    runNo=`echo $name | sed "s/\(\w*\)\.\(.*\)/\1/g"`
    runName=`echo $name | sed "s/\(\w*\)\.\(.*\)/\2/g"`
    targetDir=$runNo/$runName
    for subdir in "pdf" "png"
    do
        mkdir -p $subdir/$targetDir
        mv *_$name.*.$subdir $subdir/$targetDir
        cd $subdir/$targetDir
        mkdir -p Iter; mv Iter_* Iter/;
        mkdir -p IterN; mv IterN_* IterN/;
        mkdir -p LRB; mv LRB_* LRB/;
        mkdir -p track; mv track_* track/;
        mkdir -p xtsamples; mv xtsamples_* xtsamples/;
        mkdir -p xtsamplesn; mv xtsamplesn_* xtsamplesn/;
        mkdir -p xtslices; mv xtslices_* xtslices/;
        mkdir -p xtslicesn; mv xtslicesn_* xtslicesn/;
        mkdir -p aaVSD; mv aaVSD_* aaVSD/;
        mkdir -p aaVST; mv aaVST_* aaVST/;
        mkdir -p dedx; mv dedx_* dedx/;
        mkdir -p DOCA; mv DOCA_* DOCA/;
        mkdir -p DriftD; mv DriftD_* DriftD/;
        mkdir -p effd; mv effd_* effd/;
        mkdir -p effx; mv effx_* effx/;
        mkdir -p ggVSX; mv ggVSX_* ggVSX/;
        mkdir -p offx; mv offx_* offx/;
        mkdir -p resd; mv resd_* resd/;
        mkdir -p resVSD; mv resVSD_* resVSD/;
        mkdir -p resVSX; mv resVSX_* resVSX/;
        mkdir -p resx; mv resx_* resx/;
        mkdir -p rmsd; mv rmsd_* rmsd/;
        mkdir -p rmsx; mv rmsx_* rmsx/;
        mkdir -p TX; mv TX_* TX/;
        mkdir -p XT; mv XT_* XT/;
        cd $homeDir
    done
done
