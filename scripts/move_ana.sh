#!/bin/bash

names=`ls xtslicesn_*.pdf | sed "s/xtslicesn_\(.*\)\.i\w*\.layer\w*\.pdf/\1/g"`
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
    echo mkdir -p $targetDir
    echo mv *_$name.* $targetDir
    echo cd $targetDir
#    echo mkdir -p pdf
#    echo mv *.pdf pdf
#    echo mkdir -p png
#    echo mv *.png png
#    for dir in pdf png;
#    do
#        cd $dir
#        mkdir -p Iter; mv Iter_* Iter/;
#        mkdir -p IterN; mv IterN_* IterN/;
#        mkdir -p LRB; mv LRB_* LRB/;
#        mkdir -p track; mv track_* track/;
#        mkdir -p xtsamples; mv xtsamples_* xtsamples/;
#        mkdir -p xtsamplesn; mv xtsamplesn_* xtsamplesn/;
#        mkdir -p xtslices; mv xtslices_* xtslices/;
#        mkdir -p xtslicesn; mv xtslicesn_* xtslicesn/;
#        cd -
#    done
#    cd $homeDir
done
