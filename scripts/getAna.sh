#!/bin/bash

keyword=""
if [ $# -gt 0 ]
then
    keyword=$1
fi

xtType=11
maxChi2=1
saveHists=0

for file in root/ana_*${keyword}*.i15.layer*.root
do
    runName=`echo $file | sed 's/root\/ana_\(\w*\)\.\(.*\)\.layer\(\w\)\.root/\2/'`
    runNo=`echo $file | sed 's/root\/ana_\(\w*\)\..*a\(\w*\)n\(\w*\).*\.i15.layer\(\w\)\.root/\1/'`
    aaCut=`echo $file | sed 's/root\/ana_\(\w*\)\..*a\(\w*\)n\(\w*\).*\.i15.layer\(\w\)\.root/\2/'`
    nHitsMax=`echo $file | sed 's/root\/ana_\(\w*\)\..*a\(\w*\)n\(\w*\).*\.i15.layer\(\w\)\.root/\3/'`
    layerID=`echo $file | sed 's/root\/ana_\(\w*\)\..*a\(\w*\)n\(\w*\).*\.i15.layer\(\w\)\.root/\4/'`
    cd $CDCS8WORKING_DIR/info
    found=false;
    for setuplist in list.C*.*;
    do
        for rid in `cat $setuplist`
        do
            if [ $rid == $runNo ]
            then
                found=true
                break;
            fi
        done
        if $found
        then
            break;
        fi
    done
    setup=`echo $setuplist | sed 's/list\.\(.*\)/\1/'`
    dir=$CDCS8WORKING_DIR/results/ResEff/$setup/$runNo
    mkdir -p $dir
    cd $dir
    ana $runNo $runName $layerID $xtType $maxChi2 $nHitsMax $aaCut $saveHists
done
cd $CDCS8WORKING_DIR
