#!/bin/bash

keyword=""
if [ $# -gt 0 ]
then
    keyword=$1
fi

xtType=11
if [ $# -gt 1 ]
then
    xtType=$2
fi
maxChi2=2
if [ $# -gt 2 ]
then
    maxChi2=$3
fi
nHitsMaxini=0
if [ $# -gt 3 ]
then
    nHitsMaxini=$4
fi
saveHists=0
if [ $# -gt 4 ]
then
    saveHists=$5
fi

for file in root/ana_*${keyword}*.layer*.root
do
    runName=`echo $file | sed 's/root\/ana_\(\w*\)\.\(.*\)\.layer\(\w\)\.root/\2/'`
    runNo=`echo $file | sed 's/root\/ana_\(\w*\)\..*a\(\w*\)n\(\w*\).*\.i\w*.layer\(\w\)\.root/\1/'`
    aaCut=`echo $file | sed 's/root\/ana_\(\w*\)\..*a\(\w*\)n\(\w*\).*\.i\w*.layer\(\w\)\.root/\2/'`
    if [ $nHitsMaxini -eq 0 ]
    then
        nHitsMax=`echo $file | sed 's/root\/ana_\(\w*\)\..*a\(\w*\)n\(\w*\).*\.i\w*.layer\(\w\)\.root/\3/'`
    else
        nHitsMax=$nHitsMaxini
    fi
    layerID=`echo $file | sed 's/root\/ana_\(\w*\)\..*a\(\w*\)n\(\w*\).*\.i\w*.layer\(\w\)\.root/\4/'`
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
    dir=$CDCS8WORKING_DIR/results/ResEff/$setup/$runNo/layer${layerID}
    mkdir -p $dir
    cd $dir
    ana $runNo $runName $layerID $xtType $maxChi2 $nHitsMax $aaCut $saveHists
done
cd $CDCS8WORKING_DIR
