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
maxslz=0.1
if [ $# -gt 3 ]
then
    maxslz=$4
fi
nHitsMaxini=0
if [ $# -gt 4 ]
then
    nHitsMaxini=$5
fi
nHitsSMinini=7
if [ $# -gt 5 ]
then
    nHitsSMinini=$6
fi
geoSetup=0
if [ $# -gt 6 ]
then
    geoSetup=$7
fi
tmaxSet=0
if [ $# -gt 7 ]
then
    tmaxSet=$8
fi
saveHists=0
if [ $# -gt 8 ]
then
    saveHists=$9
fi
saveEventTree=0
if [ $# -gt 9 ]
then
    saveEventTree=$10
fi

for file in root/ana_*${keyword}*.layer*.root
do
    runName=`echo $file | sed 's/root\/ana_\(\w*\)\.\(.*\)\.layer\(\w\)\.root/\2/'`
    runNo=`echo $file | sed 's/root\/ana_\(\w*\)\..*a\(\w*\)n\(\w*\).*\.i\w*.layer\(\w\)\.root/\1/'`
    aaCut=`echo $file | sed 's/root\/ana_\(\w*\)\..*a\(\w*\)n\(\w*\).*\.i\w*.layer\(\w\)\.root/\2/'`
    if grep -q "_$runNo.$runName\." $CDCS8WORKING_DIR/bestRuns
    then
        echo "found $runNo.$runName"
    else
        echo "  cannot find $runNo.$runName"
        continue
    fi
    if [ $runNo -eq 115 -o $runNo -eq 117 ]
    then
        geoSetup=1
    else
        geoSetup=0
    fi
    if [ $nHitsMaxini -eq 0 ]
    then
        nHitsMax=`echo $file | sed 's/root\/ana_\(\w*\)\..*a\(\w*\)n\(\w*\).*\.i\w*.layer\(\w\)\.root/\3/'`
    else
        nHitsMax=$nHitsMaxini
    fi
    layerID=`echo $file | sed 's/root\/ana_\(\w*\)\..*a\(\w*\)n\(\w*\).*\.i\w*.layer\(\w\)\.root/\4/'`
    cd $CDCS8WORKING_DIR/Input
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
    GAS=`echo $setup | sed 's/\(\w*\)\.\(\w*\)/\1/'`
    HV=`echo $setup | sed 's/\(\w*\)\.\(\w*\)/\2/'`
    if [ $GAS = "C2H6" -a $HV -lt 2000 ]
    then
        nHitsSMin=5
    elif [ $GAS = "C4H10" -a $HV -lt 1650 ]
    then
        nHitsSMin=5
    elif [ $GAS = "CH4" -a $HV -lt 1850 ]
    then
        nHitsSMin=5
    else
        nHitsSMin=$nHitsSMinini
    fi
    dir=$CDCS8WORKING_DIR/result/ResEff/$setup/$runNo
    mkdir -p $dir
    cd $dir
    echo ana $runNo $runName $layerID $xtType $maxChi2 $maxslz $nHitsMax $nHitsSMin $aaCut $geoSetup $tmaxSet $saveHists $saveEventTree
    ana $runNo $runName $layerID $xtType $maxChi2 $maxslz $nHitsMax $nHitsSMin $aaCut $geoSetup $tmaxSet $saveHists $saveEventTree
done
cd $CDCS8WORKING_DIR
