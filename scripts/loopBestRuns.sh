#!/bin/bash
if [ $# -lt 2 ]
then
    echo "Please provide RunListFile and a tag"
    exit 0
fi

bestRuns=$1
runTag=$2
averageEtrack=1 # 1 means take averaged Etrack; 0 means using Etrack VS DOCA function
if [ $# -gt 2 ]
then
    averageEtrack=$3
fi

NRUNS=`cat $CDCS8WORKING_DIR/$bestRuns | wc -l`
for (( i=1; i<=$NRUNS; i++ ))
do
    runNo=`gawk -v gawk_keyword=$i 'BEGIN{i=1;}{if (i==gawk_keyword) print $3; i++;}' $bestRuns`
    runName=`gawk -v gawk_keyword=$i 'BEGIN{i=1;}{if (i==gawk_keyword) print $4"."$5; i++;}' $bestRuns`
    aaCut=`echo $runName | sed 's/.*a\(\w*\)n\(\w*\)\.i\w*/\1/'`
    CHI2MAX=2
    NMAX=`echo $runName |  sed 's/.*a\(\w*\)n\(\w*\)\.i\w*/\2/'`
    SLZMAX="0.1"

    found=false;
    for setuplist in $CDCS8WORKING_DIR/Input/list.C*.*;
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
    setup=`echo $setuplist | sed 's/.*\/list\.\(.*\)/\1/'`
    GAS=`echo $setup | sed 's/\(\w*\)\.\(\w*\)/\1/'`
    HV=`echo $setup | sed 's/\(\w*\)\.\(\w*\)/\2/'`
    if [ $GAS = "C2H6" -a $HV -lt 2000 ]
    then
        NMIN=5
    elif [ $GAS = "C4H10" -a $HV -lt 1650 ]
    then
        NMIN=5
    elif [ $GAS = "CH4" -a $HV -lt 1850 ]
    then
        NMIN=5
    else
        NMIN=7
    fi

    doIter.res.sh $runNo $runTag 0 400 1 5 $runName $aaCut $averageEtrack $CHI2MAX $NMIN $NMAX $SLZMAX > $CDCS8WORKING_DIR/res.$runNo.0426resxAV 2>&1 &
    pids="$!"
    wait $pids || { echo "there were errors in combining $runNo $currunname $ilayer" >&2; exit 1; }
done
