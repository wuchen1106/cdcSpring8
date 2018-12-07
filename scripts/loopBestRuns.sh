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
CHI2MAX=2
if [ $# -gt 3 ]
then
    CHI2MAX=$4
fi
SLZMAX="0.1"
if [ $# -gt 4 ]
then
    SLZMAX=$5
fi

NRUNS=`cat $CDCS8WORKING_DIR/$bestRuns | wc -l`
for (( i=1; i<=$NRUNS; i++ ))
do
    runNo=`gawk -v gawk_keyword=$i 'BEGIN{i=1;}{if (i==gawk_keyword) print $3; i++;}' $bestRuns`
    runName=`gawk -v gawk_keyword=$i 'BEGIN{i=1;}{if (i==gawk_keyword) print $4"."$5; i++;}' $bestRuns`
    runFile=`gawk -v gawk_keyword=$i 'BEGIN{i=1;}{if (i==gawk_keyword) print $6; i++;}' $bestRuns`
    LAYER=`echo $runFile | sed 's/.*layer\(\w\).*/\1/'`
    NMAX=`echo $runName |  sed 's/.*a\(\w*\)n\(\w*\)\.i\w*/\2/'`
    runTagthis=${runTag}tl$LAYER

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

    doIter.res.sh $runNo $runTagthis 0 400 1 5 $runName $averageEtrack $CHI2MAX $NMIN $NMAX $SLZMAX $LAYER > $CDCS8WORKING_DIR/res.$runNo.$runTagthis 2>&1 &
    pids="$!"
    wait $pids || { echo "there were errors in doIter.res.sh $runNo $currunname $ilayer" >&2; exit 1; }
done
