#!/bin/bash

#MSGDIR=/dev/null
MSGDIR=$PWD/root

#runNo=117
runNo=1012
testlayer=5
runname=sym
nHitsMax=16
t0shift=0
tmin=-20
tmax=800
#geoSetup=1 # 1: finger; 0: normal
geoSetup=0 # 1: finger; 0: normal
sumCut=30
aaCut=45
iEntryStart=0
iEntryStop=10000
debug=-1
suf="i0"

QUEUE="h"

#for nHitsMax in 10 11 12 13 14 15 16
#do
    for i in 3 4 5 6 1 2 7 8;
#    for i in 5
    do
        testlayer=$i
#        for (( j=0; j<600000; j+=10000 ))
        for (( j=0; j<490000; j+=10000 ))
        do
            iEntryStart=$j
            iEntryStop=$((j+9999))
            suffix="${suf}n${nHitsMax}.$iEntryStart-$iEntryStop"
            log="$PWD/root/t_${runNo}.${runname}.${suffix}.layer${testlayer}.log"
            err="$PWD/root/t_${runNo}.${runname}.${suffix}.layer${testlayer}.err"
            bsub -o $MSGDIR -e $MSGDIR -q $QUEUE "tracking $runNo $testlayer $runname $nHitsMax $t0shift $tmin $tmax $geoSetup $sumCut $aaCut $iEntryStart $iEntryStop $debug $suffix > $log 2> $err"
        done
#        iEntryStart=600000
#        iEntryStop=605460
        iEntryStart=490000
        iEntryStop=493189
        suffix="${suf}n${nHitsMax}.$iEntryStart-$iEntryStop"
        log="$PWD/root/t_${runNo}.${runname}.${suffix}.layer${testlayer}.log"
        err="$PWD/root/t_${runNo}.${runname}.${suffix}.layer${testlayer}.err"
        bsub -o $MSGDIR -e $MSGDIR -q $QUEUE "tracking $runNo $testlayer $runname $nHitsMax $t0shift $tmin $tmax $geoSetup $sumCut $aaCut $iEntryStart $iEntryStop $debug $suffix > $log 2> $err"
    done
#don400000e
