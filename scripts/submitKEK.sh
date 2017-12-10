#!/bin/bash

#MSGDIR=/dev/null
MSGDIR=$PWD/root

runNo=1012
N=493190
runname="syml4.i0"
prerunname="syml4.im1"
#geoSetup=1 # 1: finger; 0: normal
geoSetup=0 # 1: finger; 0: normal

nHitsMax=16
t0shift=0
tmin=-20
tmax=800
sumCut=30
aaCut=45
debug=-1

QUEUE="h"

for i in 1 2 3 4 5 6 7 8;
do
    testlayer=$i
    for (( j=0; j<N; j+=10000 ))
    do
        iEntryStart=$j
        iEntryStop=$((j+9999))
        if (( iEntryStop>N ))
        then
            iEntryStop=$((N-1))
        fi
        trunname="${runname}.$iEntryStart-$iEntryStop"
        log="$PWD/root/t_${runNo}.${trunname}.layer${testlayer}.log"
        err="$PWD/root/t_${runNo}.${trunname}.layer${testlayer}.err"
        bsub -o $MSGDIR -e $MSGDIR -q $QUEUE "tracking $runNo $testlayer $trunname $nHitsMax $t0shift $tmin $tmax $geoSetup $sumCut $aaCut $iEntryStart $iEntryStop $debug $prerunname > $log 2> $err"
    done
done
