#!/bin/bash

script="submitKEK.sh"
StartName="Garfield"
runNo="100002"
nEvents="90457"
runName="0111"
IterStart=1
IterEnd=3
layers="4"

geoType=0 # 0 for general; 1 for finger
inputType=1 # 1 for MC; 0 for data
workType=0 # 0, fr/l_0; 1, even/odd; -1, even/odd reversed; others, all layers
nHitsMax=13
t0shift=0
tmin=-20
tmax=800
sumCut=-20
aaCut=20
debug=-1

XTTYPE=0 # 2 for symmetrical; 0 for no constraints
DEBUG=-1
SAVEHISTS=0

for (( iter=IterStart; iter<=IterEnd; iter++ ))
do
    im1=$((iter-1))
    if [ $iter == 1 ]
    then
        name="${runName}.i${iter}"
        prename="${StartName}"
    else
        name="${runName}.i${iter}"
        prename="${runName}.i${im1}"
    fi
    eval $script -g $geoType -i $inputType -w $workType -n $nHitsMax -t $t0shift -a $tmin -z $tmax -s $sumCut -q $aaCut -d $debug $runNo $nEvents $prename $name $layers

    for (( i=0; i<3600; i++ )) #3600*10 sec = 10 hours running
    do
        sleep 10
        echo -n "$i "
        finished=true
        thefile="temp"
        for file in root/t_${runNo}.${name}.*.log
        do
            if ! tail -n 3 $file | grep -q "Good Events"
            then
                finished=false
                thefile=$file
            fi
        done
        if $finished
        then
            echo "Iteration $iter finished"
            break
        else
            echo -n $thefile:
            tail -n 1 $thefile
        fi
    done

    cd root/
    for ilayer in $layers
    do
        combine $runNo $name $ilayer
    done
    rm t_${runNo}.${name}.*-*.*
    cd ..
    getXT $runNo $prename $name $XTTYPE $geoType $SAVEHISTS $inputType $DEBUG
done
