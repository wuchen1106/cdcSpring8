#!/bin/bash

script="submitKEK.sh"
StartName="Garfield"
runNo="1012"
nEvents="493189"
runName="0115"
IterStart=1
IterEnd=10
layers="4"

geoType=0 # 0 for general; 1 for finger
inputType=0 # 1 for MC; 0 for data
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

    # one layer or more
    if [ $iter -gt 3 ]
    then
        layers="3 4 5 6"
    elif [ $iter -gt 6 ]
        layers="2 3 4 5 6"
    elif [ $iter -gt 9 ]
        layers="1 2 3 4 5 6 7"
    else
        layers="4"
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

#   upgrading wireposition?
    if [ $iter -gt 2 ]
    then
        getOffset $runNo $prename $name $geoType $DEBUG
        XTTYPE=1 # 1 for offset loading
    else
        cp info/wire-position.${runNo}.${prename}.root info/wire-position.${runNo}.${name}.root
        XTTYPE=2 # 2 for symmetric without offset
    fi
    getXT $runNo $prename $name $XTTYPE $geoType $SAVEHISTS $inputType $DEBUG
done
