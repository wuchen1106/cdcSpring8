#!/bin/bash

script="submitKEK.sh"
StartName="Garfield"
runNo="1012"
nEvents="493189"
runName="0115"
IterStart=1
IterEnd=20
layers="4"
wires=""

geoType=0 # 0 for general; 1 for finger
inputType=0 # 1 for MC; 0 for data
workType=0 # 0, fr/l_0; 1, even/odd; -1, even/odd reversed; others, all layers
nHitsMax=13
t0shift=0
tmin=-10
tmax=340
sumCut=-20
aaCut=20
debug=-1

stepSize=0 # maximum step size for each movement in wire position calibration; 0 means no limit
scale=1 # move scale*offset on wiremap for fitting in the next round
XTTYPE=1 # 2 for symmetrical; 1 for offset loading; 0 for no constraints.
WPTYPE=0 # 0 for changing wiremap; 1 for not changing it;
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

#    # t range
#    if [ $iter -gt 6 ]
#    then
#        tmin=-10
#        tmax=600
#    elif [ $iter -gt 4 ]
#    then
#        tmin=-10
#        tmax=360
#    elif [ $iter -gt 2 ]
#    then
#        tmin=-10
#        tmax=350
#    else
#        tmin=-10
#        tmax=340
#    fi

    # one layer or more
    if [ $iter -gt 16 ]
    then
        layers="1 2 3 4 5 6 7 8"
        wires=""
        stepSize="0.02" # move 10 micrometer per iteration at most
        scale="0.8"
    elif [ $iter -gt 6 ]
    then
        layers="1 2 3 4 5 6 7 8"
        wires="4 5 6"
        stepSize="0.01" # move 10 micrometer per iteration at most
        scale="0.5"
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
        combine $runNo $name $ilayer &
        pids+=" $!"
    done
    wait $pids || { echo "there were errors" >&2; exit 1; }
    rm t_${runNo}.${name}.*-*.*
    cd ..

#   upgrading wireposition?
    if [ $iter -gt 6 ]
    then
        WPTYPE=1 # 1 for changing wiremap
    else
        WPTYPE=0 # 0 for not changing wiremap
        cd info
        ln -s wire-position.${runNo}.${StartName}.root wire-position.${runNo}.${name}.root
        cd ..
    fi
    getOffset $runNo $prename $name $geoType $WPTYPE $scale $stepSize $DEBUG $wires
    getXT $runNo $prename $name $XTTYPE $geoType $SAVEHISTS $inputType $DEBUG
done
