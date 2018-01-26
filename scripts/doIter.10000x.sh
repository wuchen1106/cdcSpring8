#!/bin/bash

StartName="Garfield"
runNo="100007"
runNocon="100007"
nEvents="90457"
runName="0125"
IterStart=1
IterEnd=70
layers="4"
wires=""

geoSetup=0 # 0 for general; 1 for finger
inputType=1 # 1 for MC; 0 for data
workType=0 # 0, fr/l_0; 1, even/odd; -1, even/odd reversed; others, all layers
nHitsMax=13
t0shift=0
#tmin=-10
#tmax=340
tmin=-10
tmax=800
sumCut=-20
aaCut=20
debug=-1

stepSize=0 # maximum step size for each movement in wire position calibration; 0 means no limit
minslz=0 # min slz cut for mean slz value in each sample of events in wiremap calibration. minslz==maxslz==0 means no cut
maxslz=0 # min slz cut for mean slz value in each sample of events in wiremap calibration. minslz==maxslz==0 means no cut
mininx=0 # min inx cut for mean inx value in each sample of events in wiremap calibration. mininx==maxinx==0 means no cut
maxinx=0 # min inx cut for mean inx value in each sample of events in wiremap calibration. mininx==maxinx==0 means no cut
maxchi2=1
scale=1 # move scale*offset on wiremap for fitting in the next round
XTTYPE=1 # 2 for symmetrical; 1 for offset loading; 0 for no constraints.
WPTYPE=0 # 0 for changing wiremap; 1 for not changing it;
UPDATEXT=1
DEBUG=-1
SAVEHISTS=0

lastxtfile=""
for (( iter=IterStart; iter<=IterEnd; iter++ ))
do
    im1=$((iter-1))
    if [ $iter == 1 ]
    then
        currunname="${runName}.i${iter}"
        prerunname="${StartName}"
    else
        currunname="${runName}.i${iter}"
        prerunname="${runName}.i${im1}"
    fi

    # t range
    if [ $iter -gt 6 ]
    then
        tmin=-10
        tmax=600
    elif [ $iter -gt 4 ]
    then
        tmin=-10
        tmax=360
    elif [ $iter -gt 2 ]
    then
        tmin=-10
        tmax=350
    else
        tmin=-10
        tmax=340
    fi

    # one layer or more
    if [ $iter -gt 67 ]
    then
        layers="1 2 3 4 5 6 7 8"
        WPTYPE=1 # 1 for changing wiremap
        UPDATEXT=1
        wires=""
        stepSize="0" # no step size limit
        scale="0.5"
        minslz="0"
        maxslz="0"
        maxinx="0"
        mininx="0"
    elif [ $iter -gt 62 ]
    then
        layers="3 4 5 6"
        WPTYPE=1 # 1 for changing wiremap
        UPDATEXT=1
        wires=""
        stepSize="0" # no step size limit
        scale="0.5"
        minslz="0"
        maxslz="0"
        maxinx="0"
        mininx="0"
    elif [ $iter -gt 59 ]
    then
        layers="7"
        WPTYPE=1 # 1 for changing wiremap
        UPDATEXT=0
        wires=""
        stepSize="0" # no step size limit
        scale="0.5"
        minslz="0"
        maxslz="0"
        maxinx="0"
        mininx="0"
    elif [ $iter -gt 56 ]
    then
        layers="2"
        WPTYPE=1 # 1 for changing wiremap
        UPDATEXT=0
        wires=""
        stepSize="0" # no step size limit
        scale="0.5"
        minslz="0"
        maxslz="0"
        maxinx="0"
        mininx="0"
    elif [ $iter -gt 53 ]
    then
        layers="1"
        WPTYPE=1 # 1 for changing wiremap
        UPDATEXT=0
        wires=""
        stepSize="0" # no step size limit
        scale="0.5"
        minslz="0"
        maxslz="0"
        maxinx="0"
        mininx="0"
    elif [ $iter -gt 50 ]
    then
        layers="8"
        WPTYPE=1 # 1 for changing wiremap
        UPDATEXT=0
        wires=""
        stepSize="0" # no step size limit
        scale="0.5"
        minslz="0"
        maxslz="0"
        maxinx="0"
        mininx="0"
    elif [ $iter -gt 45 ]
    then
        layers="1 2 3 4 5 6 7 8"
        WPTYPE=1 # 1 for changing wiremap
        UPDATEXT=1
        wires=""
        stepSize="0" # no step size limit
        scale="0.5"
        minslz="0"
        maxslz="0"
        maxinx="0"
        mininx="0"
    elif [ $iter -gt 40 ]
    then
        layers="3 4 5 6"
        WPTYPE=1 # 1 for changing wiremap
        UPDATEXT=1
        wires=""
        stepSize="0" # no step size limit
        scale="0.5"
        minslz="0"
        maxslz="0"
        maxinx="0"
        mininx="0"
    elif [ $iter -gt 35 ]
    then
        layers="1 2 7 8"
        WPTYPE=1 # 1 for changing wiremap
        UPDATEXT=0
        wires=""
        stepSize="0" # no step size limit
        scale="0.5"
        minslz="0"
        maxslz="0"
        maxinx="0"
        mininx="0"
    elif [ $iter -gt 32 ]
    then
        layers="8"
        WPTYPE=1 # 1 for changing wiremap
        UPDATEXT=0
        wires=""
        stepSize="0" # no step size limit
        scale="0.5"
        minslz="0"
        maxslz="0"
        maxinx="0"
        mininx="0"
    elif [ $iter -gt 29 ]
    then
        layers="1"
        WPTYPE=1 # 1 for changing wiremap
        UPDATEXT=0
        wires=""
        stepSize="0" # no step size limit
        scale="0.5"
        minslz="0"
        maxslz="0"
        maxinx="0"
        mininx="0"
    elif [ $iter -gt 26 ]
    then
        layers="7"
        WPTYPE=1 # 1 for changing wiremap
        UPDATEXT=0
        wires=""
        stepSize="0" # no step size limit
        scale="0.5"
        minslz="0"
        maxslz="0"
        maxinx="0"
        mininx="0"
    elif [ $iter -gt 23 ]
    then
        layers="2"
        WPTYPE=1 # 1 for changing wiremap
        UPDATEXT=0
        wires=""
        stepSize="0" # no step size limit
        scale="0.5"
        minslz="0"
        maxslz="0"
        maxinx="0"
        mininx="0"
    elif [ $iter -gt 20 ]
    then
        layers="3 4 5 6"
        WPTYPE=1 # 1 for changing wiremap
        UPDATEXT=1
        wires=""
        stepSize="0" # no step size limit
        scale="0.5"
        minslz="0"
        maxslz="0"
        maxinx="0"
        mininx="0"
    elif [ $iter -gt 17 ]
    then
        layers="6"
        WPTYPE=1 # 1 for changing wiremap
        UPDATEXT=0
        wires=""
        stepSize="0" # no step size limit
        scale="0.5"
        minslz="0"
        maxslz="0"
        maxinx="0"
        mininx="0"
    elif [ $iter -gt 14 ]
    then
        layers="3"
        WPTYPE=1 # 1 for changing wiremap
        UPDATEXT=0
        wires=""
        stepSize="0" # no step size limit
        scale="0.5"
        minslz="0"
        maxslz="0"
        maxinx="0"
        mininx="0"
    elif [ $iter -gt 11 ]
    then
        layers="5"
        WPTYPE=1 # 1 for changing wiremap
        UPDATEXT=0
        wires=""
        stepSize="0" # no step size limit
        scale="0.5"
        minslz="0"
        maxslz="0"
        maxinx="0"
        mininx="0"
    elif [ $iter -gt 8 ]
    then
        layers="4"
        WPTYPE=1 # 1 for changing wiremap
        UPDATEXT=1
        wires=""
        stepSize="0" # no step size limit
        scale="0.5"
        minslz="0"
        maxslz="0"
        maxinx="0"
        mininx="0"
    else
        layers="4"
        WPTYPE=0 # 0 for not changing wiremap
        UPDATEXT=1
    fi

    Njobs=0
    for testlayer in $layers;
    do
        for (( j=0; j<nEvents; j+=5000 ))
        do
            ((Njobs++))
            iEntryStart=$j
            iEntryStop=$((j+4999))
            if (( iEntryStop>=nEvents ))
            then
                iEntryStop=$((nEvents-1))
            fi
            temprunname="${currunname}.$iEntryStart-$iEntryStop"
            theConf="$CDCS8WORKING_DIR/Conf/${runNocon}.layer${testlayer}.${iEntryStart}-${iEntryStop}.conf"
            echo "$runNo $testlayer $prerunname $temprunname $nHitsMax $t0shift $tmin $tmax $geoSetup $sumCut $aaCut $iEntryStart $iEntryStop $workType $inputType $debug" > $theConf # send the trigger info to the job
        done
    done
    echo "Iteration $iter, $Njobs to be finished!"

    for (( i=0; i<3600; i++ )) #3600*10 sec = 10 hours running
    do
        sleep 10
        echo -n "$i "
        finished=true
        thefile=""
        NjobsFinished=0
        for file in root/t_${runNo}.${currunname}.*.log
        do
            if ! tail -n 3 $file | grep -q "Good Events"
            then
                finished=false
                thefile=$file
            else
                ((NjobsFinished++))
            fi
        done
        if [ $NjobsFinished -lt $Njobs ]
        then
            finished=false
        fi
        if $finished
        then
            echo "Iteration $iter finished"
            break
        else
            echo -n "$NjobsFinished/$Njobs jobs finihsed, "
            if [ -z $thefile ]
            then
                echo ""
            else
                echo -n $thefile:
                tail -n 1 $thefile
            fi
        fi
    done

    cd root/
    pids=""
    for ilayer in $layers
    do
        combine $runNo $currunname $ilayer &
        pids+=" $!"
    done
    wait $pids || { echo "there were errors" >&2; exit 1; }
    rm -f t_${runNo}.${currunname}.*-*.*
    cd ..

#   upgrading wireposition?
    if [ ! $WPTYPE -eq 1 ] # ! 1 for not changing wiremap
    then
        cd info
        ln -s wire-position.${runNo}.${StartName}.root wire-position.${runNo}.${currunname}.root
        cd ..
    fi
    getOffset $runNo $prerunname $currunname $geoSetup $WPTYPE $scale $stepSize $minslz $maxslz $mininx $maxinx $maxchi2 $DEBUG $wires
    if [ ! $UPDATEXT -eq 1 ] # ! 1 for not updating xt
    then
        if [ -z $lastxtfile ]
        then
            lastxtfile=xt.${runNo}.${prerunname}.root
        fi
        cd info
        ln -s $lastxtfile xt.${runNo}.${currunname}.root
        cd ..
    else
        getXT $runNo $prerunname $currunname $XTTYPE $geoSetup $SAVEHISTS $inputType $maxchi2 $DEBUG
        lastxtfile=xt.${runNo}.${currunname}.root
    fi
done
