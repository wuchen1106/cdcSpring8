#!/bin/bash

StartName="original"
runNo="117"
runNocon="117"
nEvents="605460"
runName="0126"
IterStart=15
IterEnd=50
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

lastxtfile="xt.117.0126.i11.root"

isReady(){ # make sure this thread is not processing any job (the conf is emtpy)!
    conf=$1
    log=$2
    echo "test" > $conf
    sleep 11
    echo "" > $conf
    if [ -e $log ]
    then
        message=`tail -n 1 $log`
        if [ ! -z "$message" ]
        then
            if echo $message | grep -q "Reset this job"
            then
                return 0
            else
                return 1
            fi
        else
            return 2
        fi
    else
        return 3
    fi
}

checkThread(){
    tempname=$1
    tlayer=$2
    iStart=$3
    iStop=$4
    conf="$CDCS8WORKING_DIR/Conf/${runNocon}.layer${tlayer}.${iStart}-${iStop}.conf"
    log="$CDCS8WORKING_DIR/Conf/${runNocon}.layer${tlayer}.${iStart}-${iStop}.log"
    if [ ! -e $conf ]
    then
        echo "    ERROR! configure file doesn't exist! $conf"
        return 6
    fi
    configure=`cat $conf`
    if [ ! -z "$configure" ] # thread with configure
    then
        if [ ! -e $log ] # but no log file?
        then
            sleep 3 # wait
            if [ ! -e $log ] # really no log file!
            then
                echo "    ERROR! configure file not empty but no log file! $conf"
                echo "        $configure"
                return 4 # job dead!?
            fi
        fi
        NPARAS=`echo $configure | gawk '{print NF;}'`
        if [ $NPARAS -eq 16 ] # configure format is correct
        then
            trunno=`echo $configure | gawk '{print $1;}'`
            ttlayer=`echo $configure | gawk '{print $2;}'`
            ttempname=`echo $configure | gawk '{print $4;}'`
            tiStart=`echo $configure | gawk '{print $12;}'`
            tiStop=`echo $configure | gawk '{print $13;}'`
            if [ $trunno == $runNo -a $ttempname == $tempname -a $ttlayer == $tlayer -a $tiStart == $iStart -a $tiStop == $iStop ]
            then
                return 1 # under processing
            else
                return 2 # used by others
            fi
        else # configure format is wrong but not reset?
            echo "    ERROR! wrong format in $conf"
            echo "        $configure"
            if [ -e $log ]; then tail $log; fi
            return 5 # job dead!?
        fi
    else  # thread with no configure, probably not processing any job
        isReady $conf $log
        if [ $? -eq 0 ]
        then
            return 0 # ready to process
        else
            return 3 # Waiting in line
        fi
    fi
}

findVacentThread(){
    for (( i=0; i<3600; i++ )) #3600*10 sec = 10 hours running
    do
        for tlayer in 4 5 3 6 2 7 1 8;
        do
            for (( j=0; j<nEvents; j+=5000 ))
            do
                iStart=$j
                iStop=$((j+4999))
                if (( iStop>=nEvents ))
                then
                    iStop=$((nEvents-1))
                fi
                conf="$CDCS8WORKING_DIR/Conf/${runNocon}.layer${tlayer}.${iStart}-${iStop}.conf"
                log="$CDCS8WORKING_DIR/Conf/${runNocon}.layer${tlayer}.${iStart}-${iStop}.log"
                if [ ! -e $conf ]
                then
                    continue
                fi
                configure=`cat $conf`
                if [ -z "$configure" ] # thread with no configure, probably not processing any job
                then
                    isReady $conf $log
                    if [ $? -eq 0 ]
                    then
                        echo $conf
                        return 0 # ready to process
                    else # waiting in line
                        continue
                    fi
                else # dealing with other jobs
                    continue
                fi
            done
        done
        sleep 10
    done
    echo "Still cannot find a vacent thread!"
    exit 1
}

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
        tmax=500
    else
        tmin=-10
        tmax=220
    fi

    # one layer or more
    if [ $iter -gt 45 ]
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
            iEntryStart=$j
            iEntryStop=$((j+4999))
            if (( iEntryStop>=nEvents ))
            then
                iEntryStop=$((nEvents-1))
            fi
            jobname="${runNo}.${currunname}.$iEntryStart-$iEntryStop.layer${testlayer}"
            ((Njobs++))
            echo "checking job \"$jobname\""
            file="root/t_${jobname}.log"
            if [ -e $file ]
            then # log file already exists??
                if tail -n 3 $file | grep -q "Good Events" # finished
                then
                    echo "Job $jobname already finished!"
                    continue # no need to work on it
                else  # probably still under process
                    echo "Job $jobname is running by someone else!"
                    tail -n 1 $file
                    continue # don't know who is running it but ignore this job anyway
                fi
            else
                echo "$file doesn't exist, so generate a new job!"
            fi
            theConf="$CDCS8WORKING_DIR/Conf/${runNocon}.layer${testlayer}.${iEntryStart}-${iEntryStop}.conf"
            temprunname="${currunname}.$iEntryStart-$iEntryStop"
            tempconfig="$runNo $testlayer $prerunname $temprunname $nHitsMax $t0shift $tmin $tmax $geoSetup $sumCut $aaCut $iEntryStart $iEntryStop $workType $inputType $debug"
            checkThread ${temprunname} ${testlayer} ${iEntryStart} ${iEntryStop}
            result=$?
            if [ $result -eq 0 ] # thread is ready
            then # go with this job
                echo "Job $jobname is ready"
                echo "$tempconfig" > $theConf # send the trigger info to the job
            elif [ $result -eq 1 ] # thread is running this job
            then # do nothing
                echo "Job $jobname is already running!"
                tail $theConf
            elif [ $result -eq 2 ] # thread is used by others
            then # notify and find a job which is emtpy
                echo "Job $jobname is used by others!"
                tail $theConf
                echo "    will find another job to use..."
                theConf=`findVacentThread`
                echo "    found new job available $theConf"
                echo "$tempconfig" > $theConf # send the trigger info to the job
            elif [ $result -eq 3 ] # thread is still waiting
            then # notify and find a job which is emtpy
                echo "Job $jobname is still waiting in queue!"
                echo "    will find another job to use..."
                theConf=`findVacentThread`
                echo "    found new job available $theConf"
                echo "$tempconfig" > $theConf # send the trigger info to the job
            else # error with this thread
                echo "Job $jobname is dead?"
                echo "    will find another job to use..."
                theConf=`findVacentThread`
                echo "    found new job available $theConf"
                echo "$tempconfig" > $theConf # send the trigger info to the job
            fi
        done
    done
    echo "Iteration $iter, $Njobs to be finished!"

    allfinished=false
    for (( i=0; i<3600; i++ )) #3600*10 sec = 10 hours running
    do
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
            allfinished=true
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
        sleep 10
    done
    if [ ! $allfinished ]
    then
        echo "Still not finished!"
        exit 1
    fi

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
