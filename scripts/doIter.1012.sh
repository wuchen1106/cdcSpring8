#!/bin/bash

threadName="job"
thread_iStart=0
thread_iStop=99

StartName="Garfield"
runNo="1012"
nEvents="493189"
nEvtPerRun="15000"
runName="0226.sm10a30n35"
IterStart=1
IterEnd=10
layers="4" # layers to be reconstructed and [analyzed (in case of layers is not 0)]
LAYERS="4" # layers to be analyzed (in case of layers is 0)
wires="" # wires to be calibrated (position)

# for tracking
geoSetup=0 # 0 for general; 1 for finger
inputType=0 # 1 for MC; 0 for data
workTypeini=0 # 0, fr/l_0; 1, even/odd; -1, even/odd reversed; others, all layers
nHitsGMaxini=13
t0shift0=0
t0shift1=0
tmin=-10
tmax=800
sumCut=-10
aaCut=30
peakType=0 # 0, only the first peak over threshold; 1, all peaks over threshold; 2, even including shaddowed peaks

# for getOffset
WPTYPE=0 # 0 for changing wiremap; 1 for not changing it;
stepSize=0 # maximum step size for each movement in wire position calibration; 0 means no limit
minslz=0 # min slz cut for mean slz value in each sample of events in wiremap calibration. minslz==maxslz==0 means no cut
maxslz=0 # min slz cut for mean slz value in each sample of events in wiremap calibration. minslz==maxslz==0 means no cut
mininx=0 # min inx cut for mean inx value in each sample of events in wiremap calibration. mininx==maxinx==0 means no cut
maxinx=0 # min inx cut for mean inx value in each sample of events in wiremap calibration. mininx==maxinx==0 means no cut
maxchi2=1
scale=1 # move scale*offset on wiremap for fitting in the next round

# for getXT
UPDATEXT=1
DEFAULTLAYER=4 # use this layer to generate fl(r)_0 and so on
XTTYPE=6 # 2 for symmetrical; 1 for symmetrical + offset loading; 0 for no constraints; 6 for symmetrical + offset loading + only first over threshold peak
NHITSMAXini=35
SAVEHISTS=0

threadLists=""
threadlistfile=threadlist.$runName.$runNo

lastxtfile=""

updateThreadLists(){
    ls Conf/${threadName}_*.log > $threadlistfile
    cat $threadlistfile
}

isReady(){ # make sure this thread is not processing any job (the conf is emtpy)!
    name="$1\.log"
    if echo $threadLists | grep -q "$name"
    then
        return 0
    else
        return 1
    fi
}

prev_ithread=$thread_iStart
prev_occupied=false
findVacentThread(){
    for (( i=0; i<3600; i++ )) #3600*10 sec = 10 hours running
    do
        if [ -e kill.$runNo.$runName ]
        then
            echo "Killed by user!"
            exit 0
        fi
        for (( ithread=prev_ithread; ithread<=thread_iStop; ithread++ ))
        do
            if $prev_occupied
            then
                prev_occupied=false
                continue;
            fi
            thethread=${threadName}_${ithread}
            conf="$CDCS8WORKING_DIR/Conf/${thethread}.conf"
            log="$CDCS8WORKING_DIR/Conf/${thethread}.log"
            if [ ! -e $conf ]
            then
                continue
            fi
            configure=`cat $conf`
            if [ -z "$configure" ] # thread with no configure, probably not processing any job
            then
                isReady $thethread
                if [ $? -eq 0 ]
                then
                    theConf=$conf
                    prev_ithread=$ithread
                    prev_occupied=true
                    return 0 # ready to process
                else # waiting in line
                    continue
                fi
            else # dealing with other jobs
                continue
            fi
        done
        prev_ithread=$thread_iStart
        prev_occupied=false
        sleep 10
        threadLists=`updateThreadLists`
        if [ ! $? -eq 0 ]
        then
            return 2 # cannot get hep_q
        fi
    done
    return 1 # cannot find any vacent slots in 10 hours
}

for (( iter=IterStart; iter<=IterEnd; iter++ ))
do
    im1=$((iter-1))
    if [ $iter == $IterStart ]
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
        tmax=800
    elif [ $iter -gt 4 ]
    then
        tmin=-10
        tmax=600
    elif [ $iter -gt 2 ]
    then
        tmin=-10
        tmax=360
    else
        tmin=-10
        tmax=340
    fi

    if [ $iter -gt 1 ]
    then
        workType=$workTypeini
    else
        workType=0
    fi

    if [ $iter -eq $IterEnd ]
    then
        NHITSMAX=0
        layers="1 2 3 4 5 6 7 8"
#        layers="0"
#        LAYERS="1 2 3 4 5 6 7 8"
    else
        nHitsGMax=$nHitsGMaxini
        NHITSMAX=$NHITSMAXini
    fi

    echo "#Iteration $iter started"
    echo "  layers = $layers"
    echo "  LAYERS = $LAYERS"
    echo "  wires = $wires"
    echo "  geoSetup = $geoSetup"
    echo "  inputType = $inputType"
    echo "  workType = $workType"
    echo "  nHitsGMax = $nHitsGMax"
    echo "  t0shift0 = $t0shift0"
    echo "  t0shift1 = $t0shift1"
    echo "  tmin = $tmin"
    echo "  tmax = $tmax"
    echo "  sumCut = $sumCut"
    echo "  aaCut = $aaCut"
    echo "  peakType = $peakType"
    echo "  stepSize = $stepSize"
    echo "  minslz = $minslz"
    echo "  maxslz = $maxslz"
    echo "  mininx = $mininx"
    echo "  maxinx = $maxinx"
    echo "  maxchi2 = $maxchi2"
    echo "  scale = $scale"
    echo "  XTTYPE = $XTTYPE"
    echo "  WPTYPE = $WPTYPE"
    echo "  UPDATEXT = $UPDATEXT"
    echo "  DEFAULTLAYER = $DEFAULTLAYER"
    echo "  NHITSMAX = $NHITSMAX"
    echo "  SAVEHISTS = $SAVEHISTS"

    threadLists=`updateThreadLists`
    if [ ! $? -eq 0 ]
    then
        echo "    ERROR in updateThreadLists!"
        exit 1
    fi
    Njobs=0
    for testlayer in $layers;
    do
        for (( iEvent=0; iEvent<nEvents; iEvent+=nEvtPerRun ))
        do
            ((Njobs++))
            iEntryStart=$iEvent
            iEntryStop=$((iEvent+nEvtPerRun-1))
            if (( iEntryStop>=nEvents ))
            then
                iEntryStop=$((nEvents-1))
            fi
            jobname="${runNo}.${currunname}.$iEntryStart-$iEntryStop.layer${testlayer}"
            echo "checking job \"$jobname\""
            file="root/t_${jobname}.log"
            if [ -e $file ]
            then # log file already exists??
                if tail -n 3 $file | grep -q "Good Events" # finished
                then
                    echo "  already finished!"
                    continue # no need to work on it
                else  # probably still under process
                    echo "  running by someone else!"
                    tail -n 1 $file
                    continue # don't know who is running it but ignore this job anyway
                fi
            else
                echo "  logfile \"$file\" doesn't exist, so generate a new job!"
            fi
            temprunname="${currunname}.$iEntryStart-$iEntryStop"
            tempconfig="$runNo $testlayer $prerunname $temprunname $nHitsGMax $t0shift0 $t0shift1 $tmin $tmax $geoSetup $sumCut $aaCut $iEntryStart $iEntryStop $workType $inputType $peakType"
            findVacentThread
            if [ $? -eq 1 ]
            then
                echo "    ERROR: cannot find a vacent thread in 10 hours!"
                exit 1
            elif [ $? -eq 2 ]
            then
                echo "    ERROR: cannot access hep_q in 10 minutes!"
                exit 1
            fi
            echo "    found new thread available $theConf"
            echo "$tempconfig" > $theConf # send the trigger info to the job
        done
    done
    echo "Starting iteration $iter, $Njobs jobs to be finished!"

    allfinished=false
    for (( i=0; i<3600; i++ )) #3600*10 sec = 10 hours running
    do
        if [ -e kill.$runNo.$runName ]
        then
            echo "Killed by user!"
            exit 0
        fi
        sleep 10
        echo -n "$i "
        finished=true
        thefile=""
        NjobsFinished=0
        for file in root/t_${runNo}.${currunname}.*.log
        do
            if [ ! -e $file ]
            then
                echo "WARNING: file \"$file\" doesn't exist"
                continue
            fi
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
            echo "#Iteration $iter finished"
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
    done
    if [ ! $allfinished ]
    then
        echo "ERROR! iteration $iter still not finished after 10 hours!"
        exit 1
    fi

    if [ -e kill.$runNo.$runName ]
    then
        echo "Killed by user!"
        exit 0
    fi

    cd root/
    combine $runNo $currunname $nEvtPerRun &
    pids="$!"
    wait $pids || { echo "there were errors in combining $runNo $currunname $ilayer" >&2; exit 1; }
    rm -f t_${runNo}.${currunname}.*-*.*
    cd ..

#   using layer0?
    if [ "$layers" == "0" ]
    then
        cd root
        for lid in $LAYERS
        do
            ln -s t_${runNo}.${currunname}.layer0.root t_${runNo}.${currunname}.layer${lid}.root
        done
        cd ..
    fi

#   upgrading wireposition?
    if [ ! $WPTYPE -eq 1 ] # ! 1 for not changing wiremap
    then
        cd info
        ln -s wire-position.${runNo}.${StartName}.root wire-position.${runNo}.${currunname}.root
        cd ..
    fi
    getOffset $runNo $prerunname $currunname $geoSetup $WPTYPE $scale $stepSize $minslz $maxslz $mininx $maxinx $maxchi2 -1 $wires
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
        getXT $runNo $prerunname $currunname $XTTYPE $geoSetup $SAVEHISTS $inputType $maxchi2 $DEFAULTLAYER $NHITSMAX
        lastxtfile=xt.${runNo}.${currunname}.root
    fi
done
