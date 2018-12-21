#!/bin/bash

if [ $# -lt 6 ]
then
    echo $0 runNo runName thread_iStart nThreads iStart iStop [aaCut isLast NHITSMAX] 
    exit 0
fi

runNo="$1"
runName="$2"
thread_iStart=$3
nThreads=$4
thread_iStop=`echo "$thread_iStart+$nThreads-1"|bc`
IterStart=$5
IterEnd=$6
isLast=false
if [ $# -gt 7 ]
then
    isLast=$8
fi

StartName="Garfield"
layers="4" # layers to be reconstructed and [analyzed (in case of layers is not 0)]
LAYERS="4" # layers to be analyzed (in case of layers is 0)
wires="" # wires to be calibrated (position)

# for tracking
geoSetup=0 # 0 for general; 1 for finger
inputType=0 # 1 for MC; 0 for data
workTypeini=0 # 0, fr/l_0; 1, even/odd; -1, even/odd reversed; others, all layers
nHitsGMaxini=15
t0shift0=0
t0shift1=0
tmin=-10
tmax=800
sumCut=-10
aaCut=0
if [ $# -gt 6 ]
then
    aaCut=$7
fi
peakType=0 # 0, only the first peak over threshold; 1, all peaks over threshold; 2, even including shaddowed peaks

# for getOffset
WPTYPE=0 # 0 for changing wiremap; 1 for not changing it;
stepSize=0 # maximum step size for each movement in wire position calibration; 0 means no limit
minslz=0 # min slz cut for mean slz value in each sample of events in wiremap calibration. minslz==maxslz==0 means no cut
maxslz=0 # min slz cut for mean slz value in each sample of events in wiremap calibration. minslz==maxslz==0 means no cut
mininx=0 # min inx cut for mean inx value in each sample of events in wiremap calibration. mininx==maxinx==0 means no cut
maxinx=0 # min inx cut for mean inx value in each sample of events in wiremap calibration. mininx==maxinx==0 means no cut
maxchi2=2
scale=1 # move scale*offset on wiremap for fitting in the next round

# for ana
UPDATEXT=1
DEFAULTLAYER=4 # use this layer to generate fl(r)_0 and so on
if [ $# -gt 9 ]
then
    XTTYPE=$10
fi
NHITSMAXini=30
if [ $# -gt 8 ]
then
    NHITSMAXini=$9
fi
SAVEHISTS=0

threadName="job"

threadLists=""
threadlistfile=threadlist.$runName.$runNo

lastxtfile=""

echo "You are going to start iteration ${IterStart}~${IterEnd} for run$runNo using threads \"$threadName\" $thread_iStart~$thread_iStop"
echo "StartName is $StartName, runName = $runName, layers for tracking \"$layers\", layers for ana \"$LAYERS\", wires to calibrate \"$wires\""
echo "Is this the last iteration? $isLast"
echo "Tracking Parameters are:"
echo "        geoSetup = $geoSetup;  0 for general; 1 for finger"
echo "        inputType = $inputType;  1 for MC; 0 for data"
echo "        workTypeini = $workTypeini;  0, fr/l_0; 1, even/odd; -1, even/odd reversed; others, all layers"
echo "        nHitsGMaxini = $nHitsGMaxini; "
echo "        t0shift0 = $t0shift0; "
echo "        t0shift1 = $t0shift1; "
echo "        tmin = $tmin; "
echo "        tmax = $tmax; "
echo "        sumCut = $sumCut; "
echo "        aaCut = $aaCut; "
echo "        peakType = $peakType;  0, only the first peak over threshold; 1, all peaks over threshold; 2, even including shaddowed peaks"
echo "getOffset Parameters are:"
echo "        WPTYPE = $WPTYPE;  0 for changing wiremap; 1 for not changing it;"
echo "        stepSize = $stepSize;  maximum step size for each movement in wire position calibration; 0 means no limit"
echo "        minslz = $minslz;  min slz cut for mean slz value in each sample of events in wiremap calibration. minslz==maxslz==0 means no cut"
echo "        maxslz = $maxslz;  min slz cut for mean slz value in each sample of events in wiremap calibration. minslz==maxslz==0 means no cut"
echo "        mininx = $mininx;  min inx cut for mean inx value in each sample of events in wiremap calibration. mininx==maxinx==0 means no cut"
echo "        maxinx = $maxinx;  min inx cut for mean inx value in each sample of events in wiremap calibration. mininx==maxinx==0 means no cut"
echo "        maxchi2 = $maxchi2; "
echo "        scale = $scale;  move scale*offset on wiremap for fitting in the next round"
echo "ana Parameters are:"
echo "        UPDATEXT = $UPDATEXT; "
echo "        DEFAULTLAYER = $DEFAULTLAYER;  use this layer to generate fl(r)_0 and so on"
echo "        XTTYPE = $XTTYPE;"
echo "        NHITSMAXini = $NHITSMAXini; "
echo "        SAVEHISTS = $SAVEHISTS; "

read -p 'You are going to do the job above, is that right? (Y/n):'
if [ ! "$REPLY" = 'Y' ] && [ ! "$REPLY" = 'y' ] && [ ! "$REPLY" = '' ]; then
    exit 0
fi

# check files needed
echo "Will use following files as input"
if [ $IterStart == 1 ]
then
    prerunname="${StartName}"
else
    im1=$((IterStart-1))
    prerunname="${runName}.i${im1}"
fi
if [ -e $CDCS8WORKING_DIR/root/h_$runNo.root ]
then
    ls -ltr $CDCS8WORKING_DIR/root/h_$runNo.root
else
    echo "$CDCS8WORKING_DIR/root/h_$runNo.root doesn't exist!"
    exit 1
fi
if [ -e $CDCS8WORKING_DIR/info/wire-position.$runNo.$prerunname.root ]
then
    ls -ltr $CDCS8WORKING_DIR/info/wire-position.$runNo.$prerunname.root
else
    echo "$CDCS8WORKING_DIR/info/wire-position.$runNo.$prerunname.root doesn't exist!"
    exit 1
fi
if [ -e $CDCS8WORKING_DIR/info/xt.$runNo.$prerunname.root ]
then
    ls -ltr $CDCS8WORKING_DIR/info/xt.$runNo.$prerunname.root
else
    echo "$CDCS8WORKING_DIR/info/xt.$runNo.$prerunname.root doesn't exist!"
    exit 1
fi

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
        sleep 2
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
    if [ -e kill.$runNo.$runName ]
    then
        echo "Killed by user!"
        exit 0
    fi

    im1=$((iter-1))
    if [ $iter == 1 ]
    then
        currunname="${runName}.i${iter}"
        prerunname="${StartName}"
    else
        currunname="${runName}.i${iter}"
        prerunname="${runName}.i${im1}"
    fi

    if [ $iter -gt 1 ]
    then
        workType=$workTypeini
    else
        workType=0
    fi

    if [ $iter -eq $IterEnd ] && $isLast
    then
        nHitsGMax=$nHitsGMaxini
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
    nThreads=0 # count number of threads available for this iteration (sometimes some threads are stopped in the middle)
    for (( ithread=thread_iStart; ithread<=thread_iStop; ithread++ ))
    do
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
                ((nThreads++))
            fi
        fi
    done
    nEvents=`GetEntries $CDCS8WORKING_DIR/root/h_$runNo.root`
    nEvtPerRun=`echo "$nEvents/($nThreads)+1" | bc`

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
            logtemp="$CDCS8WORKING_DIR/root/t_${runNo}.${temprunname}.layer${testlayer}.log"
            errtemp="$CDCS8WORKING_DIR/root/t_${runNo}.${temprunname}.layer${testlayer}.err"
            tempconfig="tracking -R $runNo -L $testlayer -n $nHitsGMax -x $t0shift0 -y $t0shift1 -l $tmin -u $tmax -g $geoSetup -s $sumCut -a $aaCut -B $iEntryStart -E $iEntryStop -w $workType -i $inputType -p $peakType $prerunname $temprunname > $logtemp 2> $errtemp"
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
        ana -R $runNo -x $XTTYPE -g $geoSetup -H $SAVEHISTS -i $inputType -c $maxchi2 -L $DEFAULTLAYER -n $NHITSMAX $prerunname $currunname
        lastxtfile=xt.${runNo}.${currunname}.root
    fi
done
