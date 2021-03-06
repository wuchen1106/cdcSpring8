#!/bin/bash

if [ $# -lt 7 ]
then
    echo $0 runNo runName thread_iStart nThreads istart istop [averageEtrack CHI2MAX NHITSSMIN NHITSMAX SLZMAX LAYER] 
    exit 0
fi

runNo="$1"
runName="$2"
thread_iStart=$3
nThreads=$4
thread_iStop=`echo "$thread_iStart+$nThreads-1"|bc`
IterStart=$5
IterEnd=$6
StartName=$7
isLast=false

layers="3 4 5 6" # layers to be reconstructed and [analyzed (in case of layers is not 0)]
if [ $# -gt 12 ]
then
    layers=${13}
fi

# for tracking
geoSetup=0 # 0 for general; 1 for finger
inputType=3 # 1 for MC; 0 for data; 2 for MC using X; 3 for MC using X and according to different test layers
workType=0 # 0, fr/l_0; 1, even/odd; -1, even/odd reversed; others, all layers
nHitsGMax=13
t0shift0=0
t0shift1=0
tmin=-10
tmax=800
sumCut=-10
aaCut=30
peakType=0 # 0, only the first peak over threshold; 1, all peaks over threshold; 2, even including shaddowed peaks

# for updateRes
averageEtrack=1 # use average Etrack
if [ $# -gt 7 ]
then
    averageEtrack=$8
fi

# for mc2fitinput_Chen
CHI2MAX=2
NHITSSMIN=7
NHITSMAX=0
SLZMAX=0.1
if [ $# -gt 8 ]
then
    CHI2MAX=$9
fi
if [ $# -gt 9 ]
then
    NHITSSMIN=${10}
fi
if [ $# -gt 10 ]
then
    NHITSMAX=${11}
fi
if [ $# -gt 11 ]
then
    SLZMAX=${12}
fi

threadName="job"

threadLists=""
threadlistfile=threadlist.$runName.$runNo

lastxtfile=""

echo "====================$0========================"
echo "You are going to start iteration ${IterStart}~${IterEnd} for run$runNo using threads \"$threadName\" $thread_iStart~$thread_iStop"
echo "StartName is $StartName, runName = $runName, layers for tracking \"$layers\""
echo "Is this the last iteration? $isLast"
echo "Tracking Parameters are:"
echo "        geoSetup = $geoSetup;  0 for general; 1 for finger"
echo "        inputType = $inputType;  1 for MC; 0 for data; 2 for MC using X; 3 for MC using X and according to different test layers"
echo "        workType = $workType;  0, fr/l_0; 1, even/odd; -1, even/odd reversed; others, all layers"
echo "        nHitsGMax = $nHitsGMax; "
echo "        t0shift0 = $t0shift0; "
echo "        t0shift1 = $t0shift1; "
echo "        tmin = $tmin; "
echo "        tmax = $tmax; "
echo "        sumCut = $sumCut; "
echo "        aaCut = $aaCut; "
echo "        peakType = $peakType;  0, only the first peak over threshold; 1, all peaks over threshold; 2, even including shaddowed peaks"
echo "updateRes Parameters are:"
echo "        averageEtrack? $averageEtrack; "
echo "mc2fitinput_Chen Parameters are:"
echo "        CHI2MAX = $CHI2MAX"
echo "        NHITSSMIN = $NHITSSMIN"
echo "        NHITSMAX = $NHITSMAX"
echo "        SLZMAX = $SLZMAX"

read -p 'You are going to do the job above, is that right? (Y/n):'
if [ ! "$REPLY" = 'Y' ] && [ ! "$REPLY" = 'y' ] && [ ! "$REPLY" = '' ]; then
    exit 0
fi

# check files needed
echo "Will use following files as input"
#if [ -e $CDCS8WORKING_DIR/root/h_$runNo.root ]
#then
#    ls -ltr $CDCS8WORKING_DIR/root/h_$runNo.root
#else
#    echo "$CDCS8WORKING_DIR/root/h_$runNo.root doesn't exist!"
#    exit 1
#fi
if [ -e $CDCS8WORKING_DIR/info/wire-position.$runNo.$StartName.root ]
then
    ls -ltr $CDCS8WORKING_DIR/info/wire-position.$runNo.$StartName.root
else
    echo "$CDCS8WORKING_DIR/info/wire-position.$runNo.$StartName.root doesn't exist!"
    exit 1
fi
if [ -e $CDCS8WORKING_DIR/info/xt.$runNo.$StartName.root ]
then
    ls -ltr $CDCS8WORKING_DIR/info/xt.$runNo.$StartName.root
else
    echo "$CDCS8WORKING_DIR/info/xt.$runNo.$StartName.root doesn't exist!"
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

im1=$((IterStart-1))
prerunname="${runName}.i${im1}"
for testLayer in $layers
do
    updateRes $runNo $StartName $prerunname 1 $averageEtrack $testLayer
done
for (( iter=IterStart; iter<=IterEnd; iter++ ))
do
    if [ -e kill.$runNo.$runName ]
    then
        echo "Killed by user!"
        exit 0
    fi

    im1=$((iter-1))
    prerunname="${runName}.i${im1}"
    currunname="${runName}.i${iter}"

    echo "#Iteration $iter started"

    threadLists=`updateThreadLists` # update the thread List for this iteration
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
    if [ "$nThreads" -eq 0 ]
    then
        echo "Cannot find any available thread!"
        exit 1
    fi

    for testLayer in $layers;
    do
        nEvents=`GetEntries $CDCS8WORKING_DIR/root/h_$runNo.root`
        nEvtPR=$((nEvents/5))
        pids=""
        files=""
        for (( iev=0; iev<nEvents; iev+=nEvtPR ))
        do
            istart=$iev
            istop=$((iev+nEvtPR-1))
            if (( istop>=nEvents ))
            then
                istop=$((nEvents-1))
            fi
            mc2fitinput_Chen root/ana_$runNo.$StartName.layer${testLayer}.root info/res.$runNo.layer${testLayer}.$prerunname.root root/h_$runNo.${istart}-${istop}.MC.root 0 0 0 8 1 1 $istart $istop $CHI2MAX $NHITSSMIN $NHITSMAX $SLZMAX &
            pids="$pids $!"
            files="$files root/h_$runNo.${istart}-${istop}.MC.root"
        done
        wait $pids || { echo "there were errors in combining $runNo $currunname $testLayer" >&2; exit 1; }
        if [ -e root/h_$runNo.layer${testLayer}.MC.root ]
        then
            rm root/h_$runNo.layer${testLayer}.MC.root
        fi
        hadd root/h_$runNo.layer${testLayer}.MC.root $files;
        rm $files;
        nEvents=`GetEntries $CDCS8WORKING_DIR/root/h_$runNo.layer${testLayer}.MC.root`
        nEvtPerRun=`echo "$nEvents/($nThreads)+1" | bc`
        Njobs=0 # to count number of jobs to be finished
        for (( iEvent=0; iEvent<nEvents; iEvent+=nEvtPerRun ))
        do
            ((Njobs++))
            iEntryStart=$iEvent
            iEntryStop=$((iEvent+nEvtPerRun-1))
            if (( iEntryStop>=nEvents ))
            then
                iEntryStop=$((nEvents-1))
            fi
            jobname="${runNo}.${currunname}.$iEntryStart-$iEntryStop.layer${testLayer}"
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
            tempconfig="    $runNo $testLayer $StartName $temprunname $nHitsGMax $t0shift0 $t0shift1 $tmin $tmax $geoSetup $sumCut $aaCut $iEntryStart $iEntryStop $workType $inputType $peakType"
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
        echo "Starting iteration $iter for layer $testLayer, $Njobs jobs to be finished!"

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
        combine $runNo $currunname $nEvtPerRun h_$runNo.layer${testLayer}.MC.root $testLayer &
        pids="$!"
        wait $pids || { echo "there were errors in combining $runNo $currunname" >&2; exit 1; }
        rm -f t_${runNo}.${currunname}.*-*.*
        cd ..
    done

#   updating
    for testLayer in $layers
    do
        updateRes $runNo $StartName $currunname 0 $averageEtrack $testLayer
    done
    mv root/h_$runNo.layer${testLayer}.MC.root root/h_$runNo.layer${testLayer}.$prerunname.MC.root
done
