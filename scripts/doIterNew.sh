#!/bin/bash
THISCMD=$(basename $0)
theDate=`date`
echo $theDate $THISCMD $@ >> cmdlog

# about thread
threadName="job"

# about this job
StartName="Garfield"
layers="4" # layers to be reconstructed and [analyzed (in case of layers is not 0)]
wires="" # wires to be calibrated (position) FIXME NOT SUPPORTED YET
isLast=false
PROGRAMMED=false # by default don't load the programmed iteration setting
CONFIGTABLEDEFAULT="$CDCS8WORKING_DIR/Para/new.dat"
CONFIGTABLE=""

# for getOffset
UpdateWireMap=false # FIXME NOT SUPPORTED YET

usages() {
cat << EOF
$THISCMD:  do the iteration!
Syntax:
    $THISCMD [options] runname
    [options]
    -h     display this help and exit
    -R [R] the run number
    -T [T] the staring thread index
    -N [N] the number of threads to ask for
    -I [I] start from iteration I
    -J [J] stop at iteration J
    -L     (false) Take the last iteration as the final step and do the default complete checkings
    -W     (false) update the wire position map
    -P     (false) use programmed iteration parameters
    -C [c] Add configure file as C (default one \"$CONFIGTABLE\" is always loaded first)
    -S [S] set the start name ($StartName)
    -l [l1 (l2 ...)] Do the training of the given layers ($layers)
    -w [w1 (w2 ...)] Do the calibration of the given wires ($wires) in the ubove given layers ($layers)

Report bugs to <wuchen@ihep.ac.cn>.
EOF
}

while getopts ':hR:T:N:I:J:LW:S:Pl:w:c:' optname
do
    case "$optname" in
    'h')
        usages
        exit 0
        ;;
    'R')
        runNo="$OPTARG"
        ;;
    'T')
        thread_iStart="$OPTARG"
        ;;
    'N')
        thread_iStop=`echo "$thread_iStart+$OPTARG-1"|bc`
        ;;
    'I')
        IterStart="$OPTARG"
        ;;
    'J')
        IterEnd="$OPTARG"
        ;;
    'L')
        isLast=true;
        ;;
    'W')
        UpdateWireMap=true;
        ;;
    'S')
        StartName="$OPTARG";
        ;;
    'P')
        PROGRAMMED=true;
        ;;
    'l')
        layers="$OPTARG";
        ;;
    'w')
        wires="$OPTARG";
        ;;
    'c')
        CONFIGTABLE="$CONFIGTABLE -C $OPTARG";
        ;;
    '?')
        echo "Unknown option $OPTARG"
        echo "Try \"$THISCMD -h\" for more infomation"
        exit -1
        ;;
    ':')
        echo "Need argument value for option $OPTARG"
        echo "Try \"$THISCMD -h\" for more infomation"
        exit -1
        ;;
    *)
        # Should not occur
        echo 'Unknown error while processing options'
        exit -1
        ;;
    esac
done

case "$(($#+1-$OPTIND))" in
0)
    echo You have to intput a run name!
    usages
    exit -1
    ;;
1)
    runName=
    runName="${@:$OPTIND:1}"
    ;;
*)
    echo 'Too many jobs. Only 1 job permitted'
    usages
    exit -1
    ;;
esac

echo "You are going to start iteration ${IterStart}~${IterEnd} for run$runNo using threads \"$threadName\" $thread_iStart~$thread_iStop"
echo "StartName is $StartName, runName = $runName, layers for tracking \"$layers\", wires to calibrate \"$wires\""
echo "Use programmed parameter setting for iteration? $PROGRAMMED"
echo "Is this the last iteration? $isLast"
echo "Configure file for general: $CONFIGTABLE"
echo "getOffset Parameters are:"
echo "        UpdateWireMap = $UpdateWireMap;"

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
if [ -e $CDCS8WORKING_DIR/root/hits/h_$runNo.root ]
then
    ls -ltr $CDCS8WORKING_DIR/root/hits/h_$runNo.root
else
    echo "$CDCS8WORKING_DIR/root/hits/h_$runNo.root doesn't exist!"
    exit 1
fi
if [ -e $CDCS8WORKING_DIR/info/wire-position.$runNo.$prerunname.root ]
then
    ls -ltr $CDCS8WORKING_DIR/info/wire-position.$runNo.$prerunname.root
else
    echo "$CDCS8WORKING_DIR/info/wire-position.$runNo.$prerunname.root doesn't exist!"
    if [ -e $CDCS8WORKING_DIR/Input/wire-position.root ]
    then
        echo "Will use $CDCS8WORKING_DIR/Input/wire-position.root by default"
    else
        echo "Even $CDCS8WORKING_DIR/Input/wire-position.root doesn't exist!"
        exit 1
    fi
fi
if [ -e $CDCS8WORKING_DIR/info/xt.$runNo.$prerunname.root ]
then
    ls -ltr $CDCS8WORKING_DIR/info/xt.$runNo.$prerunname.root
else
    echo "$CDCS8WORKING_DIR/info/xt.$runNo.$prerunname.root doesn't exist!"
    echo "Will use default garfield xt"
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

prev_ithread=$thread_iStart
prev_occupied=false
threadLists=""
threadlistfile=threadlist.$runName.$runNo
lastxtfile=""
declare -a sourcefiles
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

    if $PROGRAMMED
    then
        if [ $iter -gt 20 ]
        then
            layers="1 2 3 4 5 6 7 8"
            wires="3 4 5 6"
            UpdateWireMap=true
        elif [ $iter -gt 5 ]
        then
            layers="1 2 3 4 5 6 7 8"
            wires="3 4 5 6"
            UpdateWireMap=true
        else
            layers=4
            wires=""
            UpdateWireMap=false
        fi
    fi

    if [ $iter -eq $IterEnd ] && $isLast
    then
        layers="1 2 3 4 5 6 7 8"
#        layers="0"
    fi

#   tune arguments
    arg_configure=""
    if [ ! -z "$CONFIGTABLE" ]
    then
        arg_configure="$CONFIGTABLE"
    fi

    echo "#Iteration $iter started"
    echo "  layers = \"$layers\""
    echo "  wires = \"$wires\""
    echo "getOffset Parameters:"
    echo "  UpdateWireMap = $UpdateWireMap"
    echo "arguments:"
    echo "  arg_configure = $arg_configure"

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
    nEvents=`GetEntries $CDCS8WORKING_DIR/root/hits/h_$runNo.root`
    nEvtPerRun=`echo "$nEvents/($nThreads)+1" | bc`

    Njobs=0
    for testlayer in $layers;
    do
        sourcefiles[testlayer]=""
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

            file="root/tracks/t_${jobname}.log"
            temprunname=`printf "${currunname}.%07d-%07d" $iEntryStart $iEntryStop`
            logtemp="$CDCS8WORKING_DIR/root/tracks/t_${runNo}.${temprunname}.layer${testlayer}.log"
            errtemp="$CDCS8WORKING_DIR/root/tracks/t_${runNo}.${temprunname}.layer${testlayer}.err"
            sourcefiles[testlayer]="${sourcefiles[testlayer]} t_${runNo}.${temprunname}.layer${testlayer}.root"
            tempconfig="Tracking -C $CONFIGTABLEDEFAULT $arg_configure -R $runNo -L $testlayer -B $iEntryStart -E $iEntryStop $prerunname $temprunname > $logtemp 2> $errtemp"
            echo $tempconfig

            rootfile="root/tracks/t_${jobname}.root"
            rootsize=0
            if [ -e $rootfile ]
            then
                rootsize=`du -s $rootfile | gawk '{print $1;}'`
            else
                echo "  ROOT file \"$rootfile\" not found yet!"
            fi

            if [ -e "$file" ]
            then # log file already exists??
                if tail -n 3 $file | grep -q "Good Events" # finished
                then
                    if [ $rootsize -eq 0 ]
                    then
                        echo "  ROOT file is empty although log file shows finished status!"
                    else
                        echo "  already finished!"
                        continue # no need to work on it
                    fi
                else  # probably still under process
                    if [ $rootsize -eq 0 ]
                    then
                        echo "  ROOT file is empty!"
                    else
                        echo "  running by someone else!"
                        tail -n 1 $file
                        continue # don't know who is running it but ignore this job anyway
                    fi
                fi
            else
                echo "  logfile \"$file\" doesn't exist, so generate a new job!"
            fi
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
        for file in root/tracks/t_${runNo}.${currunname}.*.log
        do
            if [ ! -e $file ]
            then
                echo "WARNING: file \"$file\" doesn't exist"
                continue
            fi
            rootfile=`echo $file | sed 's/.log/.root/g'`
            rootsize=0
            if [ -e $rootfile ]
            then
                rootsize=`du -s $rootfile | gawk '{print $1;}'`
            else
                echo "  ROOT file not found yet!"
            fi

            if ! tail -n 3 $file | grep -q "Good Events"
            then
                finished=false
                thefile=$file
            else
                if [ $rootsize -eq 0 ]
                then
                    echo "  ROOT file is empty"
                    finished=false
                    thefile=$file

                    jobname=`echo $file | sed 's/.*\<t_\(.*\)\.log/\1/g'`
                    echo "resubmitting job \"$jobname\""

                    temprunname=`echo $jobname | sed 's/\w*\.\(.*\)\.layer\w*/\1/g'`
                    testlayer=`echo $jobname | sed 's/.*layer\(\w*\)\..*/\1/'`
                    iEntryStart=`echo $jobname | sed 's/.*\.\(\w*\)-\(\w*\)\..*/\1/g'`
                    iEntryStop=`echo $jobname | sed 's/.*\.\(\w*\)-\(\w*\)\..*/\2/g'`
                    logtemp="$CDCS8WORKING_DIR/root/tracks/t_${runNo}.${temprunname}.layer${testlayer}.log"
                    errtemp="$CDCS8WORKING_DIR/root/tracks/t_${runNo}.${temprunname}.layer${testlayer}.err"
                    tempconfig="Tracking -C $CONFIGTABLEDEFAULT $arg_configure -R $runNo -L $testlayer -B $iEntryStart -E $iEntryStop $prerunname $temprunname > $logtemp 2> $errtemp"
                    echo $tempconfig

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
                else
                    ((NjobsFinished++))
                fi
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

    cd root/tracks
    #combine $runNo $currunname $nEvtPerRun &
    for testlayer in $layers;
    do
        hadd -f t_${runNo}.${currunname}.layer$testlayer.root ${sourcefiles[testlayer]} &
        pids="$pids $!"
    done
    wait $pids || { echo "there were errors in combining $runNo $currunname $ilayer" >&2; exit 1; }
    rm -f t_${runNo}.${currunname}.*-*.*
    cd ../..

#   using layer0?
    if [ "$layers" == "0" ]
    then
        cd root/tracks
        for lid in $layers
        do
            ln -s t_${runNo}.${currunname}.layer0.root t_${runNo}.${currunname}.layer${lid}.root
        done
        cd ../..
    fi

    GetXT -C $CONFIGTABLEDEFAULT $arg_configure -R $runNo $prerunname $currunname
done
