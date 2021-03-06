#!/bin/bash
THISCMD=$(basename $0)
theDate=`date`
echo $theDate $THISCMD $@ >> cmdlog

# about thread
threadName="job"

# about this job
StartName="Garfield"
layers="4" # layers to be reconstructed and [analyzed (in case of layers is not 0)]
wires="" # wires to be calibrated (position)
isLast=false
PROGRAMMED=false # by default don't load the programmed iteration setting
CONFIGTABLEDEFAULT="$CDCS8WORKING_DIR/Para/default.dat"
CONFIGTABLE_wireposition="$CDCS8WORKING_DIR/Para/wirepos.temp.dat"
CONFIGTABLE=""

# for tracking
geoSetup=0 # 0 for general; 1 for finger; 2 for tilted run
inputType=0 # 1 for MC; 0 for data
workTypeini=0 # 0, fr/l_0; 1, even/odd; -1, even/odd reversed; others, all layers
nHitsGMax=15
sumCut=-10
aaCut=0
peakType=0 # 0, only the first peak over threshold; 1, all peaks over threshold; 2, even including shaddowed peaks
blindLayer=-1 # don't use this layer for tracking (if it's not -1)

t0shift0=0
t0shift1=0
tmin=-10
tmax=800

# for getOffset
UpdateWireMap=false
OneWirePerIter=false
stepSize=0 # maximum step size for each movement in wire position calibration; 0 means no limit
minslz=0 # min slz cut for mean slz value in each sample of events in wiremap calibration. minslz==maxslz==0 means no cut
maxslz=0 # min slz cut for mean slz value in each sample of events in wiremap calibration. minslz==maxslz==0 means no cut
mininx=0 # min inx cut for mean inx value in each sample of events in wiremap calibration. mininx==maxinx==0 means no cut
maxinx=0 # min inx cut for mean inx value in each sample of events in wiremap calibration. mininx==maxinx==0 means no cut
calib_maxslz=0.025 # maximum slz for each event considered as offset estimation
calib_maxFD=7.5 # maximum FD for each event considered as offset estimation; also defining the region for a golden hit
calib_minFD=0.5 # minimum FD for each event considered as offset estimation; also defining the region for a golden hit
calib_allgolden=false # requiring all golden hits in getting offset (since the XT might not be good in "bad" regions that are not aligned well)
maxchi2=2
scale=0.5 # move scale*offset on wiremap for fitting in the next round

# for ana
DONTUPDATEXT=false
DEFAULTLAYER=4 # use this layer to generate fl(r)_0 and so on
XTTYPE=055
NHITSMAX=30
SAVEHISTS=0
tmaxSet=800
allgolden=false # resubmitting all golden hits in getting new xt

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
    -H [h] save histograms at this level ($SAVEHISTS)
    -L     (false) Take the last iteration as the final step and do the default complete checkings
    -W     (false) update the wire position map
    -U     (false) keep the xt curves unchanged.
    -O     (false) update one wire position per iteration step (the one with the largest offset)
    -P     (false) use programmed iteration parameters
    -G     (false) require all golden hits for getting XT 
    -A     (false) require all golden hits for getting wire offset
    -S [S] set the start name ($StartName)
    -D [D] set this layer ($DEFAULTLAYER) as the default layer to save in XT file as fr/l_0
    -l [l1 (l2 ...)] Do the training of the given layers ($layers)
    -w [w1 (w2 ...)] Do the calibration of the given wires ($wires) in the ubove given layers ($layers)
    -c [c] Add configure file as C (default one \"$CONFIGTABLE\" is always loaded first)
    -g [g] geometry setup ($geoSetup). 0 ordinary scintillator; 1 finger scintillator
    -t [t] work type ($workTypeini) for tracking. 0, fr/l_0; 1, even/odd; -1, even/odd reversed; others, all layers
    -a [a] aa cut ($aaCut)
    -b [b] Blind this layer from tracking ($blindLayer)
    -s [s] sum cut ($sumCut)
    -p [p] set peak type ($peakType) for tracking. 0, only the first peak over threshold; 1, all peaks over threshold; 2, even including shaddowed peaks
    -x [XYZ] xt type. XYZ ($XTTYPE) means polX for center, polY for middle and polZ for tail. If X is 0 then let middle function fit the center region.
    -n [n] maximum number ($NHITSMAX) of hits to be used in ana
    -m [m] maximum number ($nHitsGMax) of good hits to be used in tracking
    -u [u] maximum drift time ($tmax) to be considered in tracking
    -d [d] minimum drift time ($tmin) to be considered in tracking
    -o [o] maximum drift time ($tmaxSet) to be considered in getting xt

Report bugs to <wuchen@ihep.ac.cn>.
EOF
}

while getopts ':hR:T:N:I:J:LH:W:US:D:PAGl:w:c:g:t:a:b:s:p:x:n:m:u:d:o:' optname
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
    'H')
        SAVEHISTS="$OPTARG";
        ;;
    'W')
        UpdateWireMap=true;
        ;;
    'U')
        DONTUPDATEXT=true;
        ;;
    'S')
        StartName="$OPTARG";
        ;;
    'D')
        DEFAULTLAYER="$OPTARG";
        layers=$DEFAULTLAYER;
        ;;
    'O')
        OneWirePerIter=true;
        ;;
    'P')
        PROGRAMMED=true;
        ;;
    'A')
        calib_allgolden=true;
        ;;
    'G')
        allgolden=true;
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
    'g')
        geoSetup="$OPTARG";
        ;;
    't')
        workTypeini="$OPTARG";
        ;;
    'i')
        inputType="$OPTARG";
        ;;
    'a')
        aaCut="$OPTARG"
        ;;
    'b')
        blindLayer="$OPTARG"
        ;;
    's')
        sumCut="$OPTARG"
        ;;
    'p')
        peakType="$OPTARG"
        ;;
    'x')
        XTTYPE="$OPTARG"
        ;;
    'n')
        NHITSMAX="$OPTARG"
        ;;
    'm')
        nHitsGMax="$OPTARG"
        ;;
    'u')
        tmax="$OPTARG"
        ;;
    'd')
        tmin="$OPTARG"
        ;;
    'o')
        tmaxSet="$OPTARG"
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
echo "Tracking Parameters are:"
echo "        geoSetup = $geoSetup;  0 for general; 1 for finger"
echo "        inputType = $inputType;  1 for MC; 0 for data"
echo "        workTypeini = $workTypeini;  0, fr/l_0; 1, even/odd; -1, even/odd reversed; others, all layers"
echo "        nHitsGMax = $nHitsGMax; "
echo "        t0shift0 = $t0shift0; "
echo "        t0shift1 = $t0shift1; "
echo "        tmin = $tmin; "
echo "        tmax = $tmax; "
echo "        sumCut = $sumCut; "
echo "        aaCut = $aaCut; "
echo "        Blind layer: $blindLayer; "
echo "        peakType = $peakType;  0, only the first peak over threshold; 1, all peaks over threshold; 2, even including shaddowed peaks"
echo "getOffset Parameters are:"
echo "        OneWirePerIter = $OneWirePerIter;"
echo "        UpdateWireMap = $UpdateWireMap;"
echo "        stepSize = $stepSize;  maximum step size for each movement in wire position calibration; 0 means no limit"
echo "        scale = $scale;  move scale*offset on wiremap for fitting in the next round"
echo "        minslz = $minslz;  min slz cut for mean slz value in each sample of events in wiremap calibration. minslz==maxslz==0 means no cut"
echo "        maxslz = $maxslz;  min slz cut for mean slz value in each sample of events in wiremap calibration. minslz==maxslz==0 means no cut"
echo "        mininx = $mininx;  min inx cut for mean inx value in each sample of events in wiremap calibration. mininx==maxinx==0 means no cut"
echo "        maxinx = $maxinx;  min inx cut for mean inx value in each sample of events in wiremap calibration. mininx==maxinx==0 means no cut"
echo "        maxchi2 = $maxchi2; "
echo "        calib_maxslz = $calib_maxslz;   maximum slz for each event considered as offset estimation"
echo "        calib_maxFD = $calib_maxFD"
echo "        calib_minFD = $calib_minFD"
echo "        calib_allgolden = $calib_allgolden"
echo "ana Parameters are:"
echo "        DONTUPDATEXT = $DONTUPDATEXT; "
echo "        DEFAULTLAYER = $DEFAULTLAYER;  use this layer to generate fl(r)_0 and so on"
echo "        XTTYPE = $XTTYPE;"
echo "        NHITSMAX = $NHITSMAX; "
echo "        SAVEHISTS = $SAVEHISTS; "
echo "        tmaxSet = $tmaxSet; "
echo "        allgolden = $allgolden"

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

    if [ $iter -gt 1 ]
    then
        workType=$workTypeini
    else
        workType=0
    fi

    if $PROGRAMMED
    then
        if [ $iter -gt 20 ]
        then
            layers="1 2 3 4 5 6 7 8"
            wires="3 4 5 6"
            UpdateWireMap=true
            OneWirePerIter=false
        elif [ $iter -gt 5 ]
        then
            layers="1 2 3 4 5 6 7 8"
            wires="3 4 5 6"
            UpdateWireMap=true
            OneWirePerIter=true
        else
            layers=$DEFAULTLAYER
            wires=""
            UpdateWireMap=false
        fi
    fi

    if [ $iter -eq $IterEnd ] && $isLast
    then
        layers="1 2 3 4 5 6 7 8"
#        layers="0"
    fi

#   update the configure file
    cat << EOF > $CONFIGTABLE_wireposition
    < ana.scale = $scale >
    < ana.stepSize = $stepSize >
    < ana.minDeltaSlz = $minslz >
    < ana.maxDeltaSlz = $maxslz >
    < ana.minDeltaInx = $mininx >
    < ana.maxDeltaInx = $maxinx >
    < ana.calib_maxslz = $calib_maxslz >
    < ana.calib_maxFD = $calib_maxFD >
    < ana.calib_minFD = $calib_minFD >
    < ana.calib_allGolden = $calib_allgolden >
EOF

#   tune arguments
    arg_configure=""
    if [ ! -z "$CONFIGTABLE" ]
    then
        arg_configure="$CONFIGTABLE"
    fi
    arg_wiremap=""
    if $UpdateWireMap
    then
        if $OneWirePerIter
        then
            arg_wiremap="-W -O -C $CONFIGTABLE_wireposition"
        else
            arg_wiremap="-W -C $CONFIGTABLE_wireposition"
        fi
    else
        arg_wiremap="-C $CONFIGTABLE_wireposition"
    fi
    arg_xtfile=""
    if $DONTUPDATEXT
    then
        if [ -z "$lastxtfile" ]
        then
            lastxtfile=xt.${runNo}.${prerunname}.root
        fi
        arg_xtfile="-X $lastxtfile"
    else
        lastxtfile=xt.${runNo}.${currunname}.root
    fi
    arg_draw=""
    if [ $iter -eq $IterEnd ] && $isLast
    then
        arg_draw="-Z"
    fi
    arg_blindLayer=""
    if [ ! $blindLayer -eq -1 ]
    then
        arg_blindLayer="-b $blindLayer"
    fi
    arg_xtopt=""
    if $allgolden
    then
        arg_xtopt="-G"
    fi

    echo "#Iteration $iter started"
    echo "  layers = \"$layers\""
    echo "  wires = \"$wires\""
    echo "Tracking Parameters are:"
    echo "  geoSetup = $geoSetup"
    echo "  inputType = $inputType"
    echo "  workTypeini = $workTypeini"
    echo "  workType = $workType"
    echo "  nHitsGMax = $nHitsGMax"
    echo "  t0shift0 = $t0shift0"
    echo "  t0shift1 = $t0shift1"
    echo "  tmin = $tmin"
    echo "  tmax = $tmax"
    echo "  sumCut = $sumCut"
    echo "  aaCut = $aaCut"
    echo "  Blind layer: $blindLayer"
    echo "  peakType = $peakType"
    echo "getOffset Parameters:"
    echo "  OneWirePerIter = $OneWirePerIter"
    echo "  UpdateWireMap = $UpdateWireMap"
    echo "  scale = $scale"
    echo "  stepSize = $stepSize"
    echo "  minslz = $minslz"
    echo "  maxslz = $maxslz"
    echo "  mininx = $mininx"
    echo "  maxinx = $maxinx"
    echo "  maxchi2 = $maxchi2"
    echo "  calib_maxslz = $calib_maxslz"
    echo "  calib_maxFD = $calib_maxFD"
    echo "  calib_minFD = $calib_minFD"
    echo "  calib_allgolden = $calib_allgolden"
    echo "ana Parameters:"
    echo "  DONTUPDATEXT = $DONTUPDATEXT"
    echo "  DEFAULTLAYER = $DEFAULTLAYER"
    echo "  XTTYPE = $XTTYPE"
    echo "  NHITSMAX = $NHITSMAX"
    echo "  SAVEHISTS = $SAVEHISTS"
    echo "  tmaxSet = $tmaxSet"
    echo "  allgolden = $allgolden"
    echo "arguments:"
    echo "  arg_configure = $arg_configure"
    echo "  arg_xtfile = $arg_xtfile"
    echo "  arg_wiremap = $arg_wiremap"
    echo "  arg_blindLayer = $arg_blindLayer"

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
            temprunname="${currunname}.$iEntryStart-$iEntryStop"
            logtemp="$CDCS8WORKING_DIR/root/tracks/t_${runNo}.${temprunname}.layer${testlayer}.log"
            errtemp="$CDCS8WORKING_DIR/root/tracks/t_${runNo}.${temprunname}.layer${testlayer}.err"
            sourcefiles[testlayer]="${sourcefiles[testlayer]} t_${runNo}.${temprunname}.layer${testlayer}.root"
            tempconfig="tracking -C $CONFIGTABLEDEFAULT $arg_configure $arg_blindLayer -R $runNo -L $testlayer -n $nHitsGMax -x $t0shift0 -y $t0shift1 -l $tmin -u $tmax -g $geoSetup -s $sumCut -a $aaCut -B $iEntryStart -E $iEntryStop -w $workType -i $inputType -p $peakType $prerunname $temprunname > $logtemp 2> $errtemp"
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
                    tempconfig="tracking -C $CONFIGTABLEDEFAULT $arg_configure $arg_blindLayer -R $runNo -L $testlayer -n $nHitsGMax -x $t0shift0 -y $t0shift1 -l $tmin -u $tmax -g $geoSetup -s $sumCut -a $aaCut -B $iEntryStart -E $iEntryStop -w $workType -i $inputType -p $peakType $prerunname $temprunname > $logtemp 2> $errtemp"
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

    ana -C $CONFIGTABLEDEFAULT $arg_configure -R $runNo -x $XTTYPE -g $geoSetup -H $SAVEHISTS -i $inputType -c $maxchi2 -L $DEFAULTLAYER -n $NHITSMAX -o $tmaxSet $arg_wiremap $arg_xtfile $arg_draw $arg_xtopt -V Ana=info $prerunname $currunname $wires
    echo "ana -C $CONFIGTABLEDEFAULT $arg_configure -R $runNo -x $XTTYPE -g $geoSetup -H $SAVEHISTS -i $inputType -c $maxchi2 -L $DEFAULTLAYER -n $NHITSMAX -o $tmaxSet $arg_wiremap $arg_xtfile $arg_draw $arg_xtopt -V Ana=info $prerunname $currunname $wires"
done
