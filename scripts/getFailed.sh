#!/bin/bash
justCheck=true
if [ $# -gt 0 ]
then
    justCheck=false
fi

for outfile in $CDCS8WORKING_DIR/Jobs/*.out
do
    if [ ! -e $outfile ]
    then
        continue;
    fi
    jobID=`gawk 'BEGIN{i=0;}{i++; if(i==2) print $0;}' $outfile | sed 's/.*job_\(\w*\)\.sh.*/\1/'`
    averageMem=`grep "Average Memory" $outfile| sed 's/.* \(.*\) MB.*/\1/'`
    log=`tail -n 1 $CDCS8WORKING_DIR/Conf/job_${jobID}.log | sed 's/.* > \(.*\.log\).*/\1/'`
    runID=`tail -n 1 $CDCS8WORKING_DIR/Conf/job_${jobID}.log | sed 's/.*tracking \(\w*\) .*/\1/'`
    runName=`tail -n 1 $CDCS8WORKING_DIR/Conf/job_${jobID}.log | sed 's/.*tracking \w* \w* [^ ]* \([^ ]*\) .*/\1/'`
    echo "$outfile"
    echo "  job_$jobID, run $runID, runName $runName"
    if [ -e $log ]
    then
        evtstart=`echo $log| sed 's/.*\.\(\w*\)-\(\w*\)\..*/\1/'`
        evtstop=`echo $log| sed 's/.*\.\(\w*\)-\(\w*\)\..*/\2/'`
        evtnow=`tail -n 1 $log`
        nProcessed=`echo "$evtnow-$evtstart" | bc -l`
        nTotal=`echo "$evtstop-$evtstart" | bc -l`
        checkID=`tail -n 1 $CDCS8WORKING_DIR/iter.${runID}.* | gawk '{print $1;}'`
        echo "  Still waiting @ $checkID! $evtstart-$evtstop, break at @ $evtnow, $nProcessed/$nTotal, average $averageMem MB"
        if ! $justCheck
        then
#            mv $outfile $CDCS8WORKING_DIR/Jobs/failed
            # mostly the job failed for memory leackage, which is difficult to fix, and probably will fail again with the same setup
            # Now we try to ignore this run when it's failed.
            echo "Good Events" >> $log
            rm $outfile
            rm $CDCS8WORKING_DIR/Conf/job_${jobID}.log
            echo "" > $CDCS8WORKING_DIR/Conf/job_${jobID}.conf
            cd $CDCS8WORKING_DIR/Jobs
            bsub -q h -o $PWD -e $PWD "./job_${jobID}.sh"
            cd ..
        fi
    else
        evtstart=0
        evtstop=0
        evtnow=0
        nProcessed=0
        nTotal=0
        echo "   Broken but skipped..."
        if ! $justCheck
        then
            mv $outfile $CDCS8WORKING_DIR/Jobs/failed
            cd $CDCS8WORKING_DIR/Conf
            rm job_${jobID}.log
            echo "" > job_${jobID}.conf
            cd $CDCS8WORKING_DIR/Jobs
            bsub -q h -o $PWD -e $PWD "./job_${jobID}.sh"
            cd $CDCS8WORKING_DIR
        fi
    fi
done
