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
    echo "$outfile"
    echo "  job_$jobID"
    echo "  average $averageMem MB"
    if ! $justCheck
    then
#       mv $outfile $CDCS8WORKING_DIR/Jobs/failed
        # mostly the job failed for memory leackage, which is difficult to fix, and probably will fail again with the same setup
        # Now we try to ignore this run when it's failed.
        rm $outfile
        rm $CDCS8WORKING_DIR/Conf/job_${jobID}.log
        echo "" > $CDCS8WORKING_DIR/Conf/job_${jobID}.conf
        cd $CDCS8WORKING_DIR/Jobs
        bsub -q h -o $PWD -e $PWD "./job_${jobID}.sh"
        cd ..
    fi
done
