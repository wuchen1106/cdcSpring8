#!/bin/bash
Name="job"
iStart=0
iStop=399

for (( i=iStart; i<=iStop; i++ ))
do
    jobname=$Name"_"$i
    touch $CDCS8WORKING_DIR/Conf/$jobname.conf
    thejob=$CDCS8WORKING_DIR/Jobs/$jobname.sh
    jobTemplate=$CDCS8WORKING_DIR/Jobs/template
    echo "#!/bin/bash" > $thejob
    echo "cd $CDCS8WORKING_DIR" >> $thejob
    cat $jobTemplate >> $thejob
    sed -i "s/_THEJOB_/$jobname/g" $thejob
    chmod +x $thejob
done
