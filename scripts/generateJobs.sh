#!/bin/bash
runNo=100007
N=90457
layers="1 2 3 4 5 6 7 8"

for i in $layers;
do
    testlayer=$i
    for (( j=0; j<N; j+=5000 ))
    do
        iEntryStart=$j
        iEntryStop=$((j+4999))
        if (( iEntryStop>=N ))
        then
            iEntryStop=$((N-1))
        fi
        jobname=$runNo.layer$i.$iEntryStart-$iEntryStop
        touch $CDCS8WORKING_DIR/Conf/$jobname.conf
        thejob=$CDCS8WORKING_DIR/Jobs/$jobname.sh
        jobTemplate=$CDCS8WORKING_DIR/Jobs/template
        echo "#!/bin/bash" > $thejob
        echo "cd $CDCS8WORKING_DIR" >> $thejob
        cat $jobTemplate >> $thejob
        sed -i "s/_THEJOB_/$jobname/g" $thejob
        chmod +x $thejob
    done
done
