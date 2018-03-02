#!/bin/bash

for (( i=0; i<3600; i++ ))
do
    newFailed=false
    for outfile in $CDCS8WORKING_DIR/Jobs/*.out
    do
        if [ ! $outfile == "$CDCS8WORKING_DIR/Jobs/\*.out" ]
        then
            newFailed=true
        fi
    done
    if $newFailed
    then
        getFailed.sh yes
    fi
    sleep 10
done
