#!/bin/bash
if [ $# -lt 2 ]
then
    exit 0
fi
logfile=$1
output=$2
temp=/tmp/trans.analog

title=`grep "^=>" $logfile | tail -n 1 | sed 's/=>//'`
echo "$title avEtrackx0 avEtrackd0 avresInix0 avresInid0 avEtrackx avEtrackd avresInix avresInid" > $output

touch $temp
grep "^==>" $logfile |sed "s/==>//g"| gawk '{if ($5==0) print $0;}' | sort -k 4 -n | sort -k 3 >> $temp
grep "^==>" $logfile |sed "s/==>//g"| gawk '{if ($5==1) print $0;}' | sort -k 4 -n | sort -k 3 >> $temp
grep "^==>" $logfile |sed "s/==>//g"| gawk '{if ($5==2) print $0;}' | sort -k 4 -n | sort -k 3 >> $temp
NRUNS=`cat $temp| wc -l`
for (( i=1; i<=$NRUNS; i++ ))
do
    runNo=`gawk -v gawk_keyword=$i 'BEGIN{i=1;}{if (i==gawk_keyword) print $1; i++;}' $temp`
    line=`gawk -v gawk_keyword=$i 'BEGIN{i=1;}{if (i==gawk_keyword) print $0; i++;}' $temp`
    resfile=$CDCS8WORKING_DIR/result/MCreso/res.$runNo.*
    if [ ! -e $resfile ]
    then
        echo "Cannot find $resfile!"
        continue
    fi
    avEtrackx0=`cat $resfile | grep "avEtrackx" | tail -n 2 | sed "s/,//" | gawk 'BEGIN{i=1;}{if (i==1) {print $2} i++;}'`
    avEtrackd0=`cat $resfile | grep "avEtrackx" | tail -n 2 | sed "s/,//" | gawk 'BEGIN{i=1;}{if (i==1) {print $4} i++;}'`
    avEtrackx=`cat $resfile | grep "avEtrackx" | tail -n 2 | sed "s/,//" | gawk 'BEGIN{i=1;}{if (i==2) {print $2} i++;}'`
    avEtrackd=`cat $resfile | grep "avEtrackx" | tail -n 2 | sed "s/,//" | gawk 'BEGIN{i=1;}{if (i==2) {print $4} i++;}'`
    avresInix0=`cat $resfile | grep "avresInix" | tail -n 2 | sed "s/,//" | gawk 'BEGIN{i=1;}{if (i==1) {print $2} i++;}'`
    avresInid0=`cat $resfile | grep "avresInix" | tail -n 2 | sed "s/,//" | gawk 'BEGIN{i=1;}{if (i==1) {print $4} i++;}'`
    avresInix=`cat $resfile | grep "avresInix" | tail -n 2 | sed "s/,//" | gawk 'BEGIN{i=1;}{if (i==2) {print $2} i++;}'`
    avresInid=`cat $resfile | grep "avresInix" | tail -n 2 | sed "s/,//" | gawk 'BEGIN{i=1;}{if (i==2) {print $4} i++;}'`
    echo "$line $avEtrackx0 $avEtrackd0 $avresInix0 $avresInid0 $avEtrackx $avEtrackd $avresInix $avresInid" >> $output
done

rm $temp
