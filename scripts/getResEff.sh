#!/bin/bash
startName="Garfield"
options="6 0 1 35"

if [ $# -lt 3 ]
then
    echo "Please use:"
    echo "   $0 runName lid list.GAS.HV"
    echo " or"
    echo "   $0 runName lid runNo1 runNo2 ..."
    exit 0
fi

name=$1
lid=$2

if [ -e $3 ]
then
    list=`cat $3`
else
    for (( i=3; i<=$#; i++ ))
    do
        list="$list ${!i}"
    done
fi

cd $CDCS8WORKING_DIR/info
for i in $list
do
    file=$CDCS8WORKING_DIR/root/ana_${i}.${name}.layer${lid}.root
	if [ -e $file ]
	then
        for (( j=1; j>0; j++ ))
        do
            sleep 2
            n=`ps -ef | grep "\<ana\>" | wc -l`
            echo "		$j: $n"
            if [ $n -lt 4 ]
            then
                echo "submit $i"
                nohup ana $i $name $lid $options &
                break
            fi
        done
    else
        echo "Cannot find file \"$file\""
	fi
done
cd -
