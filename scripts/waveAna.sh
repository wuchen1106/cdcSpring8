#!/bin/bash
startName="Garfield"
bin=$1

which $bin
if [ ! $? -eq 0 ]
then
    echo "Cannot find \"$bin\"!"
    exit 0
fi

if [ $# -lt 2 ]
then
    echo "Please use:"
    echo "   $0 BIN list.GAS.HV"
    echo " or"
    echo "   $0 BIN runNo1 runNo2 ..."
    exit 0
fi

if [ -e $2 ]
then
    list=`cat $2`
else
    for (( i=2; i<=$#; i++ ))
    do
        list="$list ${!i}"
    done
fi

for i in $list
do
	file=`printf $CDCS8WORKING_DIR/root/run_%06d_built.root $i`
	if [ -e $file ]
	then
		for (( j=1; j>0; j++ ))
		do
			sleep 2
			n=`ps -ef | grep $bin | wc -l`
			echo "		$j: $n"
			if [ $n -lt 4 ]
			then
				echo "submit $i"
				nohup $bin $i &
				break
			fi
		done
        if [ ! -e $CDCS8WORKING_DIR/info/wire-position.${i}.${startName}.root ]
        then
            cd $CDCS8WORKING_DIR/info
            ln -s wire-position.root wire-position.${i}.${startName}.root
            cd -
        fi
        if [ ! -e $CDCS8WORKING_DIR/info/xt.${i}.${startName}.root ]
        then
            cd $CDCS8WORKING_DIR/info
            ln -s xt.C4H10.1800.root xt.${i}.${startName}.root
            cd -
        fi
    else
        echo "Cannot find file \"$file\""
	fi
done
