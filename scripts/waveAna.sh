#!/bin/bash
startName="Garfield"
bin=$1

which $bin
if [ ! $? -eq 0 ]
then
    bin="NONE"
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

cd $CDCS8WORKING_DIR/info
for i in $list
do
    if [ $bin == "getP" ]
    then
        file=`printf $CDCS8WORKING_DIR/root/run_%06d_built.root $i`
    else
        file=$CDCS8WORKING_DIR/root/p_$i.root
    fi
	if [ -e $file ]
	then
	    if [ ! $bin == "NONE" ]
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
		fi
        if [ -e $CDCS8WORKING_DIR/info/wire-position.${i}.${startName}.root ]
        then
            rm $CDCS8WORKING_DIR/info/wire-position.${i}.${startName}.root
        fi
        if [ -e $CDCS8WORKING_DIR/info/xt.${i}.${startName}.root ]
        then
            rm $CDCS8WORKING_DIR/info/xt.${i}.${startName}.root
        fi
        for runlist in list.*.*
        do
            if grep -q "\<$i\>" $runlist
            then
                gas=`echo $runlist | sed 's/list.\(.*\)/\1/'`
                break
            fi
        done
        ln -s wire-position.root wire-position.${i}.${startName}.root
        ln -s xt.${gas}.root xt.${i}.${startName}.root
    else
        echo "Cannot find file \"$file\""
	fi
done
cd -
