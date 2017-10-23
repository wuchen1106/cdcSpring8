#!/bin/bash

for i in 2006 2013 2020 2047 2046 2035 2057
do
	if [ $i -lt 100 ]
	then
		s=`du -s ../root/run_0000"$i"_built.root|sed 's/\(\w*\) *.*/\1/'`
	elif [ $i -lt 1000 ]
	then
		s=`du -s ../root/run_000"$i"_built.root|sed 's/\(\w*\) *.*/\1/'`
	else
		s=`du -s ../root/run_00"$i"_built.root|sed 's/\(\w*\) *.*/\1/'`
	fi
	echo "$i: $s KB"
	if [ $s -gt 8 ]
	then
		echo "	larger than 8 KB"
		for (( j=1; j>0; j++ ))
		do
			sleep 2
			n=`ps -ef | grep getH | wc -l`
			echo "		$j: $n"
			if [ $n -lt 3 ]
			then
				echo "submit $i"
				nohup ./getH $i &
				break
			fi
		done
	fi
done

#nohup ./getD 234 &
