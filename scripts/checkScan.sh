#!/bin/bash

for list in Input/list.C*.*
do
    setup=`echo $list|sed 's/info\/list\.\(C\w*\)\.\(\w*\)/\1\.\2/'`
    gas=`echo $list|sed 's/info\/list\.\(C\w*\)\.\(\w*\)/\1/'`
    HV=`echo $list|sed 's/info\/list\.\(C\w*\)\.\(\w*\)/\2/'`
    for run in `cat $list`
    do
        for file in root/ana_$run.*.i1.layer4.root
        do
            if [ ! -e $file ]
            then
                continue
            fi
            runname=`echo $file | sed 's/root\/ana_\w*\.\(.*\).i1.layer4.root/\1/'`
            for (( iter=20; iter>=0; iter-- ))
            do
                thefile=root/ana_$run.$runname.i$iter.layer4.root;
                if [ -e $thefile ]
                then
                    break;
                fi
            done
            if [ ! -e $thefile ]
            then
                continue
            fi
# if you want to remove intermidiate steps
#            if [ $iter -gt 2 ]
#            then
#                for (( i=2; i<$iter; i++ ))
#                do
#                    rm root/*_$run.$runname.i$i.layer4.root
#                done
#            fi
            echo "$gas $HV $run $runname $thefile"
        done
    done
done
