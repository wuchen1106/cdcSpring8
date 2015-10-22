#!/bin/bash
##for i in 1 3 5 7;
#for i in 1 2 3 4 5 6 7;
#do
##	nohup ./finding 117 $i layer$i.garfold.0 &
##	nohup ./myminuit 117 $i layer$i.garfold.0 &
##	nohup ./finding 117 $i layer$i.garf.0 &
##	nohup ./myminuit 117 $i layer$i.garf.0.difchi2 &
##	nohup ./finding 117 $i layer$i.garf.4.k3 -4 &
##	nohup ./myminuit 117 $i layer$i.garf.4.k3 &
##	nohup ./finding 117 $i layer$i.garf.8 -8 &
##	nohup ./myminuit 117 $i layer$i.garf.8 &
##	nohup ./finding 117 $i layer$i.garf.8.j4 -8 &
##	nohup ./myminuit 117 $i layer$i.garf.8.j4 &
##	nohup ./myminuit 117 $i layer$i.garf.8.j3.5hits &
#done

runName="t0fit.1xt.8"
for i in 1 2 3 4 5 6 7 8;
do
#	nohup ./finding2 117 $runName $i &
	nohup ./myminuit2 117 $i $runName &
done
tail -f nohup.out
echo -e "\a"
