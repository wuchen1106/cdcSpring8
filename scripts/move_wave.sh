#!/bin/bash
prePWD=$PWD
pdf=$PWD/results/WaveAna/pdf
png=$PWD/results/WaveAna/png
mv *.pdf $pdf
mv *.png $png
cd $pdf
./move.sh
cd $png
./move.sh
cd $prePWD
