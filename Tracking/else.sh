#!/bin/bash

runNo=117
runName="t0fit.1xt.efuc.1"

cd ../Ana
./get_xt ${runNo} ${runName}
mkdir ${runNo}.${runName}
mkdir ${runNo}.${runName}/xt
mv xt*.png ${runNo}.${runName}/xt
cd ../info
rm xt.${runNo}.root
ln -s xt.${runNo}.${runName}.root xt.${runNo}.root
cd ../Tracking
