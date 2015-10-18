#!/bin/bash
runNo=$1
runName=$2
#mkdir ${runNo}.${runName}
#mkdir ${runNo}.${runName}/dt
#mkdir ${runNo}.${runName}/xt
#./get_xt $runNo $runName
#mv *.png ${runNo}.${runName}/dt
./check_xt $runNo $runName
mv *.png ${runNo}.${runName}/xt
#./get_wire $runNo $runName
echo -e "\a"
