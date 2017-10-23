#!/bin/bash
runNo=118
inName="mergepeaks"
outName="newmerge"
echo "i/I c/I f/I s/I n/I cx cz" >  ${runNo}.${outName}.chi2
./finding4 $runNo 4 $inName 0 0 1000 15 -50 1000 0 >> ${runNo}.${outName}.chi2
#vi ${runNo}.${outName}.chi2
txt2root  ${runNo}.${outName}.chi2  ${runNo}.${outName}.chi2.root
./ana_finding_chi2  ${runNo}.${outName}.chi2.root >  ${runNo}.${outName}.minchi2
txt2root ${runNo}.${outName}.minchi2 ${runNo}.${outName}.minchi2.root
