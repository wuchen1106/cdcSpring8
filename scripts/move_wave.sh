#!/bin/bash
prePWD=$PWD
pdf=$PWD/result/WaveAna/pdf
png=$PWD/result/WaveAna/png
mv *.pdf $pdf
mv *.png $png
cd $pdf
mv *.aa.* aa/
mv *.at.* at/
mv *.st.* st/
mv *.sum.* sum/
mv *.tdc.* tdc/
mv *.wf.* wf/
cd $png
mv *.aa.* aa/
mv *.at.* at/
mv *.st.* st/
mv *.sum.* sum/
mv *.tdc.* tdc/
mv *.wf.* wf/
cd $prePWD
