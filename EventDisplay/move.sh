#!/bin/bash
if [ $# -lt 1 ]; then exit; fi
runNo=$1
mkdir l9/$runNo
mkdir s9/$runNo
mv *.l9.pdf l9/$runNo
mv *.s9.pdf s9/$runNo
