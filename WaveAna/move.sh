#!/bin/bash
mv *.pdf pdf
mv *.png png
cd pdf
./move.sh
cd ../png
./move.sh
cd ..
