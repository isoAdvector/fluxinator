#!/bin/bash

scanNames="gravityModelScan_1deltaT gravityModelScan_100deltaT gravityModelScan_1000deltaT"

source ../cfdTools/cfdTools

for scanName in $scanNames; do
    echo $scanName
    rm -rf $scanName
    parmScanner -file ${scanName}.csv -exec Allrun -parallel &
done
