#!/bin/bash

collectedFigName=tiltedBox
#scanNames="originalBehaviour reconstructGhf_isoAlpha reconstructGhf_plicRDF hardcodedHeight"
scanNames="gravityModelScan_1deltaT gravityModelScan_100deltaT gravityModelScan_1000deltaT"

mkdir -p figures

for scanName in $scanNames; do
    echo $scanName
    nFigs=$(ls ${scanName}/*/*.png | wc -l)
    inFiles=$(eval ls ${scanName}/*/*.png)
    montage ${inFiles} -tile ${nFigs}x1 -geometry +0+0 figures/${scanName}.png
done

outFile=${collectedFigName}.png
rm $outFile
nFigs=$(ls figures/*.png | wc -l)
inFiles=$(eval ls -tr figures/*.png)
montage $inFiles -tile 1x${nFigs} -geometry +0+0 $outFile
xdg-open $outFile &
