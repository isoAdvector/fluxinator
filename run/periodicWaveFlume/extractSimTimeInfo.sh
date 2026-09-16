#!/bin/sh

for n in $(ls -d ${1}/0*); 
do
    echo $n
    cd $n
    grep "ClockTime = "  log.*oam | cut -d' ' -f8 | tail -1 > ClockTime
    grep "^Time = "  log.*oam | wc -l > nTimeSteps
    cd -
done
