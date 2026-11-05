#!/bin/bash
# summarise.sh <case>: max over time of the per-step errors per variant
for f in $1/postProcessing/convectionTest/*.dat; do printf "%-18s %s\n" "$(basename $f .dat)" "$(awk 'NR>2{for(i=2;i<=7;i++) if($i>m[i]) m[i]=$i; s+=$5; n++} END{printf "maxAll %.2e  water %.2e  interface %.2e (rms mean %.2e)  air %.2e  momentum %.2e", m[2],m[3],m[4],s/n,m[6],m[7]}' $f)"; done
