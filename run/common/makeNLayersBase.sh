#!/bin/bash
# Copies a base case and adds an nLayers entry, for the studies that vary it: parmScanner sets
# existing entries only. All other cases use the default, so that nobody copies the entry into a
# case on another mesh, where nLayers 1 is not enough (Part I).
# Usage: makeNLayersBase.sh <baseCase> <newCase>
rm -rf "${2:?}"
cp -r "${1:?}" "$2"
sed -i '/^gravityModel/a nLayers                 1;' "$2/constant/gravityProperties"
