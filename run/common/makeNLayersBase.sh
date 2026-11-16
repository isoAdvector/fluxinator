#!/bin/bash
# Copyright (C) 2026 Johan Roenby
# SPDX-License-Identifier: GPL-3.0-or-later
# Copies a base case and adds an nLayers entry, for the studies that vary it: parmScanner sets
# existing entries only. All other cases use the default, so that nobody copies the entry into a
# case on another mesh, where nLayers 1 is not enough (Part I).
# Usage: makeNLayersBase.sh <baseCase> <newCase>
type dictator > /dev/null 2>&1 || { echo "Source dictator first (github.com/isoAdvector/dictator)" >&2; exit 1; }
rm -rf "${2:?}"
cp -r "${1:?}" "$2"
dictator "$2/system/fvSchemes" interfaceForces.gravity.nLayers 1 > /dev/null || exit 1
