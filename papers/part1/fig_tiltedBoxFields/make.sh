#!/bin/sh
# Copyright (C) 2026 Johan Roenby
# SPDX-License-Identifier: GPL-3.0-or-later
# Velocity fields in the tilted box after 1, 100 and 1000 time steps.
# Inputs: run/tiltedBox/gravityModels_{1step,100steps,1000steps}/000[0-3]/000[0-3].png
# (written by each case, run/tiltedBox/baseCase/generateFigure.py).
# Output: tiltedBoxFields.pdf
cd "${0%/*}" || exit 1
pdflatex -interaction=nonstopmode tiltedBoxFields.tex > log.make 2>&1 || { tail log.make; exit 1; }
rm -f tiltedBoxFields.aux tiltedBoxFields.log
