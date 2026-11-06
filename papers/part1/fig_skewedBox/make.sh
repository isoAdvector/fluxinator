#!/bin/sh
# Copyright (C) 2026 Johan Roenby
# SPDX-License-Identifier: GPL-3.0-or-later
# Velocity fields in the skewed box after 10 time steps, corrected and uncorrected schemes.
# Inputs: run/skewedBox/schemes_10steps/000[0-5]/000[0-5].png (written in each case by
# pvbatch generateFigure.py, see run/skewedBox/Allrun). Output: skewedBox.pdf
cd "${0%/*}" || exit 1
pdflatex -interaction=nonstopmode skewedBox.tex > log.make 2>&1 || { tail log.make; exit 1; }
rm -f skewedBox.aux skewedBox.log
