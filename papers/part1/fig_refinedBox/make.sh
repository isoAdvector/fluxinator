#!/bin/sh
# Velocity magnitude in the box with a refined region, z_w = 0.49, t = 0.1 s.
# Inputs: run/refinedBox/waterLevel_10steps/{iIF,cH,g1,g2}_hyd_z0.4900/*.png (written in
# each case by pvbatch generateFigure.py, see run/refinedBox/Allrun). Output: refinedBox.pdf
cd "${0%/*}" || exit 1
pdflatex -interaction=nonstopmode refinedBox.tex > log.make 2>&1 || { tail log.make; exit 1; }
rm -f refinedBox.aux refinedBox.log
