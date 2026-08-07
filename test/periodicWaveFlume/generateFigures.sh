#!/bin/sh

scan=gravityModelScan

mkdir $scan/figs
rm $scan/figs/*.png

python3 compareResults.py --scanDir $scan --caseNames 0000 0001 0002 0036 --figName $scan/figs/00
python3 compareResults.py --scanDir $scan --caseNames 0003 0004 0005 0036 --figName $scan/figs/01
python3 compareResults.py --scanDir $scan --caseNames 0006 0007 0008 0036 --figName $scan/figs/02

python3 compareResults.py --scanDir $scan --caseNames 0009 0010 0011 0036 --figName $scan/figs/03
python3 compareResults.py --scanDir $scan --caseNames 0012 0013 0014 0036 --figName $scan/figs/04
python3 compareResults.py --scanDir $scan --caseNames 0015 0016 0017 0036 --figName $scan/figs/05

python3 compareResults.py --scanDir $scan --caseNames 0018 0019 0020 0036 --figName $scan/figs/06
python3 compareResults.py --scanDir $scan --caseNames 0021 0022 0023 0036 --figName $scan/figs/07
python3 compareResults.py --scanDir $scan --caseNames 0024 0025 0026 0036 --figName $scan/figs/08

python3 compareResults.py --scanDir $scan --caseNames 0027 0028 0029 0036 --figName $scan/figs/09
python3 compareResults.py --scanDir $scan --caseNames 0030 0031 0032 0036 --figName $scan/figs/10
python3 compareResults.py --scanDir $scan --caseNames 0033 0034 0035 0036 --figName $scan/figs/11

img2pdf $scan/figs/*.png -o ${scan}.pdf

#Attempt to see show each solver separately. Too many curves on top of each other
#rm $scan/figs/*.png
#python3 compareResults.py --scanDir $scan --caseNames 0000 0003 0006 0009 0012 0015 0018 0021 0024 --figName $scan/figs/00
#python3 compareResults.py --scanDir $scan --caseNames 0001 0004 0007 0010 0013 0016 0019 0022 0025 --figName $scan/figs/01
#python3 compareResults.py --scanDir $scan --caseNames 0002 0005 0008 0011 0014 0017 0020 0023 0026 --figName $scan/figs/02
#convert $scan/figs/*.png ${scan}_solvers.pdf
