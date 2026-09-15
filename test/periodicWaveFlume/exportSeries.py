#!/usr/bin/env python
# Exports time series excerpts for the report pages as a JS data file.
# Usage: exportSeries.py <scanDir> <varName> <out.js>
import json
import os
import re
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from plotCrestHeight import find_extrema, read_parms, wave_heights  # noqa: E402

U_CREST = 0.44453480326804


def r2(a, d=2):
    return [round(float(x), d) for x in a]


def export_case(caseDir):
    t, h = np.loadtxt(os.path.join(caseDir, 'postProcessing/interfaceHeight1/0/height.dat'),
                      comments='#', unpack=True, usecols=(0, 2))
    tu, u = np.loadtxt(os.path.join(caseDir, 'postProcessing/fieldMinMax1/0/fieldMinMax.dat'),
                       comments='#', unpack=True, usecols=(0, 2))
    out = {'end': round(float(t[-1]), 3)}
    thr = 0.1*(h.max() - h.min())
    crests = find_extrema(t, h, +1, thr)
    troughs = find_extrema(t, h, -1, thr)
    H = wave_heights(crests, troughs)
    out['crests'] = [[round(float(a), 3), round(float(b)*1e3, 2)] for a, b in crests]
    out['troughs'] = [[round(float(a), 3), round(float(b)*1e3, 2)] for a, b in troughs]
    out['H'] = [[round(float(a), 3), round(float(b)*1e3, 2)] for a, b in H]

    # Surface elevation over the first and last 4 s, 5 ms resolution [mm]
    for key, (a, b) in (('first', (0.0, 4.0)), ('last', (26.0, 30.0))):
        tt = np.arange(a, b + 1e-9, 5e-3)
        if t[-1] < b:
            out[key] = None
            continue
        out[key] = r2(np.interp(tt, t, h)*1e3)

    # Crest profiles: +-0.3 s around each crest, 5 ms resolution [mm]
    tw = np.arange(-0.3, 0.3 + 1e-9, 5e-3)
    out['crestProfiles'] = [r2(np.interp(tc + tw, t, h)*1e3) for tc, _ in crests
                            if tc - 0.3 > t[0] and tc + 0.3 < t[-1]]

    # max|U|/u_crest in 0.1 s bins: bin max and bin median
    edges = np.arange(0, 30 + 1e-9, 0.1)
    idx = np.digitize(tu, edges) - 1
    umaxBin, umedBin = [], []
    for i in range(len(edges) - 1):
        v = u[idx == i]
        if len(v) == 0:
            umaxBin.append(None)
            umedBin.append(None)
        else:
            umaxBin.append(round(float(v.max()/U_CREST), 3))
            umedBin.append(round(float(np.median(v)/U_CREST), 3))
    out['umaxBin'] = umaxBin
    out['umedBin'] = umedBin

    # Mean time step in 0.1 s bins [ms], from the number of steps written per bin
    # (the time column is printed with too few digits to difference it directly)
    counts = np.bincount(np.digitize(t, edges) - 1, minlength=len(edges) - 1)[:len(edges) - 1]
    out['dtBin'] = [round(100.0/int(n), 4) if n > 0 else None for n in counts]
    return out


def main():
    scanDir, varName, outFile = sys.argv[1:4]
    series = {}
    for c in sorted(os.listdir(scanDir)):
        if not re.fullmatch(r'\d{4}', c):
            continue
        caseDir = os.path.join(scanDir, c)
        if not os.path.isfile(os.path.join(caseDir, 'parmScannerParms')):
            continue
        series[c] = export_case(caseDir)
        print(c, read_parms(caseDir).get('application'))
    with open(outFile, 'w') as f:
        f.write(f'window.{varName} = ' + json.dumps(series, separators=(',', ':')) + ';\n')
    print('wrote', outFile, os.path.getsize(outFile)//1024, 'kB')


if __name__ == '__main__':
    main()
