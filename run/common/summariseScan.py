#!/usr/bin/env python3
# Copyright (C) 2026 Johan Roenby
# SPDX-License-Identifier: GPL-3.0-or-later
#
# Prints a summary per case of a parmScanner scan of the hydrostatic box
# tests: the scanned parameters, the initial residual of the first p_rgh solve,
# and max|U| after the first and the last time step and its peak over the run,
# with 8 significant digits.
#
# Usage (from the test directory holding the scan):
#   python3 ../summariseScan.py <scanDir>

import os
import re
import sys


def first_prgh_residual(logFile):
    if not os.path.isfile(logFile):
        return None
    with open(logFile) as f:
        for line in f:
            if 'Solving for p_rgh' in line:
                m = re.search(r'Initial residual = ([^,]+),', line)
                return m.group(1) if m else None
    return None


def max_u(caseDir):
    f = os.path.join(caseDir, 'postProcessing', 'fieldMinMax1', '0', 'fieldMinMax.dat')
    if not os.path.isfile(f):
        return None
    rows = []
    with open(f) as fh:
        for line in fh:
            if line.startswith('#') or not line.strip():
                continue
            cols = line.split()
            rows.append((cols[0], float(cols[2])))
    return rows


def main():
    if len(sys.argv) != 2:
        sys.exit('Usage: summariseScan.py <scanDir>')
    scanDir = sys.argv[1]
    with open(os.path.join(scanDir, 'parmScannerParms')) as f:
        lines = [line.rstrip() for line in f if line.strip()]

    print(f'Scan {scanDir}')
    for line in lines:
        caseDir, parms = line.split(':', 1)
        parms = ' '.join(parms.split())
        m = re.search(r'application:\s*(\S+)', parms)
        application = m.group(1) if m else 'interIsoGoam'
        residual = first_prgh_residual(os.path.join(caseDir, 'log.' + application))
        rows = max_u(caseDir)
        if rows:
            peak = max(r[1] for r in rows)
            u = (f'max|U| t={rows[0][0]}: {rows[0][1]:.7e}  t={rows[-1][0]}: {rows[-1][1]:.7e}'
                 f'  peak: {peak:.7e}')
        else:
            u = 'no fieldMinMax output'
        print(f'{os.path.basename(caseDir)}  {parms}')
        print(f'      first p_rgh residual {residual}  {u}')


if __name__ == '__main__':
    main()
