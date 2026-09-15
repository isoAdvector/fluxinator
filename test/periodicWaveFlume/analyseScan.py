#!/usr/bin/env python
#
# Computes accuracy and efficiency metrics for all 4-digit cases in a scan
# directory, compared against a reference case (by default the analytic
# streamFunFoam run). Prints a table and optionally writes it as CSV.
#
# Metrics (H = wave height = crest - following trough, see plotCrestHeight.py):
#   avgMaxCo    time average of the domain max Courant number (from the log)
#   end         last time in height.dat (< endTime means the run stopped early)
#   clockTime   last ClockTime in the log [s], steps = number of time steps
#   msPerStep   clockTime/steps [ms]
#   H1err       error of the first wave height relative to the exact H [%]
#   dH10, dH20  wave height error at t ~ 10 s and t ~ 20 s [%]
#   dHend       error of the last wave height [%]
#   dHfit       wave height error at the last wave from a quadratic fit through
#               all per-wave errors [%]. Less sensitive to wave-to-wave scatter
#               than dHend.
#   Hscatter    rms deviation of the per-wave errors from that fit [%-points]
#   crestEnd    error of the last crest height [%]
#   troughEnd   error of the last trough height [%] (+ means shallower trough)
#   crestPeaks  mean number of peaks (prominence > 0.5 mm) within +-0.3 s of
#               each crest after the first; a smooth crest gives 1
#   crestScatter rms deviation of the per-crest height errors from a quadratic
#               fit through all crests [%-points]; wave-to-wave irregularity
#   shapeErr    rms surface elevation error over the last 4 s after removing the
#               best time shift, divided by H [%] (amplitude and shape, no phase)
#   phase       time of last crest minus the reference crest time [s]
#               (- means crests arrive early)
#   rmsEta      rms of surface elevation error over the record, divided by H
#   Umean, Uend time mean and final max|U| divided by the exact crest velocity
#   Umed, U95   median and 95th percentile of max|U|/u_crest for t > 2 s
#
# Usage example:
#   ./analyseScan.py --scanDir gravityModelScan_fixedDt --csv fixedDtMetrics.csv

import argparse
import csv
import os
import re

import numpy as np
from scipy.signal import find_peaks

from plotCrestHeight import find_extrema, read_log, read_parms, wave_heights

U_CREST = 0.44453480326804  # Exact crest velocity of the stream function wave
PERIOD = 2.0                # Wave period [s]


def load(scanDir, case):
    caseDir = os.path.join(scanDir, case)
    t, h = np.loadtxt(os.path.join(caseDir, 'postProcessing/interfaceHeight1/0/height.dat'),
                      comments='#', unpack=True, usecols=(0, 2))
    threshold = 0.1*(h.max() - h.min())
    crests = find_extrema(t, h, +1, threshold)
    troughs = find_extrema(t, h, -1, threshold)
    tu, maxU = np.loadtxt(os.path.join(caseDir, 'postProcessing/fieldMinMax1/0/fieldMinMax.dat'),
                          comments='#', unpack=True, usecols=(0, 2))
    parms = read_parms(caseDir)
    application = parms.get('application')
    clockTime, meanMaxCo = read_log(caseDir, application)
    nSteps = None
    logFile = os.path.join(caseDir, 'log.' + str(application))
    if os.path.isfile(logFile):
        with open(logFile) as file:
            nSteps = sum(1 for line in file if line.startswith('Time = '))
    return dict(case=case, parms=parms, t=t, h=h, crests=crests, troughs=troughs,
                H=wave_heights(crests, troughs), tu=tu, maxU=maxU, clockTime=clockTime,
                meanMaxCo=meanMaxCo, nSteps=nSteps)


def crest_peaks(t, h, crests, halfWidth=0.3, prominence=5e-4):
    counts = []
    for tc, _ in crests[1:]:
        m = abs(t - tc) < halfWidth
        if m.sum() < 5:
            continue
        # Pad with the window minimum so peaks at the window edge are not counted
        hw = np.concatenate(([h[m].min()], h[m], [h[m].min()]))
        peaks, _ = find_peaks(hw, prominence=prominence)
        counts.append(len(peaks))
    return float(np.mean(counts)) if counts else None


def fit_scatter(x, y):
    # rms deviation of y from a quadratic fit in x
    return float(np.sqrt(np.mean((y - np.polyval(np.polyfit(x, y, 2), x))**2)))


def shape_error(t, h, ref, tEnd, window=2*PERIOD, dt=4e-3, maxShift=0.5):
    tt = np.arange(tEnd - window, tEnd, dt)
    hh = np.interp(tt, t, h)
    best = np.inf
    for s in np.arange(-maxShift, maxShift + 1e-9, 1e-3):
        err = hh - np.interp(tt - s, ref['t'], ref['h'])
        best = min(best, np.sqrt(np.mean(err**2)))
    return float(best)


def metrics(r, ref, endTime):
    Hexact = ref['H'][:, 1].mean()
    crestExact = ref['crests'][:, 1].mean()
    troughExact = ref['troughs'][:, 1].mean()

    co = r['meanMaxCo']
    diverged = co is not None and (not np.isfinite(co) or co > 1e3)
    row = dict(case=r['case'])
    row.update(r['parms'])
    row.update(avgMaxCo=None if (co is None or diverged) else co,
               diverged=diverged, end=r['t'][-1], clockTime=r['clockTime'],
               steps=r['nSteps'])
    if r['clockTime'] is not None and r['nSteps']:
        row['msPerStep'] = r['clockTime']/r['nSteps']*1e3

    H = r['H']
    if len(H) < 3 or r['t'][-1] < endTime - 1:
        return row  # Incomplete run: accuracy metrics not meaningful

    Hrel = (H[:, 1] - Hexact)/Hexact*100

    def at(time):
        return Hrel[np.argmin(abs(H[:, 0] - time))]

    fit = np.polyfit(H[:, 0], Hrel, 2)
    n = min(len(r['crests']), len(ref['crests']))
    etaRef = np.interp(r['t'], ref['t'], ref['h'])
    late = r['tu'] > 2
    row.update(H1err=Hrel[0], dH10=at(10), dH20=at(20), dHend=Hrel[-1],
               dHfit=np.polyval(fit, H[-1, 0]),
               Hscatter=np.sqrt(np.mean((Hrel - np.polyval(fit, H[:, 0]))**2)),
               crestEnd=(r['crests'][-1, 1] - crestExact)/crestExact*100,
               troughEnd=(r['troughs'][-1, 1] - troughExact)/abs(troughExact)*100,
               crestPeaks=crest_peaks(r['t'], r['h'], r['crests']),
               crestScatter=fit_scatter(r['crests'][:, 0],
                                        (r['crests'][:, 1] - crestExact)/crestExact*100),
               shapeErr=100*shape_error(r['t'], r['h'], ref, endTime)/Hexact,
               phase=r['crests'][n-1, 0] - ref['crests'][n-1, 0],
               rmsEta=np.sqrt(np.mean((r['h'] - etaRef)**2))/Hexact,
               Umean=r['maxU'].mean()/U_CREST, Uend=r['maxU'][-1]/U_CREST,
               Umed=np.median(r['maxU'][late])/U_CREST,
               U95=np.percentile(r['maxU'][late], 95)/U_CREST)
    return row


def fmt(value, spec):
    if value is None:
        return '-'
    return format(value, spec)


def main():
    parser = argparse.ArgumentParser(description='Accuracy and efficiency metrics for a scan.')
    parser.add_argument('--scanDir', required=True, help='Directory where cases are found')
    parser.add_argument('--refCase', default='0036',
                        help='Reference (analytic) case. Default: 0036')
    parser.add_argument('--csv', default=None, help='Also write the table to this CSV file')
    args = parser.parse_args()

    cases = sorted(c for c in os.listdir(args.scanDir)
                   if re.fullmatch(r'\d{4}', c)
                   and os.path.isfile(os.path.join(args.scanDir, c, 'parmScannerParms')))
    data = {c: load(args.scanDir, c) for c in cases}
    ref = data[args.refCase]
    endTime = ref['t'][-1]

    refRow = metrics(ref, ref, endTime)
    print(f"Reference {args.refCase}: H = {ref['H'][:, 1].mean():.6f}, "
          f"crest = {ref['crests'][:, 1].mean():.6f}, trough = {ref['troughs'][:, 1].mean():.6f}, "
          f"period = {np.diff(ref['crests'][:, 0]).mean():.5f} s, "
          f"max|U|/u_crest = {ref['maxU'].mean()/U_CREST:.3f}, crestPeaks = {refRow['crestPeaks']:.2f}, "
          f"crestScatter = {refRow['crestScatter']:.4f}")

    rows = [metrics(data[c], ref, endTime) for c in cases if c != args.refCase]

    # Scan parameters other than the solver, in order of appearance
    parmNames = []
    for r in rows:
        parmNames += [p for p in data[r['case']]['parms'] if p not in parmNames and p != 'application']

    columns = [('case', 's'), ('application', '12s')] + [(p, '>9s') for p in parmNames] + [
        ('avgMaxCo', '8.3f'), ('end', '6.2f'), ('clockTime', '9.0f'), ('steps', '7d'),
        ('msPerStep', '9.2f'), ('H1err', '6.2f'), ('dH10', '7.2f'), ('dH20', '7.2f'),
        ('dHend', '7.2f'), ('dHfit', '7.2f'), ('Hscatter', '8.2f'), ('crestEnd', '8.2f'),
        ('troughEnd', '9.2f'), ('crestPeaks', '10.2f'), ('crestScatter', '12.2f'),
        ('shapeErr', '8.2f'), ('phase', '7.3f'), ('rmsEta', '7.3f'), ('Umean', '6.2f'),
        ('Uend', '5.2f'), ('Umed', '5.2f'), ('U95', '5.2f')]
    widths = {name: max(len(name), len(format(0 if spec[-1] in 'fd' else '', spec)))
              for name, spec in columns}
    print("  ".join(name.rjust(widths[name]) for name, _ in columns))
    for r in rows:
        cells = []
        for name, spec in columns:
            if name == 'avgMaxCo' and r['diverged']:
                cell = 'diverged'
            else:
                cell = fmt(r.get(name), spec)
            cells.append(cell.rjust(widths[name]))
        print("  ".join(cells))

    if args.csv:
        with open(args.csv, 'w', newline='') as file:
            writer = csv.DictWriter(file, fieldnames=[name for name, _ in columns] + ['diverged'])
            writer.writeheader()
            for r in rows:
                writer.writerow({k: r.get(k) for k in writer.fieldnames})
        print(f"Wrote {args.csv}")


if __name__ == "__main__":
    main()
