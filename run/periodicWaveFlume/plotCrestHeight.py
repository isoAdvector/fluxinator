#!/usr/bin/env python
#
# Plots relative crest and trough heights, (h - h1)/h1, as functions of time,
# where h1 is the height of the first crest (or trough). Heights are taken
# from the hL column (elevation relative to the wave gauge location) in
# postProcessing/interfaceHeight1/0/height.dat.
#
# Crests and troughs are shown in separate panels sharing the time axis.
# Optionally (--waveHeight) a third panel shows the relative wave height,
# (H - H1)/H1, where H = crest height - height of the following trough.
#
# Cases are selected either by name (--caseNames) or by [name=]pattern filters
# matched against each case's parmScannerParms (or both).
#
# Usage examples:
#   ./plotCrestHeight.py --scanDir gravityModelScan --caseNames 0000 0001 0002
#   ./plotCrestHeight.py --scanDir gravityModelScan nx=760 "application=inter*Foam" maxCo=0.2
#
# The functions in this file are also used by plotWaveHeightBySolver.py.

import argparse
import fnmatch
import os
import re

import matplotlib.pyplot as plt
import numpy as np

plt.rcParams.update({
    "font.family": "monospace",
    "font.monospace": "DejaVu Sans Mono",
    "font.size": 20,
    'lines.linewidth': 1.5
})


def refine_extremum(t, h, i):
    # Parabola through the three samples around index i (non-uniform times)
    # to get sub-sample estimates of extremum time and height.
    if i == 0 or i == len(h) - 1:
        return t[i], h[i]
    tt = t[i-1:i+2] - t[i]
    a, b, c = np.polyfit(tt, h[i-1:i+2], 2)
    if a == 0:
        return t[i], h[i]
    ts = -b/(2*a)
    if abs(ts) > max(abs(tt[0]), abs(tt[2])):
        return t[i], h[i]
    return t[i] + ts, c - b*b/(4*a)


def find_extrema(t, h, sign, threshold):
    # sign = +1 for crests, -1 for troughs.
    # The signal is split into segments where sign*(h - mean) > 0, and the
    # extremum of each segment is taken. Segments whose extremum lies at the
    # ends of the record (incomplete crests/troughs) or whose amplitude is
    # below threshold (noise around the mean level) are skipped.
    d = sign*(h - np.mean(h))
    idx = []
    i = 0
    n = len(d)
    while i < n:
        if d[i] > 0:
            j = i
            while j < n and d[j] > 0:
                j += 1
            k = i + np.argmax(d[i:j])
            if 0 < k < n - 1 and d[k] > threshold:
                idx.append(k)
            i = j
        else:
            i += 1
    ext = np.array([refine_extremum(t, h, k) for k in idx])
    return ext.reshape(-1, 2)


def wave_heights(crests, troughs):
    # Pair each crest with the first trough following it (before the next
    # crest). Returns (mean time of crest and trough, crest - trough).
    H = []
    for n, (tc, hc) in enumerate(crests):
        tNext = crests[n+1, 0] if n + 1 < len(crests) else np.inf
        later = troughs[(troughs[:, 0] > tc) & (troughs[:, 0] < tNext)]
        if len(later) > 0:
            tt, ht = later[0]
            H.append((0.5*(tc + tt), hc - ht))
    return np.array(H).reshape(-1, 2)


def read_parms(caseDir):
    # Returns scanned parameters as an ordered dict {name: value}
    parmFile = os.path.join(caseDir, "parmScannerParms")
    if not os.path.isfile(parmFile):
        return {}
    with open(parmFile, 'r') as file:
        parms = file.read().strip()
    parms = parms.split(':', 1)[1]  # Strip leading "scanDir/case:"
    return dict(re.findall(r'([A-Za-z]\w*):\s*(\S+)', parms))


def read_log(caseDir, application):
    # Returns (last ClockTime, time average of max Courant number) from
    # log.<application>, with None for values that are not available.
    # The max Courant number of a time step is taken from the last
    # "Courant Number mean: ... max: ..." line before its "Time = " line.
    logFile = os.path.join(caseDir, "log." + str(application))
    if application is None or not os.path.isfile(logFile):
        return None, None
    clockTime = None
    maxCo = []
    lastMaxCo = None
    with open(logFile, 'r') as file:
        for line in file:
            if line.startswith('Courant Number'):
                match = re.search(r'max:\s*(\S+)', line)
                if match:
                    lastMaxCo = float(match.group(1))
            elif line.startswith('Time = '):
                if lastMaxCo is not None:
                    maxCo.append(lastMaxCo)
                lastMaxCo = None
            elif 'ClockTime' in line:
                match = re.search(r'ClockTime = (\d+\.?\d*) s', line)
                if match:
                    clockTime = float(match.group(1))
    return clockTime, (np.mean(maxCo) if maxCo else None)


def value_matches(value, pattern):
    # pattern may be a comma separated list of alternatives, each of which is
    # a shell-style wildcard pattern or a number (compared numerically)
    for p in pattern.split(','):
        try:
            if float(value) == float(p):
                return True
        except ValueError:
            pass
        if fnmatch.fnmatchcase(value, p):
            return True
    return False


def filter_matches(parms, f):
    # "name=pattern" matches that parameter only; a bare "pattern" matches
    # if any parameter value matches
    if '=' in f:
        name, pattern = f.split('=', 1)
        return name in parms and value_matches(parms[name], pattern)
    return any(value_matches(v, f) for v in parms.values())


def relative(y):
    return (y - y[0])/y[0]


def add_common_arguments(parser):
    parser.add_argument('filters', nargs='*', metavar='[name=]pattern',
                        help='Only plot cases whose parmScannerParms match all filters. '
                             'name=pattern matches the named parameter (e.g. maxCo, deltaT); '
                             'a bare pattern matches any parameter value. '
                             'Patterns may contain shell wildcards (*, ?) and comma '
                             'separated alternatives. Numbers are compared numerically.')
    parser.add_argument('--scanDir', required=True, help='Directory where cases are found')
    parser.add_argument('--caseNames', nargs='+', default=None,
                        help='List of case names (default: all 4-digit case directories in scanDir)')
    parser.add_argument('--figName', default=None,
                        help='Save figure as figName.png instead of showing it')
    parser.add_argument('--tMin', type=float, default=-1, help='Start of time window')
    parser.add_argument('--tMax', type=float, default=1e30, help='End of time window')
    parser.add_argument('--column', type=int, default=2,
                        help='Column in height.dat to use (1: hB, 2: hL). Default: 2')
    parser.add_argument('--threshold', type=float, default=0.1,
                        help='Minimum |h - mean| of an extremum as fraction of signal range. Default: 0.1')


def select_cases(args, parser):
    caseNames = args.caseNames
    if caseNames is None:
        # Only regular scan cases named by 4 digits, e.g. 0026 (not 0026_special)
        caseNames = sorted(d for d in os.listdir(args.scanDir)
                           if re.fullmatch(r'\d{4}', d)
                           and os.path.isfile(os.path.join(args.scanDir, d, 'parmScannerParms')))
    selected = [caseName for caseName in caseNames
                if all(filter_matches(read_parms(os.path.join(args.scanDir, caseName)), f)
                       for f in args.filters)]
    if not selected:
        parser.error("no cases match the given filters")
    return selected


def load_cases(args, caseNames):
    # Read data and detect extrema for the given cases
    cases = []
    for caseName in caseNames:
        caseDir = os.path.join(args.scanDir, caseName)
        t, h = np.loadtxt(os.path.join(caseDir, 'postProcessing/interfaceHeight1/0/height.dat'),
                          comments="#", unpack=True, usecols=(0, args.column))
        mask = (t >= args.tMin) & (t <= args.tMax)
        t, h = t[mask], h[mask]
        threshold = args.threshold*(h.max() - h.min())
        crests = find_extrema(t, h, +1, threshold)
        troughs = find_extrema(t, h, -1, threshold)
        parms = read_parms(caseDir)
        clockTime, meanMaxCo = read_log(caseDir, parms.get('application'))
        cases.append(dict(name=caseName, parms=parms, crests=crests,
                          troughs=troughs, H=wave_heights(crests, troughs),
                          clockTime=clockTime, meanMaxCo=meanMaxCo))
    return cases


def legend_columns(cases, parameters):
    # Legend columns: case name, the given parameters, the time averaged max
    # Courant number (only for fixed deltaT runs, placed after deltaT) and
    # the clock time
    columns = ['case'] + list(parameters)
    if any('deltaT' in c['parms'] for c in cases):
        i = columns.index('deltaT') + 1 if 'deltaT' in columns else len(columns)
        columns.insert(i, 'Co')
    columns.append('time')
    return columns


def legend_cell(c, column):
    if column == 'case':
        return c['name']
    if column == 'Co':
        co = c['meanMaxCo']
        if co is None:
            return '-'
        # A diverging run gives huge (or inf/nan) Courant numbers
        if not np.isfinite(co) or co > 1e3:
            return 'diverged'
        return f"{co:.3f}" if co < 10 else f"{co:.3g}"
    if column == 'time':
        return f"{c['clockTime']:.0f} s" if c['clockTime'] is not None else '-'
    return c['parms'].get(column, '-')


def add_aligned_legend(axis, handles, cases, columns, allCases=None, **kwargs):
    # Legend with one row per case and columns padded to equal width (the font
    # is monospace). The first row is a header naming the columns. Column
    # widths are computed from allCases (default: cases), so legends in
    # different panels can be aligned with each other.
    headers = {'Co': 'avg(maxCo)', 'time': 'clockTime'}
    header = [headers.get(col, col) for col in columns]
    widthRows = [header] + [[legend_cell(c, col) for col in columns]
                            for c in (allCases if allCases is not None else cases)]
    widths = [max(len(row[i]) for row in widthRows) for i in range(len(columns))]

    def fmt(row):
        return "  ".join(v.ljust(w) for v, w in zip(row, widths)).rstrip()

    labels = [fmt([legend_cell(c, col) for col in columns]) for c in cases]
    headerHandle = plt.Line2D([], [], linestyle='none')
    return axis.legend([headerHandle] + list(handles), [fmt(header)] + labels, **kwargs)


def split_parameters(cases):
    # Returns (parameters shared by all cases, parameters that differ)
    allNames = []
    for c in cases:
        allNames += [p for p in c['parms'] if p not in allNames]
    common = [p for p in allNames
              if len(set(c['parms'].get(p) for c in cases)) == 1]
    differing = [p for p in allNames if p not in common]
    return common, differing


def show_or_save(fig, figName):
    fig.tight_layout()
    if figName:
        plt.savefig(figName + '.png', bbox_inches='tight', dpi=300)
    else:
        plt.show()


def main():
    parser = argparse.ArgumentParser(
        description='Plot relative crest and trough heights vs time.',
        epilog='Example: ./plotCrestHeight.py --scanDir gravityModelScan nx=760 '
               '"application=inter*Foam" maxCo=0.1,0.2')
    add_common_arguments(parser)
    parser.add_argument('--waveHeight', action='store_true',
                        help='Add a panel with relative wave height (crest - trough)')
    args = parser.parse_args()

    cases = load_cases(args, select_cases(args, parser))

    # Parameters shared by all cases go in the title, differing ones in the legend
    common, differing = split_parameters(cases)

    nPanels = 3 if args.waveHeight else 2
    fig, ax = plt.subplots(nPanels, sharex=True)
    fig.set_size_inches(16, 5*nPanels)

    handles, legendCases = [], []
    for count, c in enumerate(cases):
        color = 'C' + str(count % 10)
        print(c['name'] + (f", avg(maxCo) = {c['meanMaxCo']:.4g}" if c['meanMaxCo'] is not None else ""))
        for axis, key, marker, what in [(ax[0], 'crests', 'o', 'crest'),
                                        (ax[1], 'troughs', 's', 'trough')] + \
                ([(ax[2], 'H', 'D', 'H')] if args.waveHeight else []):
            data = c[key]
            if len(data) == 0:
                continue
            rel = relative(data[:, 1])
            line, = axis.plot(data[:, 0], rel, '-', marker=marker, color=color)
            if axis is ax[0]:
                handles.append(line)
                legendCases.append(c)
            print(f"  {what:7s}: {len(data)}, first = {data[0, 1]:.6g}, last rel = {rel[-1]:.4g}")

    ylabels = [r'crest $(h-h_1)/h_1$', r'trough $(h-h_1)/h_1$', r'$(H-H_1)/H_1$']
    for i in range(nPanels):
        ax[i].axhline(0, color='k', linewidth=1)
        ax[i].set_ylabel(ylabels[i])
        ax[i].grid(True)
    ax[-1].set_xlabel(r'$t$ [s]')

    title = "Stream Function Wave"
    if common:
        title += " (" + ", ".join(p + ": " + cases[0]['parms'][p] for p in common) + ")"
    ax[0].set_title(title)
    add_aligned_legend(ax[0], handles, legendCases, legend_columns(cases, differing),
                       fontsize=14, loc='best')

    show_or_save(fig, args.figName)


if __name__ == "__main__":
    main()
