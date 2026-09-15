#!/usr/bin/env python
#
# Plots relative wave height, (H - H1)/H1, as a function of time with one
# panel per solver (default: interFoam, interIsoFoam, interIsoGoam from top
# to bottom). H is crest height minus the height of the following trough, and
# H1 is the first wave height of each run (see plotCrestHeight.py).
#
# Within each panel there is one line per selected run of that solver. Runs
# with the same remaining parameters (e.g. nx, deltaT) get the same colour in
# all panels, so the solvers can be compared directly. The y-axes are shared.
#
# Cases are selected with the same --caseNames and [name=]pattern filters as
# in plotCrestHeight.py.
#
# Usage examples:
#   ./plotWaveHeightBySolver.py --scanDir gravityModelScan maxCo=0.1
#   ./plotWaveHeightBySolver.py --scanDir gravityModelScan_fixedDt nx=190,380 deltaT=0.005

import argparse

import matplotlib.pyplot as plt

from plotCrestHeight import (add_aligned_legend, add_common_arguments, legend_columns,
                             load_cases, relative, select_cases, show_or_save,
                             split_parameters)

parser = argparse.ArgumentParser(
    description='Plot relative wave height vs time, one panel per solver.',
    epilog='Example: ./plotWaveHeightBySolver.py --scanDir gravityModelScan_fixedDt '
           'nx=190,380 deltaT=0.005')
add_common_arguments(parser)
parser.add_argument('--solvers', nargs='+',
                    default=['interFoam', 'interIsoFoam', 'interIsoGoam'],
                    help='Solvers to plot, one panel each, from top to bottom. '
                         'Default: interFoam interIsoFoam interIsoGoam')
args = parser.parse_args()

selected = select_cases(args, parser)
cases = [c for c in load_cases(args, selected)
         if c['parms'].get('application') in args.solvers]
if not cases:
    parser.error("no cases with the requested solvers match the given filters")

# Shared parameters go in the title. Differing ones, except the solver which
# is given by the panel, label the lines and define the colours.
common, differing = split_parameters(cases)
differing = [p for p in differing if p != 'application']
common = [p for p in common if p != 'application']


def run_key(c):
    return tuple(c['parms'].get(p) for p in differing)


keys = []
for c in cases:
    if run_key(c) not in keys:
        keys.append(run_key(c))
colors = {k: 'C' + str(n % 10) for n, k in enumerate(keys)}

nPanels = len(args.solvers)
fig, ax = plt.subplots(nPanels, sharex=True, sharey=True, squeeze=False)
ax = ax[:, 0]
fig.set_size_inches(16, 5*nPanels)

columns = legend_columns(cases, differing)

for axis, solver in zip(ax, args.solvers):
    print(solver)
    handles, legendCases = [], []
    for c in cases:
        if c['parms'].get('application') != solver or len(c['H']) == 0:
            continue
        rel = relative(c['H'][:, 1])
        line, = axis.plot(c['H'][:, 0], rel, '-D', color=colors[run_key(c)])
        handles.append(line)
        legendCases.append(c)
        co = f", avg(maxCo) = {c['meanMaxCo']:.4g}" if c['meanMaxCo'] is not None else ""
        print(f"  {c['name']}: {len(rel)} waves, H1 = {c['H'][0, 1]:.6g}, last rel = {rel[-1]:.4g}{co}")

    axis.axhline(0, color='k', linewidth=1)
    axis.set_ylabel(solver + '\n' + r'$(H-H_1)/H_1$')
    axis.grid(True)
    if handles:
        # Widths from all cases, so the legends of all panels line up
        add_aligned_legend(axis, handles, legendCases, columns, allCases=cases,
                           fontsize=14, loc='lower left')
ax[-1].set_xlabel(r'$t$ [s]')

title = "Stream Function Wave"
if common:
    title += " (" + ", ".join(p + ": " + cases[0]['parms'][p] for p in common) + ")"
ax[0].set_title(title)

show_or_save(fig, args.figName)
