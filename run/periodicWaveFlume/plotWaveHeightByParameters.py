#!/usr/bin/env python
#
# Plots relative wave height, (H - H1)/H1, as a function of time with all
# solvers in the same panel and one panel per combination of the remaining
# parameters that differ between the selected runs (e.g. nx, deltaT or
# maxCo). H is crest height minus the height of the following trough, and H1
# is the first wave height of each run (see plotCrestHeight.py).
#
# Each solver has the same colour in all panels. The y-axes are shared.
#
# Cases are selected with the same --caseNames and [name=]pattern filters as
# in plotCrestHeight.py.
#
# Usage examples:
#   ./plotWaveHeightByParameters.py --scanDir gravityModelScan nx=380 maxCo=0.1,0.5
#   ./plotWaveHeightByParameters.py --scanDir gravityModelScan_fixedDt nx=190,380 deltaT=0.005

import argparse

import matplotlib.pyplot as plt

from plotCrestHeight import (add_aligned_legend, add_common_arguments, legend_columns,
                             load_cases, relative, select_cases, show_or_save,
                             split_parameters)

parser = argparse.ArgumentParser(
    description='Plot relative wave height vs time, all solvers in one panel, '
                'one panel per parameter combination.',
    epilog='Example: ./plotWaveHeightByParameters.py --scanDir gravityModelScan_fixedDt '
           'nx=190,380 deltaT=0.005')
add_common_arguments(parser)
parser.add_argument('--solvers', nargs='+',
                    default=['interFoam', 'interIsoFoam', 'interIsoGoam'],
                    help='Solvers to plot, in legend order. '
                         'Default: interFoam interIsoFoam interIsoGoam')
args = parser.parse_args()

selected = select_cases(args, parser)
cases = [c for c in load_cases(args, selected)
         if c['parms'].get('application') in args.solvers]
if not cases:
    parser.error("no cases with the requested solvers match the given filters")

# Shared parameters go in the figure title. Differing ones, except the solver,
# define the panels.
common, differing = split_parameters(cases)
panelParms = [p for p in differing if p != 'application']
common = [p for p in common if p != 'application']


def panel_key(c):
    return tuple(c['parms'].get(p) for p in panelParms)


panelKeys = []
for c in cases:
    if panel_key(c) not in panelKeys:
        panelKeys.append(panel_key(c))

colors = {s: 'C' + str(n % 10) for n, s in enumerate(args.solvers)}
columns = legend_columns(cases, ['application'])

nPanels = len(panelKeys)
fig, ax = plt.subplots(nPanels, sharex=True, sharey=True, squeeze=False)
ax = ax[:, 0]
fig.set_size_inches(16, 5*nPanels)

for axis, key in zip(ax, panelKeys):
    panelTitle = ", ".join(p + ": " + str(v) for p, v in zip(panelParms, key))
    print(panelTitle if panelTitle else "all cases")
    handles, legendCases = [], []
    # Loop over solvers first so lines and legend follow the --solvers order
    for solver in args.solvers:
        for c in cases:
            if (panel_key(c) != key or c['parms'].get('application') != solver
                    or len(c['H']) == 0):
                continue
            rel = relative(c['H'][:, 1])
            line, = axis.plot(c['H'][:, 0], rel, '-D', color=colors[solver])
            handles.append(line)
            legendCases.append(c)
            co = f", avg(maxCo) = {c['meanMaxCo']:.4g}" if c['meanMaxCo'] is not None else ""
            print(f"  {c['name']} {solver}: {len(rel)} waves, H1 = {c['H'][0, 1]:.6g}, "
                  f"last rel = {rel[-1]:.4g}{co}")

    axis.axhline(0, color='k', linewidth=1)
    axis.set_ylabel(r'$(H-H_1)/H_1$')
    axis.grid(True)
    if panelTitle:
        axis.set_title(panelTitle, fontsize=18)
    if handles:
        # Widths from all cases, so the legends of all panels line up
        add_aligned_legend(axis, handles, legendCases, columns, allCases=cases,
                           fontsize=14, loc='lower left')
ax[-1].set_xlabel(r'$t$ [s]')

title = "Stream Function Wave"
if common:
    title += " (" + ", ".join(p + ": " + cases[0]['parms'][p] for p in common) + ")"
fig.suptitle(title)

show_or_save(fig, args.figName)
