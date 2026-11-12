# Copyright (C) 2026 Johan Roenby
# SPDX-License-Identifier: GPL-3.0-or-later
"""Maximum velocity magnitude [m/s] after 10 time steps (t = 0.1 s) at the water level
z_w = 0.49 m, in the layout of the non-orthogonal tables of the paper: corrected schemes with
zero initial p_rgh and 0 or 2 non-orthogonal correctors, with a hydrostatic initial p_rgh,
and uncorrected schemes. Below it, the range over the 15 scanned levels at which the
interface cuts cells (all but z_w = 0.48 and 0.52, where it lies on coarse cell faces; at
0.50 it lies on fine faces in the refined region but cuts the coarse cells outside it).
Part I, Table tab:nonOrthRefined. Input: run/refinedBox/waterLevel_10steps/<config>_z<zw>.
Output: tab_nonOrthRefined.txt. Usage: python3 make.py"""
import glob, os, re, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..'))
from common import case, write_table  # noqa: E402
H = case('refinedBox')
LEVEL = '0.4900'
ON_COARSE_FACES = {'0.4800', '0.5200'}
ROWS = [('interFoam, basic', 'iF_corr', 'iF_nc2', 'iF_hyd', 'iF_unc'),
        ('interIsoFoam, basic', 'iIF_corr', 'iIF_nc2', 'iIF_hyd', 'iIF_unc'),
        ('interIsoGoam, constantHeight', 'cH_zero', 'cH_nc2', 'cH_hyd', 'cH_unc'),
        ('interIsoGoam, geometric, nL = 1', 'g1_zero', 'g1_nc2', 'g1_hyd', 'g1_unc'),
        ('interIsoGoam, geometric, nL = 2', 'g2_zero', 'g2_nc2', 'g2_hyd', 'g2_unc')]

def last(config):
    out = {}
    for d in glob.glob(f'{H}/waterLevel_10steps/{config}_z*'):
        z = re.search(r'_z([0-9.]+)$', d).group(1)
        u = np.loadtxt(f'{d}/postProcessing/fieldMinMax1/0/fieldMinMax.dat', comments='#', usecols=(2,), ndmin=1)
        out[z] = u[-1]
    return out

HEAD = f"{'':32s}{'zero, nc=0':>22s}{'zero, nc=2':>22s}{'hydrostatic, nc=0':>22s}{'uncorrected':>22s}"
out = [f'max|U| [m/s] at t = 0.1 s, z_w = {float(LEVEL):.2f} m', HEAD]
ranges = ['', 'min - max over the 15 levels at which the interface cuts cells', HEAD]
for label, *configs in ROWS:
    at, span = [], []
    for c in configs:
        d = last(c)
        cut = [v for z, v in d.items() if z not in ON_COARSE_FACES]
        at.append(f'{d[LEVEL]:.1e}')
        span.append(f'{min(cut):.1e} - {max(cut):.1e}')
    out.append(f'{label:32s}' + ''.join(f'{x:>22s}' for x in at))
    ranges.append(f'{label:32s}' + ''.join(f'{x:>22s}' for x in span))
write_table(HERE, 'tab_nonOrthRefined', out + ranges)
