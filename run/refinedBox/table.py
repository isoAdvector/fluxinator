"""Maximum velocity magnitude [m/s] after 10 time steps (t = 0.1 s), maximum over the
scanned water levels, in the layout of the non-orthogonal tables of the paper:
corrected schemes with zero initial p_rgh and 0 or 2 non-orthogonal correctors, with a
hydrostatic initial p_rgh, and uncorrected schemes. Each entry gives the maximum over
the 14 levels at which the interface cuts cells / over the three levels at which it
lies on mesh faces (z_w = 0.48, 0.50, 0.52; at 0.50 only in the refined region).
Usage: python3 table.py"""
import glob, os, re
import numpy as np
H = os.path.dirname(os.path.abspath(__file__))
ALIGNED = {'0.4800', '0.5000', '0.5200'}
ROWS = [('interFoam, basic', 'iF_corr', 'iF_nc2', 'iF_hyd', 'iF_unc'),
        ('interIsoFoam, basic', 'iIF_corr', 'iIF_nc2', 'iIF_hyd', 'iIF_unc'),
        ('interIsoGoam, constantHeight', 'cH_zero', 'cH_nc2', 'cH_hyd', 'cH_unc'),
        ('interIsoGoam, geometric, nL = 1', 'g1_zero', 'g1_nc2', 'g1_hyd', 'g1_unc'),
        ('interIsoGoam, geometric, nL = 2', 'g2_zero', 'g2_nc2', 'g2_hyd', 'g2_unc')]

def last(config):
    out = {}
    for d in glob.glob(f'{H}/scan/{config}_z*'):
        z = re.search(r'_z([0-9.]+)$', d).group(1)
        u = np.loadtxt(f'{d}/postProcessing/fieldMinMax1/0/fieldMinMax.dat', comments='#', usecols=(2,), ndmin=1)
        out[z] = u[-1]
    return out

print('max|U| [m/s] at t = 0.1 s, max over levels cutting cells / levels on faces')
print(f"{'':32s}{'zero, nc=0':>22s}{'zero, nc=2':>22s}{'hydrostatic, nc=0':>22s}{'uncorrected':>22s}")
for label, *configs in ROWS:
    cells = []
    for c in configs:
        d = last(c)
        cut = max(v for z, v in d.items() if z not in ALIGNED)
        on = max(v for z, v in d.items() if z in ALIGNED)
        cells.append(f'{cut:.1e} / {on:.1e}')
    print(f'{label:32s}' + ''.join(f'{x:>22s}' for x in cells))
