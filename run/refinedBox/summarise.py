"""max|U| after the first and the tenth step, and the first p_rgh residual, per
configuration (rows) and water level z_w (columns), for scan/<cfg>_z<zw>
(configurations: makeCases.sh). Usage: python3 summarise.py"""
import glob, os, re, sys
H = os.path.dirname(os.path.abspath(__file__))
CFG = ['iF_corr', 'iF_nc2', 'iF_hyd', 'iF_unc', 'iIF_corr', 'iIF_nc2', 'iIF_hyd', 'iIF_unc',
       'cH_zero', 'cH_nc2', 'cH_hyd', 'cH_unc', 'g1_zero', 'g1_nc2', 'g1_hyd', 'g1_unc',
       'g2_zero', 'g2_nc2', 'g2_hyd', 'g2_unc', 'g3_hyd']
def read(c):
    f = f'{c}/postProcessing/fieldMinMax1/0/fieldMinMax.dat'
    rows = [l.split() for l in open(f) if not l.startswith('#') and l.strip()] if os.path.exists(f) else []
    app = 'interFoam' if '/iF_' in c else 'interIsoGoam'
    L = open(f'{c}/log.{app}').read() if os.path.exists(f'{c}/log.{app}') else ''
    m = re.search(r'Solving for p_rgh, Initial residual = ([^,]+)', L)
    return [float(r[2]) for r in rows], (float(m.group(1)) if m else float('nan'))
Z = sorted({re.search(r'_z([0-9.]+)$', d).group(1) for d in glob.glob(f'{H}/scan/*_z*')})
data = {(c, z): read(f'{H}/scan/{c}_z{z}') for c in CFG for z in Z}
for title, k in (('max|U| after step 1 [m/s]', 0), ('max|U| after step 10 [m/s]', -1)):
    print(f'\n{title}\n{"z_w":9s}' + ''.join(f'{z[1:]:>8s}' for z in Z))
    for c in CFG:
        print(f'{c:9s}' + ''.join(f'{data[(c, z)][0][k]:8.0e}' if data[(c, z)][0] else '     n/a' for z in Z))
print(f'\nfirst p_rgh initial residual (hydrostatic starts)\n{"z_w":9s}' + ''.join(f'{z[1:]:>8s}' for z in Z))
for c in [c for c in CFG if 'hyd' in c]:
    print(f'{c:9s}' + ''.join(f'{data[(c, z)][1]:8.0e}' for z in Z))
