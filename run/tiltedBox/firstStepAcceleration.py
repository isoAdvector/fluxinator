"""First-step spurious acceleration max|U|/(dt g) of the standard (basic) gravity term in the
tilted box, per mesh, from paperData/meshDeltaTScan.txt. Writes paperData/firstStepAcceleration.txt.
Usage (in run/tiltedBox): python3 firstStepAcceleration.py"""
import re
txt = open('paperData/meshDeltaTScan.txt').read().split('\n')
rows = {}
for i, l in enumerate(txt):
    m = re.match(r'\d{4}\s+nx: (\d+).*application: (\S+) deltaT: (\S+)', l)
    if m:
        u = float(re.search(r'max\|U\| t=\S+: (\S+)', txt[i + 1]).group(1))
        rows.setdefault(int(m.group(1)), []).append((m.group(2), float(m.group(3)), u))
out = ['First time step in the tilted box, standard (basic) gravity term: max|U|/(dt g), g = 9.81 m/s^2',
       'Source: run/tiltedBox/meshDeltaTScan (paperData/meshDeltaTScan.txt), script firstStepAcceleration.py',
       f'{"mesh":>10s} {"interFoam":>22s} {"interIsoFoam (basic)":>22s}   (min ... max over dt = 1e-4 ... 1 s)']
for n in sorted(rows):
    s = []
    for app in ('interFoam', 'interIsoGoam'):
        a = [u/dt/9.81 for ap, dt, u in rows[n] if ap == app]
        s.append(f'{min(a):.4f} ... {max(a):.4f}')
    out.append(f'{n}x{n:<7d} {s[0]:>22s} {s[1]:>22s}')
open('paperData/firstStepAcceleration.txt', 'w').write('\n'.join(out) + '\n')
print('\n'.join(out))
