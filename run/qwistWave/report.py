"""Tables for the Qwist & Christensen wave runs (runs/<variant>/m<k>_dt<dt>).

Stability: last time reached, max Courant number at t = 0, median and max over
the run (log). Kinematics (kinematics.py) in two windows:
  16-20 s  ~ their measuring point: a crest in 10L-11L of their flume has
           propagated ~9 L (18 s) from the end of the inlet relaxation zone
  35-40 s  their last 5 s of 40 s (here 17.5-20 L of propagation)
CFL_Q = their CFL definition: 0.09 at their dt on each mesh, scaling with dt.
Usage: python3 report.py  -> results.txt"""
import os, re, sys, numpy as np
from multiprocessing import Pool
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import kinematics as K
VARIANTS = ['iIF_E_on', 'iIF_E_off', 'iIF_B_on', 'iIF_B_off', 'iIG_E_on', 'iIG_E_off', 'iIG_B_on', 'iIG_B_off']
DT0 = {1: 0.02, 2: 0.01, 3: 0.005, 4: 0.0025, 5: 0.002}
HDY = {1: 1.5, 2: 3, 3: 6, 4: 12, 5: 15}
WINS = [(16, 20), (35, 40)]

def runs():
    out = []
    for v in VARIANTS:
        d = os.path.join(HERE, 'runs', v)
        if not os.path.isdir(d): continue
        for r in os.listdir(d):
            m = re.fullmatch(r'm(\d)_dt([0-9.e-]+)', r)
            if m: out.append((v, int(m.group(1)), float(m.group(2)), os.path.join(d, r)))
    return out

def logstats(case):
    try: txt = open(os.path.join(case, 'log.interIsoGoam')).read()
    except OSError: return None
    co = [float(x) for x in re.findall(r'^Courant Number mean: \S+ max: (\S+)', txt, re.M)]
    ts = re.findall(r'^Time = (\S+)', txt, re.M)
    return dict(tLast=float(ts[-1]) if ts else 0.0, end='\nEnd' in txt, co0=co[0] if co else np.nan,
                coMed=np.median(co) if co else np.nan, coMax=np.max(co) if co else np.nan)

def job(r):
    v, m, dt, case = r
    s = logstats(case)
    try: k = K.analyse(case, WINS)
    except Exception: k = {w: None for w in WINS}
    return r, s, k

def cell(k, key, sc=100, fm='{:6.2f}'):
    return ' ' * 6 if k is None else fm.format(sc*k[key])

if __name__ == '__main__':
    R = runs()
    with Pool(32) as p: res = p.map(job, R)
    by = {(v, m, dt): (s, k) for (v, m, dt, _), s, k in res}
    L = [__doc__.split('Usage')[0].strip(), '',
         'emax/emean: max/mean | |U| - |U_exact| | over wet crest-column cells, % of the exact crest-surface',
         'speed (1.04 m/s); r0 = u_x/u_x,exact in the top full water cell; dH wave height error [%];',
         'phase = crest position error [L]. "x@t" = run stopped (diverged) at time t.', '']
    def status(s):
        if s is None: return 'missing'
        return 'ok' if s['end'] else f"x@{s['tLast']:.1f}"
    # A: mesh convergence at their dt
    for w in WINS:
        L.append(f'== A: their meshes and dt (CFL_Q 0.09), window {w[0]}-{w[1]} s')
        L.append(f'{"variant":11s}{"H/dy":>5s} {"status":>8s} {"Co0":>6s} {"CoMed":>6s} {"CoMax":>8s}  {"emax":>6s} {"emean":>6s} {"r0":>6s} {"dH":>7s} {"phase":>8s}')
        for v in VARIANTS:
            for m in DT0:
                if (v, m, DT0[m]) not in by: continue
                s, k = by[(v, m, DT0[m])]; a = k.get(w) if k else None
                L.append(f'{v:11s}{HDY[m]:5g} {status(s):>8s} {s["co0"]:6.3f} {s["coMed"]:6.3f} {s["coMax"]:8.3g}  '
                         f'{cell(a,"emax")} {cell(a,"emean")} {cell(a,"r0",1,"{:6.3f}")} {cell(a,"dH",100,"{:+7.2f}")} {cell(a,"phase",1,"{:+8.4f}")}')
            L.append('')
    # B: dt scan
    for m in (3, 4):
        for w in WINS:
            L.append(f'== B: dt scan on m{m} (H/dy {HDY[m]:g}), window {w[0]}-{w[1]} s')
            L.append(f'{"variant":11s}{"CFL_Q":>6s} {"dt":>8s} {"status":>8s} {"Co0":>6s} {"CoMed":>6s} {"CoMax":>8s}  {"emax":>6s} {"emean":>6s} {"r0":>6s} {"dH":>7s} {"phase":>8s}')
            for v in VARIANTS:
                for dt in sorted(dt for (vv, mm, dt) in by if vv == v and mm == m):
                    s, k = by[(v, m, dt)]; a = k.get(w) if k else None
                    L.append(f'{v:11s}{0.09*dt/DT0[m]:6.3f} {dt:8g} {status(s):>8s} {s["co0"]:6.3f} {s["coMed"]:6.3f} {s["coMax"]:8.3g}  '
                             f'{cell(a,"emax")} {cell(a,"emean")} {cell(a,"r0",1,"{:6.3f}")} {cell(a,"dH",100,"{:+7.2f}")} {cell(a,"phase",1,"{:+8.4f}")}')
                L.append('')
    open(os.path.join(HERE, 'results.txt'), 'w').write('\n'.join(L) + '\n')
    print('\n'.join(L))
