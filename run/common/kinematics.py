# Copyright (C) 2026 Johan Roenby
# SPDX-License-Identifier: GPL-3.0-or-later
"""Crest kinematics of a periodic stream-function wave against the exact solution.

Works for any uniform 2D blockMesh case of the periodic flume (one wavelength,
x periodic, z from the bed z1 to the top z2) whose constant/out holds the
stream-function coefficients (calcStreamFunctionWave output).

At each written time:
  surface  h(x) = z1 + sum(alpha dz) per column; crest and trough by parabolic
           interpolation; wave height H = crest - trough
  phase    crest position minus the exact crest x_c(0) + c t, as a fraction of L
  profile  the column nearest to the crest. The exact solution is evaluated at
           the column's offset from the solver's crest (theta = k (x_i - x_c)),
           i.e. the theoretical profile is matched to the numerical crest phase,
           as in Qwist & Christensen (2023, Sect. 4.2.3).
           Wet cells: alpha >= 0.5. e = |U| - |U_exact| per wet cell, scaled by
           the exact crest-surface speed U_ref = u(theta = 0, eta_crest).
           emax = max |e|/U_ref, emean = mean |e|/U_ref  (velocity magnitude
           profile errors in the spirit of their Fig. 16; their exact definition is
           not stated). r0 = u_x/u_x,exact in the top full water cell
           (alpha > 0.99); rTop = the same in the top wet cell.
Window averages are means over the snapshots in the window.

Usage: python3 kinematics.py <case> [<case> ...] [--win t0,t1 ...]
"""
import os, re, sys, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import foamio

G = 9.81

class Wave:
    def __init__(self, case):
        *v, B, E = foamio.entries(f'{case}/constant/out', 'H', 'd', 'k', 'omega', 'uBar', 'B', 'E')
        self.H, self.d, self.k, self.omega, self.uBar = map(float, v)
        self.c = self.omega/self.k; self.L = 2*np.pi/self.k
        self.B, self.E = foamio.as_array(B), foamio.as_array(E)
        self.j = np.arange(1, self.B.size + 1)

    def eta(self, theta):        # elevation above the bed
        jj = np.arange(self.E.size)
        return np.sum(2*self.E[:, None]*np.cos(np.outer(jj, np.atleast_1d(theta))), axis=0)

    def uv(self, theta, y):      # y above the bed
        th, y = np.broadcast_arrays(np.atleast_1d(theta), np.atleast_1d(y))
        jk = self.j[:, None]*self.k
        a = (self.j*self.B)[:, None]/np.cosh(jk*self.d)
        u = np.sqrt(G/self.k)*np.sum(a*np.cosh(jk*y)*np.cos(self.j[:, None]*th), axis=0) + self.omega/self.k - self.uBar
        v = np.sqrt(G/self.k)*np.sum(a*np.sinh(jk*y)*np.sin(self.j[:, None]*th), axis=0)
        return u, v

class Mesh:
    def __init__(self, case):
        v = foamio.entries(f'{case}/system/blockMeshDict', 'L', 'z1', 'z2', 'nx', 'nz')
        self.L, self.z1, self.z2 = map(float, v[:3])
        self.nx, self.nz = int(v[3]), int(v[4])
        self.dx, self.dz = self.L/self.nx, (self.z2 - self.z1)/self.nz

def _parab(h, i, sign):
    n = h.size; y0, y1, y2 = h[(i - 1) % n], h[i], h[(i + 1) % n]
    den = y0 - 2*y1 + y2; s = 0.5*(y0 - y2)/den if den != 0 else 0.0
    return i + 0.5 + s, y1 - 0.25*(y0 - y2)*s

def times(case):
    return sorted((d for d in os.listdir(case) if re.fullmatch(r'\d+(\.\d+)?(e-?\d+)?', d)), key=float)

def snapshot(case, t, M, W, xc0):
    A = foamio.read_vol(f'{case}/{t}/alpha.water').reshape(M.nz, M.nx)
    U = foamio.read_vol(f'{case}/{t}/U').reshape(M.nz, M.nx, 3)
    h = M.z1 + A.sum(0)*M.dz
    ic = int(np.argmax(h)); sc, eta_c = _parab(h, ic, 1)
    it = int(np.argmin(h)); _, eta_t = _parab(h, it, -1)
    xc = sc*M.dx
    out = dict(t=float(t), xc=xc, H=eta_c - eta_t)
    if xc0 is not None:
        dxc = (xc - (xc0 + W.c*float(t))) % M.L
        out['phase'] = (dxc - M.L if dxc > M.L/2 else dxc)/M.L
    xr = (ic + 0.5)*M.dx - xc; xr -= M.L*np.round(xr/M.L); th = W.k*xr
    z = M.z1 + (np.arange(M.nz) + 0.5)*M.dz
    wet = A[:, ic] >= 0.5
    ue, ve = W.uv(th, z[wet] - M.z1)
    Uref = W.uv(0.0, W.eta(0.0))[0][0]
    un, vn = U[wet, ic, 0], U[wet, ic, 2]
    e = np.hypot(un, vn) - np.hypot(ue, ve)
    out['emax'], out['emean'] = np.max(np.abs(e))/Uref, np.mean(np.abs(e))/Uref
    kW = np.where(A[:, ic] > 0.99)[0].max(); kT = np.where(wet)[0].max()
    uW = W.uv(th, z[kW] - M.z1)[0][0]; uT = W.uv(th, z[kT] - M.z1)[0][0]
    out['r0'], out['rTop'] = U[kW, ic, 0]/uW, U[kT, ic, 0]/uT
    out['prof'] = (z[wet], un, ue)
    return out

def analyse(case, windows):
    M, W = Mesh(case), Wave(case)
    ts = times(case)
    xc0 = None
    if ts and float(ts[0]) == 0:
        try: xc0 = snapshot(case, ts[0], M, W, None)['xc']
        except Exception: pass
    if xc0 is None:     # initial crest from streamFunFoam's phase: theta = k x - omega t + 2 pi phase
        for f in sorted(os.listdir(case)):
            if f.startswith('log.streamFunFoam'):
                m = re.search(r'Wave phase:\s*([-+0-9.eE]+)', open(f'{case}/{f}').read())
                if m: xc0 = (-float(m.group(1))*M.L) % M.L
    res = {}
    for (t0, t1) in windows:
        rows = []
        for t in ts:
            if not t0 < float(t) <= t1 + 1e-9: continue
            try: rows.append(snapshot(case, t, M, W, xc0))
            except Exception: continue
        if not rows: res[(t0, t1)] = None; continue
        agg = {k: np.mean([r[k] for r in rows]) for k in ('emax', 'emean', 'r0', 'rTop', 'H')}
        agg['phase'] = np.mean([r['phase'] for r in rows]) if xc0 is not None else np.nan
        agg['dH'] = agg['H']/W.H - 1; agg['n'] = len(rows); agg['last'] = rows[-1]
        res[(t0, t1)] = agg
    return res

def fmt(a):
    if a is None: return 'no snapshots'
    return (f"emax {100*a['emax']:6.2f}  emean {100*a['emean']:5.2f}  r0 {a['r0']:5.3f}  rTop {a['rTop']:5.3f}"
            f"  dH {100*a['dH']:+6.2f} %  phase {a['phase']:+.4f} L  (n={a['n']})")

if __name__ == '__main__':
    args = sys.argv[1:]; wins = []
    while '--win' in args:
        i = args.index('--win'); wins.append(tuple(float(x) for x in args[i + 1].split(','))); del args[i:i + 2]
    wins = wins or [(35, 40)]
    for case in args:
        for w, a in analyse(case, wins).items():
            print(f'{case}  {w[0]:g}-{w[1]:g} s  {fmt(a)}')
