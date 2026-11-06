# Copyright (C) 2026 Johan Roenby
# SPDX-License-Identifier: GPL-3.0-or-later
"""Height-matched crest kinematics: compare each run with the exact stream-function
wave of the run's own mean wave height in the window (same T, calcStreamFunctionWave),
so that wave damping does not masquerade as a smaller velocity overshoot.
r0 = u_x/u_x,exact in the crest's top full water cell, emax = max | |U|-|U_exact| |
over the wet crest-column cells [% of the exact crest-surface speed of that wave].
Celerity: least-squares slope of the unwrapped crest position over all snapshots with
t > 2 s (first period excluded), relative to the exact celerity of the initial wave (dc)
and of the stream-function wave with the run's mean height over that time (dcH).
Usage: python3 hmatch.py <case> ... [--win t0,t1 ...] [--cel]"""
import os, re, sys, subprocess, tempfile, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import kinematics as K
TMP = tempfile.mkdtemp(prefix='hmatch_')

def wave_of_height(case, H):
    d = os.path.join(TMP, f'{abs(hash((case, round(H, 7))))}'); os.makedirs(d + '/constant', exist_ok=True)
    wp = re.sub(r'\n\s*H\s+[^;]*;', f'\n    H {H:.7f};', open(case + '/constant/waveProperties').read())
    open(d + '/constant/waveProperties', 'w').write(wp)
    env = ""
    subprocess.run(['bash', '-c', env + f'calcStreamFunctionWave {d}/constant/waveProperties {d}/constant/out'], capture_output=True)
    return K.Wave(d)

def analyse(case, win):
    M, W0 = K.Mesh(case), K.Wave(case)
    ts = [t for t in K.times(case) if win[0] < float(t) <= win[1] + 1e-9]
    s0 = [K.snapshot(case, t, M, W0, None) for t in ts]
    if not s0: return None
    H = np.mean([s['H'] for s in s0]); W = wave_of_height(case, H)
    s = [K.snapshot(case, t, M, W, None) for t in ts]
    return dict(dH=H/W0.H - 1, r0full=np.mean([x['r0'] for x in s0]), r0=np.mean([x['r0'] for x in s]),
                emax=np.mean([x['emax'] for x in s]), n=len(s))

def analyse_input(case, win):
    """As analyse, but against the input (initial) wave: no height matching."""
    M, W0 = K.Mesh(case), K.Wave(case)
    ts = [t for t in K.times(case) if win[0] < float(t) <= win[1] + 1e-9]
    s = [K.snapshot(case, t, M, W0, None) for t in ts]
    if not s: return None
    return dict(dH=np.mean([x['H'] for x in s])/W0.H - 1, r0=np.mean([x['r0'] for x in s]),
                emax=np.mean([x['emax'] for x in s]), n=len(s))

def celerity(case, t0=2.0):
    M, W0 = K.Mesh(case), K.Wave(case)
    ts = [t for t in K.times(case) if float(t) > t0]
    s = [K.snapshot(case, t, M, W0, None) for t in ts]
    if len(s) < 3: return None
    tt = np.array([float(t) for t in ts]); xc = np.unwrap(np.array([x['xc'] for x in s])*2*np.pi/M.L)*M.L/(2*np.pi)
    c = np.polyfit(tt, xc, 1)[0]; H = np.mean([x['H'] for x in s]); W = wave_of_height(case, H)
    return dict(c=c, dc=c/W0.c - 1, dcH=c/W.c - 1, H=H/W0.H - 1)

if __name__ == '__main__':
    args = sys.argv[1:]; wins = []
    while '--win' in args:
        i = args.index('--win'); wins.append(tuple(float(x) for x in args[i + 1].split(','))); del args[i:i + 2]
    wins = wins or [(8, 10)]
    cel = '--cel' in args
    if cel: args.remove('--cel')
    for c in args:
        if cel:
            a = celerity(c)
            print(f'{c:44s} celerity ' + ('n/a' if a is None else
                  f'{a["c"]:.4f} m/s  vs initial wave {100*a["dc"]:+.2f} %  vs same-height wave {100*a["dcH"]:+.2f} %  (mean dH {100*a["H"]:+.1f} %)'))
        for w in wins:
            a = analyse(c, w)
            print(f'{c:44s} {w[0]:g}-{w[1]:g} s  ' + ('no snapshots' if a is None else
                  f'dH {100*a["dH"]:+6.1f} %   r0 vs full H {a["r0full"]:.3f}   height-matched: r0 {a["r0"]:.3f}  emax {100*a["emax"]:5.1f} %  (n={a["n"]})'))
