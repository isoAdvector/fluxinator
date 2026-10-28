"""Jesper's wave (H 0.3, T 2, d 1), mesh m3 (108x40, 6 cells per H): velocity profiles
under the crest and wave shape at t = 1T, 2T, 5T, 10T, 20T for implicit convection
(today's interIsoGoam), explicit + Lax-Wendroff, explicit + interfaceUpwind, against the
exact stream-function wave of the initial height. Rows: his dt (0.005) and 2x dt (0.01).
The domain is one wavelength and periodic, so the exact wave at t = nT is the initial
wave at its initial position: any offset is a phase (celerity) error.
Output figures/jesper_profiles.png, figures/jesper_shape.png"""
import os, numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
import kinematics as K
INK, INK2, GRID = '#1f1f1e', '#5f5e5a', '#e4e2dc'
SOL = [('implicit (today)', 'iIG_E_off', '#eb6834'), ('explicit + Lax–Wendroff', 'LW', '#2a78d6'),
       ('explicit + interfaceUpwind', 'IU', '#1baf7a')]
ROWS = [('m3_dt0.005', 'Δt = 0.005 s (his)'), ('m3_dt0.01', 'Δt = 0.01 s (2×)')]
NT = [1, 2, 5, 10, 20]

def crest_x(c, M, t):
    A = K.foamio.read_vol(f'{c}/{t}/alpha.water').reshape(M.nz, M.nx)
    h = A.sum(0)*M.dz; s, _ = K._parab(h, int(np.argmax(h)), 1); return s*M.dx, M.z1 + h

def style(ax):
    ax.grid(True, color=GRID, lw=0.6); [s.set_visible(False) for s in (ax.spines['top'], ax.spines['right'])]

# --- profiles under the crest
fig, axs = plt.subplots(len(ROWS), len(NT), figsize=(16, 9), sharey=True)
for r, (m, lab) in enumerate(ROWS):
    for jc, n in enumerate(NT):
        ax = axs[r, jc]; t = f'{2*n:g}'
        for name, d, col in SOL:
            c = f'runs/{d}/{m}'; M, W = K.Mesh(c), K.Wave(c)
            s = K.snapshot(c, t, M, W, None); z, un, _ = s['prof']
            ax.plot(un, z, 'o-', color=col, lw=1.6, ms=3, label=f'{name}  (H {100*(s["H"]/W.H - 1):+.0f} %)')
        ze = np.linspace(M.z1, M.z1 + W.eta(0.0)[0], 200)
        ax.plot(W.uv(0.0, ze - M.z1)[0], ze, '--', color=INK, lw=1.4, label='exact (initial wave)')
        ax.set_title(f'{lab},  t = {n}T = {2*n} s', fontsize=9.5, color=INK); style(ax)
        ax.set_ylim(-1.0, 0.3); ax.legend(fontsize=6.5, frameon=False, loc='lower right')
        if jc == 0: ax.set_ylabel('z above still water level [m]')
        if r == len(ROWS) - 1: ax.set_xlabel('u under the crest [m/s]')
fig.suptitle("Jesper's wave on m3 (6 cells per wave height): horizontal velocity under the crest", fontsize=11, color=INK)
fig.tight_layout(); fig.savefig('figures/jesper_profiles.png', dpi=120)

# --- wave shape
fig, axs = plt.subplots(len(ROWS), len(NT), figsize=(16, 8), sharey=True, sharex=True)
for r, (m, lab) in enumerate(ROWS):
    c0 = f'runs/{SOL[0][1]}/{m}'; M, W = K.Mesh(c0), K.Wave(c0)
    x = (np.arange(M.nx) + 0.5)*M.dx
    x0, _ = crest_x(f'runs/IU/{m}', M, '0')
    xr = (x - x0 + M.L/2) % M.L - M.L/2; o = np.argsort(xr)
    xe = np.linspace(-M.L/2, M.L/2, 400)
    for jc, n in enumerate(NT):
        ax = axs[r, jc]; t = f'{2*n:g}'
        ax.plot(xe, W.eta(W.k*xe) - W.d, '--', color=INK, lw=1.4, label='exact (initial wave)')
        for name, d, col in SOL:
            c = f'runs/{d}/{m}'; _, eta = crest_x(c, M, t)
            ax.plot(xr[o], eta[o], '-', color=col, lw=1.6, label=name)
        ax.set_title(f'{lab},  t = {n}T', fontsize=9.5, color=INK); style(ax)
        if jc == 0: ax.set_ylabel('η [m]')
        if r == len(ROWS) - 1: ax.set_xlabel('x − initial crest position [m]')
        ins = ax.inset_axes([0.64, 0.56, 0.35, 0.42])
        ins.plot(xe, W.eta(W.k*xe) - W.d, '--', color=INK, lw=1.2)
        for name, d, col in SOL:
            _, eta = crest_x(f'runs/{d}/{m}', M, t); ins.plot(xr[o], eta[o], '-', color=col, lw=1.4, marker='.', ms=3)
        ins.set_xlim(-0.8, 0.8); ins.set_ylim(0.12, 0.26); ins.tick_params(labelsize=6); ins.grid(True, color=GRID, lw=0.5)
axs[0, 0].legend(fontsize=7, frameon=False, loc='upper left')
fig.suptitle("Jesper's wave on m3: surface elevation at whole periods vs the exact wave (periodic one-wavelength domain; inset: crest)", fontsize=11, color=INK)
fig.tight_layout(); fig.savefig('figures/jesper_shape.png', dpi=120)
print('done')
