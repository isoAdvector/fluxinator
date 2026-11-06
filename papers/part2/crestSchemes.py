# Copyright (C) 2026 Johan Roenby
# SPDX-License-Identifier: GPL-3.0-or-later
"""Face-value schemes in the crest model of crestFaces.py.
Reads data/crestFaces_dt<dt>.pkl and evaluates, per scheme, the velocity error that explicit
convection puts into each cell in one step,
    e_i = - sum_f s_if (F_f U_f - M_f) / m_i^{n+1},
with U_f from the mass-weighted cell velocities U^n and the exact mass flux F_f.
Usage: python3 crestSchemes.py [dt]"""
import os
import sys
import pickle
import numpy as np
from crestFaces import DX, DZ, ZB, XE, K, XC, NX, RW, RA

DT = float(sys.argv[1]) if len(sys.argv) > 1 else 0.005
HERE = os.path.dirname(os.path.abspath(__file__))
R = pickle.load(open(os.path.join(HERE, "data", "crestFaces_dt%g.pkl" % DT), "rb"))
EPS = 1e-6
al0, m0, m1 = R["al0"], R["m0"], R["m1"]
U0 = R["P0"]/m0[..., None]
U1 = R["P1"]/m1[..., None]
NZ = al0.shape[0]

# Faces as lists: owner P, neighbour N (cell indices), normal, F, M (vector), Fw, Mw, Fa, Ma
fP, fN, fd, F, M, Fw, Mw, Fa, Ma = [], [], [], [], [], [], [], [], []
V, H = R["V"], R["H"]
for j in range(NZ):
    for i in range(NX):                       # vertical face at XE[i]: P = (j, i-1), N = (j, i)
        fP.append((j, (i - 1) % NX)); fN.append((j, i)); fd.append(0)
        Fw.append(V["Fw"][j, i]); Fa.append(V["Fa"][j, i])
        Mw.append((V["Mxw"][j, i], V["Mzw"][j, i])); Ma.append((V["Mxa"][j, i], V["Mza"][j, i]))
for j in range(1, NZ):
    for i in range(NX):                       # horizontal face at ZB[j]: P = (j-1, i), N = (j, i)
        fP.append((j - 1, i)); fN.append((j, i)); fd.append(1)
        Fw.append(H["Fw"][j, i]); Fa.append(H["Fa"][j, i])
        Mw.append((H["Mxw"][j, i], H["Mzw"][j, i])); Ma.append((H["Mxa"][j, i], H["Mza"][j, i]))
fP, fN, fd = np.array(fP), np.array(fN), np.array(fd)
Fw, Fa, Mw, Ma = np.array(Fw), np.array(Fa), np.array(Mw), np.array(Ma)
F, M = Fw + Fa, Mw + Ma
nf = len(F)
UP, UN = U0[fP[:, 0], fP[:, 1]], U0[fN[:, 0], fN[:, 1]]
aP, aN = al0[fP[:, 0], fP[:, 1]], al0[fN[:, 0], fN[:, 1]]
pos = F > 0

d = np.where(fd[:, None] == 0, np.array([DX, 0.0]), np.array([0.0, DZ]))
NB = np.full((NZ, NX, 4, 2), -1)               # cell neighbours (E, W, T, B), -1 outside the band
for j in range(NZ):
    for i in range(NX):
        NB[j, i] = [(j, (i + 1) % NX), (j, (i - 1) % NX), (j + 1, i) if j < NZ - 1 else (-1, -1),
                    (j - 1, i) if j > 0 else (-1, -1)]
mw0 = RW*al0*DX*DZ; ma0 = m0 - mw0             # phase masses at t^n


def gradient(U):
    """Cell gradient (central differences = Gauss linear on the uniform mesh), g[j, i, a, b] = d_a U_b."""
    g = np.zeros((NZ, NX, 2, 2))
    g[:, :, 0, :] = (np.roll(U, -1, 1) - np.roll(U, 1, 1))/(2*DX)
    g[1:-1, :, 1, :] = (U[2:] - U[:-2])/(2*DZ)
    return g


def limitedLinearV(U, k=1.0):
    UP, UN = U[fP[:, 0], fP[:, 1]], U[fN[:, 0], fN[:, 1]]
    grad = gradient(U)
    gP, gN = grad[fP[:, 0], fP[:, 1]], grad[fN[:, 0], fN[:, 1]]
    gfV = UN - UP
    gf = np.sum(gfV*gfV, 1)
    gc = np.where(pos[:, None], np.einsum("fa,fab->fb", d, gP), np.einsum("fa,fab->fb", d, gN))
    gcf = np.sum(gfV*gc, 1)
    with np.errstate(divide="ignore", invalid="ignore"):
        r = np.where(np.abs(gcf) >= 1000*np.abs(gf), 2*1000*np.sign(gcf)*np.sign(gf) - 1, 2*gcf/gf - 1)
    lim = np.clip(2/k*r, 0, 1)
    w = lim*0.5 + (1 - lim)*pos
    return w[:, None]*UP + (1 - w[:, None])*UN


def phase_velocities(U, limit=False):
    """Estimate of the water and air velocity in each cell. The minority phase (by mass) takes the
    velocity of the neighbour richest in that phase; the majority phase follows from the cell
    momentum, m U = m_w u_w + m_a u_a."""
    uw, ua = U.copy(), U.copy()
    for j in range(NZ):
        for i in range(NX):
            if not (EPS < al0[j, i] < 1 - EPS):
                continue
            nb = [tuple(n) for n in NB[j, i] if n[0] >= 0]
            if mw0[j, i] >= ma0[j, i]:
                n = min(nb, key=lambda c: al0[c])
                if al0[n] < al0[j, i]:
                    ua[j, i] = U[n]
                    uw[j, i] = (m0[j, i]*U[j, i] - ma0[j, i]*ua[j, i])/mw0[j, i]
            else:
                n = max(nb, key=lambda c: al0[c])
                if al0[n] > al0[j, i]:
                    uw[j, i] = U[n]
                    ua[j, i] = (m0[j, i]*U[j, i] - mw0[j, i]*uw[j, i])/ma0[j, i]
    if limit:                                   # keep both within the local range of U
        lo, hi = U.copy(), U.copy()
        for n in range(4):
            jj, ii = NB[..., n, 0], NB[..., n, 1]
            ok = (jj >= 0)[..., None]
            Un = U[np.clip(jj, 0, None), ii]
            lo = np.where(ok, np.minimum(lo, Un), lo); hi = np.where(ok, np.maximum(hi, Un), hi)
        uw, ua = np.clip(uw, lo, hi), np.clip(ua, lo, hi)
    return uw, ua


def fifo_labels(Fk, uk, mk):
    """Velocity carried by the outflow of one phase from each cell: the first m_k^n carry the
    cell's u_k, any excess (fluid that entered and left in the same step) carries the mass-weighted
    velocity of the inflow of that phase, from the upwind neighbours' u_k."""
    out = np.zeros((NZ, NX)); fin = np.zeros((NZ, NX)); pin = np.zeros((NZ, NX, 2))
    p = Fk > 0
    np.add.at(out, (fP[p, 0], fP[p, 1]), Fk[p]); np.add.at(out, (fN[~p, 0], fN[~p, 1]), -Fk[~p])
    np.add.at(fin, (fN[p, 0], fN[p, 1]), Fk[p]); np.add.at(fin, (fP[~p, 0], fP[~p, 1]), -Fk[~p])
    np.add.at(pin, (fN[p, 0], fN[p, 1]), Fk[p, None]*uk[fP[p, 0], fP[p, 1]])
    np.add.at(pin, (fP[~p, 0], fP[~p, 1]), -Fk[~p, None]*uk[fN[~p, 0], fN[~p, 1]])
    exc = np.clip(out - mk, 0, None)
    with np.errstate(divide="ignore", invalid="ignore"):
        uin = np.where(fin[..., None] > 0, pin/fin[..., None], uk)
        lab = np.where(out[..., None] > 0, ((out - exc)[..., None]*uk + exc[..., None]*uin)/out[..., None], uk)
    return lab


def face_flux(U, scheme):
    """Numerical momentum through each face (along +normal) over the step."""
    UP, UN = U[fP[:, 0], fP[:, 1]], U[fN[:, 0], fN[:, 1]]
    if scheme == "upwind":
        return F[:, None]*np.where(pos[:, None], UP, UN)
    if scheme == "linear":
        return F[:, None]*0.5*(UP + UN)
    if scheme == "limitedLinearV":
        return F[:, None]*limitedLinearV(U)
    if scheme == "interfaceUpwind":
        return F[:, None]*np.where(pure[:, None], limitedLinearV(U), np.where(pos[:, None], UP, UN))
    if scheme == "phaseUpwind":
        uw, ua = phase_velocities(U)
        uwP, uwN = uw[fP[:, 0], fP[:, 1]], uw[fN[:, 0], fN[:, 1]]
        uaP, uaN = ua[fP[:, 0], fP[:, 1]], ua[fN[:, 0], fN[:, 1]]
        ph = Fw[:, None]*np.where(Fw[:, None] > 0, uwP, uwN) + Fa[:, None]*np.where(Fa[:, None] > 0, uaP, uaN)
        return np.where(pure[:, None], F[:, None]*limitedLinearV(U), ph)
    if scheme in ("upwindFIFO", "phaseUpwindFIFO", "phaseUpwindFIFOlim"):
        if scheme == "upwindFIFO":
            phases = [(F, U, m0)]
        else:
            uw, ua = phase_velocities(U, limit=scheme.endswith("lim"))
            phases = [(Fw, uw, mw0), (Fa, ua, ma0)]
        flux = np.zeros((nf, 2))
        for Fk, uk, mk in phases:
            lab = fifo_labels(Fk, uk, mk)
            flux += Fk[:, None]*np.where(Fk[:, None] > 0, lab[fP[:, 0], fP[:, 1]], lab[fN[:, 0], fN[:, 1]])
        if scheme == "upwindFIFO":
            return flux
        return np.where(pure[:, None], F[:, None]*limitedLinearV(U), flux)
    if scheme == "ideal":
        return M.copy()
    raise ValueError(scheme)


pure = ((aP <= EPS) & (aN <= EPS)) | ((aP >= 1 - EPS) & (aN >= 1 - EPS))
with np.errstate(divide="ignore", invalid="ignore"):
    ideal = np.where(np.abs(F)[:, None] > 0, M/F[:, None], 0.0)
NAMES = ["upwind", "linear", "limitedLinearV", "interfaceUpwind", "phaseUpwind", "upwindFIFO",
         "phaseUpwindFIFO", "phaseUpwindFIFOlim", "ideal"]


def step(U, scheme):
    acc = np.zeros((NZ, NX, 2))
    flux = face_flux(U, scheme)
    np.add.at(acc, (fP[:, 0], fP[:, 1]), flux)
    np.add.at(acc, (fN[:, 0], fN[:, 1]), -flux)
    return (m0[..., None]*U - acc)/m1[..., None]


def cell_error(scheme):
    exact = step(U0, "ideal")
    return step(U0, scheme) - exact


def gain(scheme, cells, h=1e-3):
    """Sum over k of |d U1_i / d U0_k| (x component) for the cells i given, by perturbation."""
    base = step(U0, scheme)[..., 0]
    near = set()
    for (j, i) in cells:
        near |= {(j, i)} | {tuple(n) for n in NB[j, i] if n[0] >= 0}
        for n in NB[j, i]:
            if n[0] >= 0:
                near |= {tuple(m) for m in NB[tuple(n)] if m[0] >= 0}
    g = {c: 0.0 for c in cells}
    for k in near:
        Up = U0.copy(); Up[k][0] += h
        dU = (step(Up, scheme)[..., 0] - base)/h
        for c in cells:
            g[c] += abs(dU[c])
    return g


ERR = {k: cell_error(k) for k in NAMES}
inner = np.zeros((NZ, NX), bool); inner[1:-1] = True
iface = inner & ((al0 > EPS) & (al0 < 1 - EPS) | (R["al1"] > EPS) & (R["al1"] < 1 - EPS))

if __name__ == "__main__":
    mfr = mw0/m0                                  # water mass fraction at t^n
    wd = iface & (mfr >= 0.5); ad = iface & (mfr < 0.5); pc = inner & ~iface
    cand = [tuple(c) for c in np.argwhere(iface)]
    print("dt %g: %d interface cells (%d water-dominated by mass), water crest speed %.3f m/s"
          % (DT, iface.sum(), wd.sum(), U0[..., 0].max()))
    print("%-16s %10s %10s %10s %10s %10s" % ("one-step error", "max wdom", "max adom", "max pure",
                                             "rms iface", "max gain"))
    for k in NAMES:
        mag = np.linalg.norm(ERR[k], axis=-1)
        g = gain(k, cand) if k != "ideal" else {0: 1.0}
        print("%-16s %10.4f %10.4f %10.4f %10.4f %10.2f" % (k, mag[wd].max(), mag[ad].max(), mag[pc].max(),
              np.sqrt(np.mean(mag[iface]**2)), max(g.values())))
