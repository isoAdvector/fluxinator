"""Model of explicit momentum convection at the faces near a wave crest.
Two-phase linear wave on the 190x40 mesh of the paper (dx = L/190, dz = 0.02, water depth 0.4,
air 0.4 under a lid). The flow is given by a stream function in each phase; the surface is the
zero streamline of the water flow in the wave frame, so it is exactly material and all volume
fluxes follow exactly from stream-function differences. The air flows against the water at the
crest (linear air potential, plus a small correction that puts its zero streamline on the surface).
For one step t^n -> t^n + dt the script computes, per face, the exact mass F_f and momentum M_f
that cross it (split into water and air), and per cell the exact mass and momentum at t^n and
t^{n+1}. The cell velocity U^n is the mass-weighted one. Explicit convection replaces M_f by
F_f U_f with U_f from U^n; its error in cell i is
    e_i = - sum_f s_if (F_f U_f - M_f) / m_i^{n+1}   (s_if = +1 for outflow).
The ideal face value is ubar_f = M_f / F_f, the mass-weighted velocity of the fluid that crosses.
Usage: python3 crestFaces.py [dt]   (default 0.005)"""
import sys
import numpy as np

RW, RA, G = 1000.0, 1.0, 9.81
L, T, D, HA, AMP = 3.7873860376294, 2.0, 0.4, 0.4, 0.0625
NX, DZ = 190, 0.02
DX = L/NX
K, C = 2*np.pi/L, L/T
OM = K*C
XC = 0.5*L                       # crest position at t^n
DT = float(sys.argv[1]) if len(sys.argv) > 1 else 0.005
ZB = np.arange(-0.16, 0.16 + 1e-9, DZ)          # row edges of the band analysed
NZ = len(ZB) - 1
XE = np.arange(NX + 1)*DX
GLX, GLW = np.polynomial.legendre.leggauss(16)
NT = 64
TX, TW = np.polynomial.legendre.leggauss(NT)
TN = 0.5*DT*(TX + 1); TWN = 0.5*DT*TW            # time nodes and weights on [0, dt]


def Sw(z): return np.sinh(K*(z + D))/np.sinh(K*D)
def Cw(z): return np.cosh(K*(z + D))/np.sinh(K*D)
def Sa(z): return np.sinh(K*(HA - z))/np.sinh(K*HA)
def Ca(z): return np.cosh(K*(HA - z))/np.sinh(K*HA)


def eta(th):
    """Surface: eta = AMP Sw(eta) cos th (water zero streamline in the wave frame)."""
    e = AMP*np.cos(th)
    for _ in range(30):
        f = e - AMP*Sw(e)*np.cos(th)
        e = e - f/(1 - AMP*K*Cw(e)*np.cos(th))
    return e


def deta(th):
    e = eta(th)
    return -AMP*Sw(e)*np.sin(th)/(1 - AMP*K*Cw(e)*np.cos(th))


def theta(x, t): return K*(x - XC - C*t)


def psi(x, z, t, water):
    th = theta(x, t)
    if water:
        return C*AMP*Sw(z)*np.cos(th)
    e = eta(th)
    dl = C*e - C*AMP*Sa(e)*np.cos(th)
    return C*AMP*Sa(z)*np.cos(th) + dl*(HA - z)/(HA - e)


def vel(x, z, t, water):
    th = theta(x, t)
    if water:
        return C*AMP*K*Cw(z)*np.cos(th), C*AMP*K*Sw(z)*np.sin(th)
    e, de = eta(th), deta(th)
    dl = C*e - C*AMP*Sa(e)*np.cos(th)
    ddl = C*de - C*AMP*(-K*Ca(e)*de*np.cos(th) - Sa(e)*np.sin(th))
    u = -C*AMP*K*Ca(z)*np.cos(th) - dl/(HA - e)
    dpsi = -C*AMP*Sa(z)*np.sin(th) + ddl*(HA - z)/(HA - e) + dl*(HA - z)*de/(HA - e)**2
    return u, -K*dpsi


def gl(a, b):
    """Gauss nodes/weights on [a, b] (arrays broadcast, last axis = nodes)."""
    a, b = np.asarray(a)[..., None], np.asarray(b)[..., None]
    return 0.5*(b - a)*GLX + 0.5*(b + a), 0.5*(b - a)*GLW


def crossings(x0, x1, t, zl):
    """x where eta(x, t) = zl inside [x0, x1] (up to two), returned as sorted breakpoints
    padded with x0 (arrays of shape (..., 2))."""
    if zl >= eta(0.0) or zl <= eta(np.pi):
        return np.stack([x0 + 0*t, x0 + 0*t], -1)
    lo, hi = 0.0, np.pi                       # eta decreasing on [0, pi]
    for _ in range(60):
        mid = 0.5*(lo + hi)
        lo, hi = (mid, hi) if eta(mid) > zl else (lo, mid)
    ts = 0.5*(lo + hi)
    out = []
    for sgn in (1, -1):
        th0 = theta(x0, t)
        thc = sgn*ts + 2*np.pi*np.ceil((th0 - sgn*ts)/(2*np.pi))         # first crossing >= th0
        xc = x0 + (thc - th0)/K
        out.append(np.where(xc < x1, xc, x0))
    return np.sort(np.stack(out, -1), -1)


def seg_x(x0, x1, t, levels):
    """Breakpoints of [x0, x1] at the surface crossings of the given levels: (..., nseg+1)."""
    bp = [np.broadcast_to(x0, np.shape(t)), np.broadcast_to(x1, np.shape(t))]
    for zl in levels:
        cr = crossings(x0, x1, t, zl)
        bp += [cr[..., 0], cr[..., 1]]
    bp = np.sort(np.stack(bp, -1), -1)
    return bp


# ---------------------------------------------------------------- vertical faces (normal +x)
def vface(xf, z0, z1):
    """Mass and momentum (x, z) crossing the face over the step, per phase: dict of arrays."""
    t = TN
    e = eta(theta(xf, t))
    ew = np.clip(e, z0, z1)
    res = {}
    for ph, (a, b) in (("w", (z0 + 0*t, ew)), ("a", (ew, z1 + 0*t))):
        water = ph == "w"
        vol = psi(xf, b, t, water) - psi(xf, a, t, water)
        zn, wn = gl(a, b)
        u, w = vel(xf, zn, t[:, None], water)
        rho = RW if water else RA
        res["F" + ph] = rho*np.sum(TWN*vol)
        res["Mx" + ph] = rho*np.sum(TWN*np.sum(wn*u*u, -1))
        res["Mz" + ph] = rho*np.sum(TWN*np.sum(wn*u*w, -1))
    return res


# ---------------------------------------------------------------- horizontal faces (normal +z)
def hface(zf, x0, x1):
    t = TN
    bp = seg_x(x0, x1, t, [zf])
    res = {k: 0.0 for k in ("Fw", "Fa", "Mxw", "Mxa", "Mzw", "Mza")}
    for s in range(bp.shape[-1] - 1):
        a, b = bp[:, s], bp[:, s + 1]
        mid = 0.5*(a + b)
        water = eta(theta(mid, t)) > zf
        for ph, wmask in (("w", water), ("a", ~water)):
            isw = ph == "w"
            rho = RW if isw else RA
            vol = -(psi(b, zf, t, isw) - psi(a, zf, t, isw))
            xn, wn = gl(a, b)
            u, w = vel(xn, zf, t[:, None], isw)
            m = wmask*(b > a)
            res["F" + ph] += rho*np.sum(TWN*m*vol)
            res["Mx" + ph] += rho*np.sum(TWN*m*np.sum(wn*w*u, -1))
            res["Mz" + ph] += rho*np.sum(TWN*m*np.sum(wn*w*w, -1))
    return res


# ---------------------------------------------------------------- cells
def cell(x0, x1, z0, z1, t):
    """Water volume, mass and momentum of the cell at time t."""
    tt = np.array([t])
    bp = seg_x(x0, x1, tt, [z0, z1])[0]
    vw, mx, mz = 0.0, 0.0, 0.0
    for s in range(len(bp) - 1):
        a, b = bp[s], bp[s + 1]
        if b <= a:
            continue
        xn, wn = gl(a, b)
        e = np.clip(eta(theta(xn, t)), z0, z1)
        vw += np.sum(wn*(e - z0))
        for isw, (za, zb) in ((True, (z0 + 0*e, e)), (False, (e, z1 + 0*e))):
            rho = RW if isw else RA
            mx += rho*np.sum(wn*(psi(xn, zb, t, isw) - psi(xn, za, t, isw)))
            zn, zw = gl(za, zb)
            _, w = vel(xn[:, None], zn, t, isw)
            mz += rho*np.sum(wn*np.sum(zw*w, -1))
    va = (x1 - x0)*(z1 - z0)
    return vw/va, RW*vw + RA*(va - vw), np.array([mx, mz])


def run():
    al0 = np.zeros((NZ, NX)); al1 = np.zeros((NZ, NX))
    m0 = np.zeros((NZ, NX)); m1 = np.zeros((NZ, NX))
    P0 = np.zeros((NZ, NX, 2)); P1 = np.zeros((NZ, NX, 2))
    for j in range(NZ):
        for i in range(NX):
            al0[j, i], m0[j, i], P0[j, i] = cell(XE[i], XE[i + 1], ZB[j], ZB[j + 1], 0.0)
            al1[j, i], m1[j, i], P1[j, i] = cell(XE[i], XE[i + 1], ZB[j], ZB[j + 1], DT)
    # vertical face i at x = XE[i] between cells i-1 and i (periodic: face 0 = face NX)
    V = {k: np.zeros((NZ, NX)) for k in ("Fw", "Fa", "Mxw", "Mxa", "Mzw", "Mza")}
    for j in range(NZ):
        for i in range(NX):
            r = vface(XE[i], ZB[j], ZB[j + 1])
            for k in V:
                V[k][j, i] = r[k]
    # horizontal face j at z = ZB[j] between rows j-1 and j (only interior faces 1..NZ-1)
    Hf = {k: np.zeros((NZ + 1, NX)) for k in V}
    for j in range(1, NZ):
        for i in range(NX):
            r = hface(ZB[j], XE[i], XE[i + 1])
            for k in Hf:
                Hf[k][j, i] = r[k]
    return dict(al0=al0, al1=al1, m0=m0, m1=m1, P0=P0, P1=P1, V=V, H=Hf)


if __name__ == "__main__":
    import os, pickle
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", "crestFaces_dt%g.pkl" % DT)
    r = run()
    pickle.dump(r, open(out, "wb"))
    # mass balance check on interior rows
    F = r["V"]["Fw"] + r["V"]["Fa"]; Fh = r["H"]["Fw"] + r["H"]["Fa"]
    div = (np.roll(F, -1, 1) - F) + (Fh[1:] - Fh[:-1])
    res = (r["m1"] - r["m0"] + div)[1:-1]
    print("mass residual max %.3e kg/m (air mass per cell %.3e)" % (np.abs(res).max(), RA*DX*DZ))
    print("written", out)
