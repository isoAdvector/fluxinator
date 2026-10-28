#!/usr/bin/env python3
#
# Closed-form theory for the first-time-step spurious currents of the standard
# ("basic") gravity term in the tilted box, checked against the OpenFOAM scans.
#
# Summary of the theory (see the paper, "Which cells are hit"):
#
#   1. Starting from rest the pressure equation removes exactly the discrete
#      gradient part of the face gravity force, so the velocity after the first
#      time step is driven only by what is left: the imbalance of that force
#      around closed loops of cells.  On a two-dimensional mesh the smallest
#      loops are the four cells around an interior vertex, and the imbalance of
#      every other loop is a sum of these, so the forcing is one scalar per
#      vertex.  Each term |d_f| F_g,f is the jump in p_rgh that would balance
#      the force, so the loop imbalance is measured in Pa: it is the amount by
#      which those jumps fail to close around the loop.
#
#   2. On a uniform square mesh of side h tilted by theta, with a flat interface
#      at height z_w, that imbalance is a universal function of a single offset
#
#          xi_v = z_w - z_v            (interface height above the vertex)
#
#          Gamma(xi) = g drho { b [A(xi+a) - A(xi-a)] - a [A(xi+b) - A(xi-b)] }
#
#      with a = (h/2)(cos t - sin t), b = (h/2)(cos t + sin t) the two corner
#      offsets of a cell and A the exact volume fraction of a square cut by a
#      horizontal plane.  Gamma vanishes for |xi| >= 2b = h(cos t + sin t) and
#      peaks where the interface passes exactly through the vertex,
#
#          Gamma_max = (1/2) g drho h (cos t - sin t) tan t.
#
#   3. The vertex offsets form a lattice.  Writing, for vertex column I,
#
#          m(I) = (I - N/2) tan t + N/2 + 1 + (z_w - z_0)/(h cos t)
#          K    = round(m),   delta = m - K,   xi = h cos t (delta - 1)
#
#      the row K is the vertex row of that column that carries the interface and
#      delta is a fractional offset.  Which columns are hit is therefore set by how
#      close I tan t + const comes to an integer - an equidistribution problem,
#      which is why increasing the tilt moves the hot cells around instead of
#      simply making the currents stronger.
#
#   4. The velocity an imbalance produces scales as dt/(rho h), so a loop lying
#      entirely in the light fluid responds rho1/rho2 times more strongly than
#      one straddling the interface.  That happens exactly for
#
#          delta in (-tan t, 0],
#
#      a window of width tan t, at whose edge Gamma is still half its peak.  The
#      worst vertex is the one whose delta falls closest to this window, and the
#      cells hit are its four neighbours - in practice the light-fluid one.
#
# Usage:
#   ./vertexOffsetTheory.py                            # closed-form identities
#   ./vertexOffsetTheory.py --validate SCAN [SCAN ...] # against OpenFOAM scans

import argparse

import numpy as np

import tiltedBoxCase as tbc

G = 9.81
RHO1 = 1000.0
RHO2 = 1.0
DT = 0.01

# Response constant: max|U| produced by a unit loop imbalance at a vertex whose
# four cells all have density rho is C0*dt/(h*rho).  A property of the 5-point
# Laplacian on this mesh, not of the case.
C0 = 0.2472


# ---------------------------------------------------------------------------
#  Exact volume fraction of a square cell cut by a horizontal plane
# ---------------------------------------------------------------------------

def alphaSquare(u, theta, h):
    """Water fraction of a square cell of side h whose centre lies a height u
    below the water plane (u = z_w - z_C), for a mesh tilted by theta [rad]."""
    p = h*np.sin(theta)                       # <= q for theta <= 45 deg
    q = h*np.cos(theta)
    b = 0.5*(p + q)
    up = np.asarray(u, float) + b             # 0 .. 2b across the cell
    a = np.zeros_like(up)
    if p < 1e-14:                             # mesh aligned with gravity
        return np.clip(up/q, 0.0, 1.0)
    m1 = (up > 0) & (up <= p); a[m1] = up[m1]**2/(2*p*q)
    m2 = (up > p) & (up <= q); a[m2] = (2*up[m2] - p)/(2*q)
    m3 = (up > q) & (up < 2*b); a[m3] = 1.0 - (2*b - up[m3])**2/(2*p*q)
    a[up >= 2*b] = 1.0
    return a


def cornerOffsets(theta, h):
    """The two corner heights a <= b of a cell of side h tilted by theta."""
    return 0.5*h*(np.cos(theta) - np.sin(theta)), 0.5*h*(np.cos(theta) + np.sin(theta))


# ---------------------------------------------------------------------------
#  Universal vertex loop imbalance
# ---------------------------------------------------------------------------

def Gamma(xi, thetaDeg, h, g=G, drho=RHO1-RHO2):
    """Loop imbalance [Pa] of the basic gravity term around the four-cell loop at a
    vertex whose offset is xi = z_w - z_vertex."""
    t = np.radians(thetaDeg)
    a, b = cornerOffsets(t, h)
    A = lambda u: alphaSquare(u, t, h)
    xi = np.atleast_1d(np.asarray(xi, float))
    return g*drho*(b*(A(xi + a) - A(xi - a)) - a*(A(xi + b) - A(xi - b)))


def GammaPeak(thetaDeg, h, g=G, drho=RHO1-RHO2):
    """Gamma(0): the interface passing exactly through the vertex."""
    t = np.radians(thetaDeg)
    return 0.5*g*drho*h*(np.cos(t) - np.sin(t))*np.tan(t)


# ---------------------------------------------------------------------------
#  Ranking the vertices
# ---------------------------------------------------------------------------

def loopAlphas(xi, theta, h):
    """Volume fractions of the four cells sharing a vertex of offset xi."""
    c, s = np.cos(theta), np.sin(theta)
    off = np.array([(-0.5, -0.5), (0.5, -0.5), (-0.5, 0.5), (0.5, 0.5)])*h
    dz = np.array([o1*(-s) + o3*c for o1, o3 in off])      # heights above the vertex
    return alphaSquare((np.atleast_1d(xi)[None, :] - dz[:, None]).ravel(),
                       theta, h).reshape(4, -1)


def amplitude(xi, thetaDeg, h, dt=DT, rho1=RHO1, rho2=RHO2, g=G):
    """Peak speed expected from a vertex of offset xi: the imbalance divided by
    the density of the loop, whose response sets how much velocity it buys."""
    t = np.radians(thetaDeg)
    a4 = np.sort(loopAlphas(xi, t, h), axis=0)
    rhoEff = rho2 + (rho1 - rho2)*a4[2]                    # third-lightest cell
    return C0*dt*np.abs(Gamma(xi, thetaDeg, h, g, rho1-rho2))/(h*rhoEff)


def lightFluidWindow(xi, thetaDeg, h):
    """Vertices whose four-cell loop lies entirely in the light fluid while the
    imbalance is still non-zero: -h(cos t + sin t) < xi <= -h cos t.  In terms
    of the fractional offset delta of columns() this is -tan t < delta <= 0."""
    t = np.radians(thetaDeg)
    return (xi > -h*(np.cos(t) + np.sin(t))) & (xi <= -h*np.cos(t))


def rankVertices(V, zw, thetaDeg, dt=DT, rho1=RHO1, rho2=RHO2, g=G):
    """Rank the interior vertices of a mesh by the spurious velocity expected
    there.  V is the (N+1, N+1, 2) grid of mesh vertices, so the offset is read
    straight from the mesh.  Returns the interior arrays (K, I indexing) of
    offset, loop imbalance and amplitude."""
    h = np.linalg.norm(V[0, 1] - V[0, 0])
    xi = zw - V[1:-1, 1:-1, 1]
    G_ = Gamma(xi.ravel(), thetaDeg, h, g, rho1-rho2).reshape(xi.shape)
    A = amplitude(xi.ravel(), thetaDeg, h, dt, rho1, rho2, g).reshape(xi.shape)
    return dict(xi=xi, Gamma=G_, A=A, h=h,
                window=lightFluidWindow(xi, thetaDeg, h))


def loopCells(I, K, N=25):
    """The four cells sharing vertex (I, K), as (i, k) cell indices."""
    return [(i, k) for i, k in ((I-1, K-1), (I, K-1), (I-1, K), (I, K))
            if 0 <= i < N and 0 <= k < N]


def windowCells(V, zw, thetaDeg, N=25, **kw):
    """Every cell touching a vertex in the light-fluid window: the cells the
    theory says can carry a strongly amplified spurious velocity."""
    r = rankVertices(V, zw, thetaDeg, **kw)
    out = set()
    for k, i in np.argwhere(r['window']):
        out.update(loopCells(i+1, k+1, N))
    return out, r


def predictFromMesh(V, zw, thetaDeg, N=25, **kw):
    """The vertex expected to carry the largest spurious velocity, its four
    cells, and the predicted max|U|."""
    r = rankVertices(V, zw, thetaDeg, **kw)
    k, i = np.unravel_index(int(np.argmax(r['A'])), r['A'].shape)
    return dict(vertex=(i+1, k+1), cells=loopCells(i+1, k+1, N),
                cell=(i, k), xi=float(r['xi'][k, i]),
                Gamma=float(r['Gamma'][k, i]), maxU=float(r['A'][k, i]), rank=r)


# ---------------------------------------------------------------------------
#  The same thing without a mesh: the offset lattice
# ---------------------------------------------------------------------------

def columns(thetaDeg, zw=0.5, N=25, L=1.0, z0=None, **kw):
    """Per vertex-column I of an N x N box rotated about its centre at height
    z0: the interface-carrying row K, the fractional offset delta, the offset, the
    loop imbalance and the expected peak speed - from geometry alone."""
    h = L/N
    t = np.radians(thetaDeg)
    c = np.cos(t)
    if z0 is None:
        z0 = L/2
    I = np.arange(1, N)
    m = (I - N/2)*np.tan(t) + N/2 + 1 + (zw - z0)/(h*c)
    K = np.rint(m).astype(int)
    delta = m - K
    xi = h*c*(delta - 1.0)
    ok = (K >= 1) & (K <= N-1)
    return dict(I=I, K=K, delta=delta, xi=xi, h=h,
                Gamma=Gamma(xi, thetaDeg, h),
                A=np.where(ok, amplitude(xi, thetaDeg, h, **kw), 0.0), ok=ok,
                airWindow=(delta > -np.tan(t)) & (delta <= 0.0))


def predict(thetaDeg, zw=0.5, N=25, L=1.0, z0=None, **kw):
    """predictFromMesh without needing a mesh."""
    d = columns(thetaDeg, zw, N, L, z0, **kw)
    j = int(np.argmax(d['A']))
    return dict(vertex=(int(d['I'][j]), int(d['K'][j])),
                cell=(int(d['I'][j]) - 1, int(d['K'][j]) - 1),
                delta=float(d['delta'][j]), xi=float(d['xi'][j]),
                Gamma=float(d['Gamma'][j]), maxU=float(d['A'][j]), columns=d)


# ---------------------------------------------------------------------------

def identities(N=25, L=1.0):
    h = L/N
    print('Closed-form identities of the universal loop imbalance')
    for th in (5.0, 12.0):
        t = np.radians(th)
        xi = np.linspace(-2.5*h, 2.5*h, 200001)
        g = Gamma(xi, th, h)
        i = int(np.argmax(g))
        sup = np.abs(xi[np.abs(g) > 1e-9]).max()/h
        edge = Gamma(-h*np.cos(t), th, h)[0]
        print(f'  theta = {th:5.1f} deg')
        print(f'    peak           {g[i]:9.4f} Pa at xi/h = {xi[i]/h:+.5f}; '
              f'closed form {GammaPeak(th, h):9.4f} Pa')
        print(f'    support        |xi|/h < {sup:.5f}   (cos+sin = {np.cos(t)+np.sin(t):.5f})')
        print(f'    light-fluid window width tan(theta) = {np.tan(t):.5f}; '
              f'Gamma at its edge {edge:9.4f} Pa = {edge/GammaPeak(th, h):+.5f} x peak')


def validate(scanDirs, N=25):
    """Compare the prediction with the OpenFOAM runs of one or more scans."""
    inWindow, inLoop, qual, ratio, nWin, rows = [], [], [], [], [], []
    for scan in scanDirs:
        for case in tbc.cases(scan):
            r = tbc.firstStep(case, N)
            mU = r['magU']
            kk, ii = np.unravel_index(int(np.argmax(mU)), mU.shape)
            cells, rk = windowCells(r['V'], r['zw'], r['theta'], N)
            p = predictFromMesh(r['V'], r['zw'], r['theta'], N)
            inWindow.append((ii, kk) in cells)
            inLoop.append((ii, kk) in p['cells'])
            qual.append(max(mU[k, i] for i, k in p['cells'])/mU.max())
            ratio.append(mU.max()/p['maxU'])
            nWin.append(int(rk['window'].sum()))
            rows.append((r['theta'], r['zw'], mU.max(), p['maxU'], (ii, kk),
                         p['vertex'], nWin[-1], inWindow[-1], qual[-1]))
    inWindow = np.array(inWindow); inLoop = np.array(inLoop)
    qual = np.array(qual); ratio = np.array(ratio); nWin = np.array(nWin)
    print(f'Prediction against {len(qual)} interFoam runs ({", ".join(scanDirs)})')
    print(f'  vertices in the light-fluid window: {nWin.min()} to {nWin.max()} '
          f'per run (N tan(theta) = {25*np.tan(np.radians(2)):.1f} to '
          f'{25*np.tan(np.radians(20)):.1f} over 2-20 deg)')
    print(f'  the peak cell touches one of them in {100*inWindow.mean():.1f}% of runs')
    print(f'  the peak cell is one of the four cells of the highest-ranked vertex '
          f'in {100*inLoop.mean():.1f}% of runs;')
    print(f'    the best of those four carries a mean {qual.mean():.3f} of the peak '
          f'velocity (>= 0.90 in {100*np.mean(qual >= 0.90):.1f}%, worst {qual.min():.3f})')
    print(f'  true max|U| / predicted: median {np.median(ratio):.3f}, '
          f'5-95% [{np.percentile(ratio, 5):.3f}, {np.percentile(ratio, 95):.3f}]')
    print()
    print(f'  {"tilt":>6} {"z_w":>7} {"max|U|":>10} {"predicted":>10} '
          f'{"peak cell":>10} {"top vertex":>11} {"#window":>8} {"in window":>10}')
    step = max(1, len(rows)//24)
    for th, zw, mu, pm, pc, vx, nw, iw, q in rows[::step]:
        print(f'  {th:6.2f} {zw:7.3f} {mu:10.5f} {pm:10.5f} '
              f'{str(pc):>10} {str(vx):>11} {nw:8d} {str(iw):>10}')


def latticeCheck(scanDirs, N=25):
    """The closed-form offset lattice against the offsets read from the meshes."""
    err = []
    for scan in scanDirs:
        for case in tbc.cases(scan):
            r = tbc.firstStep(case, N)
            d = columns(r['theta'], r['zw'], N)
            xiMesh = r['zw'] - r['V'][1:-1, 1:-1, 1]
            for j in range(len(d['I'])):
                if d['ok'][j]:
                    err.append(abs(d['xi'][j] - xiMesh[d['K'][j]-1, d['I'][j]-1]))
    print(f'Closed-form offset lattice against the mesh vertices of the same runs: '
          f'max error {max(err):.1e} m')


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--validate', nargs='+', metavar='SCAN',
                    default=['tiltScan_1deltaT', 'waterLevelScan_1deltaT'],
                    help='parmScanner scan directories of one-step interFoam runs')
    ap.add_argument('-N', type=int, default=25, help='cells per side (default 25)')
    a = ap.parse_args()
    identities(a.N)
    print()
    latticeCheck(a.validate, a.N)
    print()
    validate(a.validate, a.N)
