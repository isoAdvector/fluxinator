#!/usr/bin/env python3
#
# An independent numpy re-implementation of the first interFoam / interIsoFoam
# time step in the tilted box, written while the theory in vertexOffsetTheory.py
# was being developed.  It reproduces OpenFOAM to about 1e-10 relative.
#
# NOTHING IN THE PAPER DEPENDS ON THIS FILE. Every number and figure there comes
# from the OpenFOAM scans.  It is kept because it makes the first step of the
# argument directly checkable: that the velocity after the first time step
# depends on the face gravity force only through its loop imbalance around the
# four-cell loops.  Demonstrating that needs the face force to be replaced by a
# different one with the same imbalance, which stock OpenFOAM cannot do.
#
# The reimplemented step is: rAU = dt/rho (no momentum predictor, no viscosity,
# zero flux so no convection), phig = -ghf*snGrad(rho)*rAUf*|Sf|,
# laplacian(rAUf, p_rgh) = div(phig) with zero flux on all walls, and
# U = rAU*reconstruct((phig - flux)/rAUf).
#
# Usage:
#   ./firstStepTwin.py [--validate tiltScan_1deltaT]

import argparse

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import splu

import tiltedBoxCase as tbc
from vertexOffsetTheory import G, RHO1, RHO2, DT, alphaSquare


def firstStep(N=25, L=1.0, thetaDeg=5.0, zw=0.5, dt=DT, rho1=RHO1, rho2=RHO2,
              g=G, z0=None, thick=1.0, faceForce=None):
    """Velocity after one time step from rest, for the basic gravity term.

    faceForce, if given, replaces the standard term: a pair of arrays
    (i-faces, k-faces) holding F_g on each internal face."""
    h = L/N
    t = np.radians(thetaDeg)
    c, s = np.cos(t), np.sin(t)
    if z0 is None:
        z0 = L/2
    e1, e3 = np.array([c, -s]), np.array([s, c])
    q = (np.arange(N) + 0.5)*h
    X, Z = np.meshgrid(q, q, indexing='xy')
    d = np.stack([X - L/2, Z - z0], -1)
    C = np.array([L/2, z0]) + d[..., 0:1]*e1 + d[..., 1:2]*e3
    alpha = alphaSquare(zw - C[..., 1], t, h)
    rho = alpha*rho1 + (1 - alpha)*rho2
    rAU = dt/rho
    cid = np.arange(N*N).reshape(N, N)
    Af = h*thick
    rows, cols, vals, faces = [], [], [], []
    b = np.zeros(N*N)
    for n, (axis, e) in enumerate(((1, e1), (0, e3))):
        P = cid.take(np.arange(N-1), axis=axis).ravel()
        Q = cid.take(np.arange(1, N), axis=axis).ravel()
        Cf = 0.5*(C.reshape(-1, 2)[P] + C.reshape(-1, 2)[Q])
        rAUf = 0.5*(rAU.ravel()[P] + rAU.ravel()[Q])
        Dc = rAUf*Af/h
        if faceForce is None:
            Fg = -(-g*Cf[:, 1])*(rho.ravel()[Q] - rho.ravel()[P])/h
        else:
            Fg = np.asarray(faceForce[n]).ravel()
        phig = Fg*rAUf*Af
        faces.append((P, Q, phig, Dc, rAUf, e, Fg))
        rows += [P, Q, P, Q]; cols += [Q, P, P, Q]; vals += [Dc, Dc, -Dc, -Dc]
        np.add.at(b, P, phig); np.add.at(b, Q, -phig)
    A = coo_matrix((np.concatenate(vals),
                    (np.concatenate(rows), np.concatenate(cols))),
                   shape=(N*N, N*N)).tolil()
    A[0, :] = 0; A[0, 0] = 1.0; b[0] = 0.0     # pin the pressure reference
    p = splu(A.tocsc()).solve(b)
    num = np.zeros((N*N, 2))
    for (P, Q, phig, Dc, rAUf, e, Fg) in faces:
        ssf = (phig - Dc*(p[Q] - p[P]))/rAUf
        np.add.at(num, P, np.outer(ssf, e))
        np.add.at(num, Q, np.outer(ssf, e))    # surfaceSum: no sign flip
    U = (rAU.ravel()[:, None]*num/(2*Af)).reshape(N, N, 2)
    Fi = faces[0][6].reshape(N, N-1)
    Fk = faces[1][6].reshape(N-1, N)
    imb = h*(Fi[:-1, :] + Fk[:, 1:] - Fi[1:, :] - Fk[:, :-1])
    return dict(alpha=alpha, rho=rho, C=C, p=p.reshape(N, N), U=U, h=h,
                magU=np.linalg.norm(U, axis=2), faceForce=[faces[0][6], faces[1][6]],
                imbalance=imb)


def staircaseFaceForce(imb, h):
    """A face force with the given vertex loop imbalance, put entirely on the
    i-faces.  It differs from the standard term by a discrete gradient."""
    nk, ni = imb.shape
    N = ni + 1
    Fi = np.zeros((N, N-1))
    for k in range(nk-1, -1, -1):
        Fi[k, :] = Fi[k+1, :] + imb[k, :]/h
    return [Fi.ravel(), np.zeros((N-1)*N)]


def main(scans, N=25):
    print('The twin against OpenFOAM')
    rel = []
    for scan in scans:
        for case in tbc.cases(scan):
            r = tbc.firstStep(case, N)
            tw = firstStep(N=N, thetaDeg=r['theta'], zw=r['zw'])
            rel.append(abs(tw['magU'].max() - r['magU'].max())/r['magU'].max())
    rel = np.array(rel)
    print(f'  max|U| over {len(rel)} runs: max relative difference {rel.max():.1e}, '
          f'median {np.median(rel):.1e}')
    print('\nThe velocity depends on the face force only through the loop imbalance')
    for th in (5.0, 12.0):
        r = firstStep(N=N, thetaDeg=th)
        alt = firstStep(N=N, thetaDeg=th,
                        faceForce=staircaseFaceForce(r['imbalance'], r['h']))
        d = np.abs(alt['magU'] - r['magU']).max()/r['magU'].max()
        print(f'  theta = {th:5.1f} deg: a staircase face force with the same vertex '
              f'loop imbalances gives the same |U| to {d:.1e} relative')


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--validate', nargs='+', metavar='SCAN',
                    default=['tiltScan_1deltaT'],
                    help='parmScanner scan directories of one-step interFoam runs')
    ap.add_argument('-N', type=int, default=25, help='cells per side (default 25)')
    a = ap.parse_args()
    main(a.validate, a.N)
