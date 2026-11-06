#!/usr/bin/env python3
# Copyright (C) 2026 Johan Roenby
# SPDX-License-Identifier: GPL-3.0-or-later
#
# Reading a tilted box case: the ascii fields, the mesh vertices of the 2-D
# plane, the tilt angle and the water level it was initialised with.
#
# The mesh is one blockMesh hex block of N x 1 x N cells, so the points are
# ordered with x fastest, then y, then z, and the N+1 by N+1 grid of vertices
# in the front plane is points[ix + (N+1)*2*iz].

import os
import re

import numpy as np


def _body(path):
    txt = open(path).read()
    txt = re.sub(r'//.*', ' ', re.sub(r'/\*.*?\*/', ' ', txt, flags=re.S))
    return txt[txt.find('}', txt.find('FoamFile')) + 1:]


def readVol(path):
    """internalField of an ascii volScalarField or volVectorField."""
    s = _body(path)
    s = s[s.find('internalField'):s.find('boundaryField')]
    m = re.search(r'(\d+)\s*\(', s)
    n = int(m.group(1))
    tail = s[m.end()-1:]
    if 'vector' in s[:m.start()].lower():
        v = re.findall(r'\(([^()]*)\)', tail)[:n]
        return np.array([[float(x) for x in w.split()] for w in v])
    inner = tail[tail.find('(')+1:tail.rfind(')')]
    return np.array([float(x) for x in inner.split()])[:n]


def readPoints(path):
    s = _body(path)
    m = re.search(r'(\d+)\s*\(', s)
    n = int(m.group(1))
    v = re.findall(r'\(([^()]*)\)', s[m.end()-1:])[:n]
    return np.array([[float(x) for x in w.split()] for w in v])


def vertexGrid(case, N=25):
    """(N+1, N+1, 2) array of the (x, z) mesh vertices, indexed [K, I]."""
    p = readPoints(os.path.join(case, 'constant', 'polyMesh', 'points'))
    idx = np.arange(N+1)[None, :] + (N+1)*2*np.arange(N+1)[:, None]
    return p[idx][..., [0, 2]]


def cellCentres(case, N=25):
    """(N, N, 2) array of the (x, z) cell centres, indexed [k, i]."""
    V = vertexGrid(case, N)
    return 0.25*(V[:-1, :-1] + V[:-1, 1:] + V[1:, :-1] + V[1:, 1:])


def _entry(path, key):
    if not os.path.isfile(path):
        return None
    m = re.search(rf'\b{key}\s+([^;]+);', _body(path))
    return m.group(1).strip() if m else None


def tiltAngle(case):
    v = _entry(os.path.join(case, 'rotAngleDict'), 'rotAngle')
    if v is None:
        raise RuntimeError('no rotAngle in ' + case)
    return float(v)


def waterLevel(case):
    """z of the plane setAlphaField was given (setAlphaFieldDict origin)."""
    v = _entry(os.path.join(case, 'system', 'setAlphaFieldDict'), 'origin')
    if v is None:
        raise RuntimeError('no origin in ' + case)
    return float(v.strip('()').split()[2])


def cases(scanDir):
    """The case directories of a parmScanner scan that ran, in scan order."""
    out = []
    for name in sorted(os.listdir(scanDir)):
        d = os.path.join(scanDir, name)
        if os.path.isdir(d) and os.path.isfile(os.path.join(d, '0.01', 'U')):
            out.append(d)
    return out


def firstStep(case, N=25):
    """Tilt, water level, alpha and |U| after the first time step."""
    U = readVol(os.path.join(case, '0.01', 'U'))
    return dict(case=case, theta=tiltAngle(case), zw=waterLevel(case),
                alpha=readVol(os.path.join(case, '0', 'alpha.water')).reshape(N, N),
                magU=np.linalg.norm(U, axis=1).reshape(N, N),
                C=cellCentres(case, N), V=vertexGrid(case, N))
