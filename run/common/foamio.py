"""Minimal OpenFOAM ASCII reader: points, faces, owner, neighbour, vol fields."""
import re, numpy as np

def _body(path):
    txt = open(path).read()
    # strip header block FoamFile { ... } and comments
    txt = re.sub(r'/\*.*?\*/', ' ', txt, flags=re.S)
    txt = re.sub(r'//.*', ' ', txt)
    i = txt.find('FoamFile')
    j = txt.find('}', i)
    return txt[j+1:]

def read_points(path):
    b = _body(path)
    m = re.search(r'(\d+)\s*\(', b)
    n = int(m.group(1))
    vals = np.array(re.findall(r'\(([^()]*)\)', b[m.end()-1:])[:n])
    pts = np.array([[float(x) for x in v.split()] for v in vals])
    assert pts.shape == (n, 3), pts.shape
    return pts

def read_faces(path):
    b = _body(path)
    m = re.search(r'(\d+)\s*\(', b)
    n = int(m.group(1))
    grp = re.findall(r'(\d+)\s*\(([^()]*)\)', b[m.end()-1:])[:n]
    faces = [np.array([int(x) for x in g[1].split()]) for g in grp]
    assert len(faces) == n
    return faces

def read_labels(path):
    b = _body(path)
    m = re.search(r'(\d+)\s*\(', b)
    n = int(m.group(1))
    vals = re.findall(r'-?\d+', b[m.end():])[:n]
    a = np.array([int(v) for v in vals])
    assert a.size == n
    return a

def read_vol(path):
    """Return (internalField array, dict of boundary patch entries)."""
    b = _body(path)
    i = b.find('internalField')
    j = b.find('boundaryField')
    s = b[i:j]
    if 'nonuniform' not in s.split('\n')[0]:
        raise ValueError('uniform internalField in ' + path)
    m = re.search(r'(\d+)\s*\(', s)
    n = int(m.group(1))
    tail = s[m.end()-1:]
    if 'vector' in s[:m.start()] or 'Vector' in s[:m.start()]:
        vals = re.findall(r'\(([^()]*)\)', tail)[:n]
        arr = np.array([[float(x) for x in v.split()] for v in vals])
    else:
        inner = tail[tail.find('(')+1:tail.rfind(')')]
        arr = np.array([float(x) for x in inner.split()])
    assert arr.shape[0] == n, (arr.shape, n, path)
    return arr

def read_boundary(path):
    """constant/polyMesh/boundary -> list of (name, nFaces, startFace)."""
    b = _body(path)
    out = []
    for m in re.finditer(r'(\w+)\s*\{([^{}]*)\}', b):
        d = m.group(2)
        nf = re.search(r'nFaces\s+(\d+)', d)
        sf = re.search(r'startFace\s+(\d+)', d)
        if nf and sf:
            out.append((m.group(1), int(nf.group(1)), int(sf.group(1))))
    return out
