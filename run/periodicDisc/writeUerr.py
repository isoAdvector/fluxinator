"""Write Uerr = U - (U0 + a t) (the spurious velocity) into every time directory of the
full-solver uniform-acceleration cases, for viewing in ParaView. Usage: python3 writeUerr.py <case> ..."""
import os, re, sys, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'common'))
import kinematics as K
U0, a = np.array([0.1, 0, 0.05]), np.array([0.5, 0, 0])
for c in sys.argv[1:]:
    open(f'{c}/case.foam', 'w').close()
    for t in K.times(c):
        if float(t) == 0 or not os.path.exists(f'{c}/{t}/U'): continue
        e = K.foamio.read_vol(f'{c}/{t}/U') - (U0 + a*float(t))
        txt = open(f'{c}/{t}/U').read()
        head = txt[:txt.index('internalField')].replace('object      U;', 'object      Uerr;')
        bnd = txt[txt.index('boundaryField'):]
        body = f'internalField   nonuniform List<vector>\n{len(e)}\n(\n' + '\n'.join(f'({x:.6e} {y:.6e} {z:.6e})' for x, y, z in e) + '\n)\n;\n\n'
        open(f'{c}/{t}/Uerr', 'w').write(head + body + bnd)
    print(c, 'done')
