#!/usr/bin/env python3
#
# Sets a hydrostatic initial p_rgh in 0/p_rgh, consistent with the initial
# volume fraction written by setAlphaField and a flat interface through the
# origin in system/setAlphaFieldDict:
#
#   p_rgh = p - rho*(g & x) = -rho*(g & x_I)   (p = 0 at the interface, hRef = 0)
#
# with rho = alpha*rho_water + (1 - alpha)*rho_air per cell. Only the internal
# field is written; fixedFluxPressure patches without a value take the adjacent
# cell values when read. Run from the case directory after setAlphaField.

import re
import sys


def read(path):
    with open(path) as f:
        return f.read()


def vector(text, key):
    m = re.search(r'\b%s\s+(?:\[[^\]]*\]\s*)?\(\s*([-\d.eE+]+)\s+([-\d.eE+]+)\s+([-\d.eE+]+)\s*\)' % key, text)
    return tuple(float(v) for v in m.groups())


def phase_rho(text, phase):
    m = re.search(r'\b%s\s*\{[^}]*\brho\s+([-\d.eE+]+)\s*;' % phase, text)
    return float(m.group(1))


g = vector(read('constant/g'), 'value')
origin = vector(read('system/setAlphaFieldDict'), 'origin')
transport = read('constant/transportProperties')
rho_water, rho_air = phase_rho(transport, 'water'), phase_rho(transport, 'air')
gxI = sum(gi*oi for gi, oi in zip(g, origin))

alpha_text = read('0/alpha.water')
m = re.search(r'internalField\s+nonuniform\s+List<scalar>\s*(\d+)\s*\(', alpha_text)
if not m:
    sys.exit('0/alpha.water has no nonuniform internalField; run setAlphaField first')
n = int(m.group(1))
alpha = [float(v) for v in alpha_text[m.end():alpha_text.index(')', m.end())].split()[:n]]

prgh = [-(a*rho_water + (1 - a)*rho_air)*gxI for a in alpha]
field = 'internalField   nonuniform List<scalar>\n%d\n(\n%s\n)\n;' % (n, '\n'.join('%.14g' % v for v in prgh))

prgh_text = read('0/p_rgh')
new_text, count = re.subn(r'internalField\s+uniform\s+[-\d.eE+]+\s*;', lambda _: field, prgh_text, count=1)
if count != 1:
    sys.exit('0/p_rgh: could not find a uniform internalField to replace')
with open('0/p_rgh', 'w') as f:
    f.write(new_text)

print('g & x_I = %g, rho_water = %g, rho_air = %g, cells = %d' % (gxI, rho_water, rho_air, n))
print('p_rgh range: %g to %g' % (min(prgh), max(prgh)))
