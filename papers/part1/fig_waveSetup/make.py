#!/usr/bin/env pvbatch
#
# Renders caseIllustration.pvsm (mesh and stream function wave) to a PNG.
#
# The state file shows three meshes side by side, each in one third of the
# domain, with the water volume fraction, the free surface and velocity
# arrows. Its readers are pointed at the given cases, which must contain the
# selected time. The 0 directory only holds the placeholder fields restored
# from 0.orig, so the default time is the first written one, t = 0.2 s.
#
# Part I, Figure fig:waveSetup. Input: run/progressiveWave/solvers_maxCo cases 0002, 0014,
# 0026 at t = 0.2 s. Output: waveSetup.png (transparent background).
#
# Usage:
#   pvbatch --force-offscreen-rendering make.py

import argparse
import os

from paraview.simple import (GetActiveViewOrCreate, GetAnimationScene, GetSources, Hide,
                             LoadState, SaveScreenshot)

HERE = os.path.dirname(os.path.abspath(__file__))
READER_IDS = ['10087', '10485', '10755']  # Reader proxy ids in caseIllustration.pvsm, left to right

parser = argparse.ArgumentParser(description='Render the case illustration state file.')
parser.add_argument('--scanDir', default=os.path.join(HERE, '..', '..', '..', 'run', 'progressiveWave', 'solvers_maxCo'), help='Study directory with the cases')
parser.add_argument('--cases', nargs=3, default=['0002', '0014', '0026'],
                    help='Cases for the left, middle and right mesh (nz 20, 40, 80 by default)')
parser.add_argument('--time', type=float, default=0.2, help='Time to show. Default: 0.2')
parser.add_argument('--out', default=os.path.join(HERE, 'waveSetup.png'), help='Output PNG file')
parser.add_argument('--opaque', action='store_true', help='Opaque background (default transparent)')
args = parser.parse_args()

scanDir = os.path.abspath(os.path.join(HERE, args.scanDir))
LoadState(os.path.join(HERE, 'caseIllustration.pvsm'),
          filenames=[{'id': rid, 'FileName': os.path.join(scanDir, case, 'periodicWaveTank.foam')}
                     for rid, case in zip(READER_IDS, args.cases)])

view = GetActiveViewOrCreate('RenderView')
scene = GetAnimationScene()
scene.UpdateAnimationUsingDataTimeSteps()
scene.AnimationTime = args.time
view.ViewTime = args.time

# The time annotation in the state does not follow the selected time
for (name, _), source in GetSources().items():
    if name == 'AnnotateTime1':
        Hide(source, view)

SaveScreenshot(args.out, view, ImageResolution=[2100, 1149],
               TransparentBackground=0 if args.opaque else 1)

# Crop the empty margins if Pillow is available
try:
    from PIL import Image, ImageChops
    image = Image.open(args.out)
    if image.mode == 'RGBA':
        box = image.getchannel('A').getbbox()
    else:
        image = image.convert('RGB')
        box = ImageChops.difference(image, Image.new('RGB', image.size, (255, 255, 255))).getbbox()
    if box:
        pad = 12
        image.crop((max(0, box[0] - pad), max(0, box[1] - pad),
                    min(image.width, box[2] + pad), min(image.height, box[3] + pad))).save(args.out)
except ImportError:
    pass
print('Wrote', args.out)
