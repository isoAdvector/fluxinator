# Render the velocity snapshot of one standing-wave case with ParaView.
#
# The look reproduces baseCase/state.pvsm: cells with white edges coloured by
# alpha.water (air black, water blue), white arrows oriented and scaled by the
# cell values of U (placed at the cell centres, no interpolation),
# and the alpha = 0.5 contour as a white line. The arrows use the same glyph
# sampling as the state (bounds-based uniform distribution, 5000 points, fixed
# seed), so the three solvers are sampled at identical points, and one common
# scale factor, so arrow lengths compare directly across the panels.
#
# Usage:
#   pvbatch snapshotImage.py <case dir> <time> <output png>
#           [--scale 0.1] [--zmin 0.36] [--zmax 0.80] [--width 2400]
#           [--reference]
#
# --reference adds a 1 m/s reference arrow in the water near the lower left.

import argparse
import os

from paraview.simple import *

paraview.simple._DisableFirstRenderCameraReset()

p = argparse.ArgumentParser()
p.add_argument("case")
p.add_argument("time", type=float)
p.add_argument("out")
p.add_argument("--scale", type=float, default=0.1,
               help="arrow length per unit velocity [m per m/s]")
p.add_argument("--zmin", type=float, default=0.36)
p.add_argument("--zmax", type=float, default=0.80)
p.add_argument("--width", type=int, default=2400)
p.add_argument("--reference", action="store_true")
args = p.parse_args()

foam = os.path.join(os.path.abspath(args.case), "case.foam")
open(foam, "a").close()

reader = OpenFOAMReader(FileName=foam)
reader.MeshRegions = ["internalMesh"]
reader.CellArrays = ["U", "alpha.water"]
reader.Createcelltopointfiltereddata = 1
reader.UpdatePipeline(args.time)

view = GetActiveViewOrCreate("RenderView")
height = int(round(args.width*(args.zmax - args.zmin)))
view.ViewSize = [args.width, height]
view.Background = [1, 1, 1]
view.UseColorPaletteForBackground = 0
view.OrientationAxesVisibility = 0
view.ViewTime = args.time

# Cells coloured by alpha.water, white edges (state.pvsm)
alphaLUT = GetColorTransferFunction("alpha.water")
alphaLUT.ColorSpace = "RGB"
alphaLUT.RGBPoints = [0.0, 0.0, 0.0, 0.0,
                      0.5, 0.865003, 0.865003, 0.865003,
                      1.0, 0.0, 0.6666666666666666, 1.0]
alphaLUT.AutomaticRescaleRangeMode = "Never"

meshDisplay = Show(reader, view)
meshDisplay.SetRepresentationType("Surface With Edges")
ColorBy(meshDisplay, ("CELLS", "alpha.water"))
meshDisplay.LookupTable = alphaLUT
meshDisplay.EdgeColor = [1.0, 1.0, 1.0]
meshDisplay.EdgeOpacity = 0.4
meshDisplay.LineWidth = 1.5
meshDisplay.SetScalarBarVisibility(view, False)

# Interface: alpha = 0.5 contour as a white line
contour = Contour(Input=reader)
contour.ContourBy = ["POINTS", "alpha.water"]
contour.Isosurfaces = [0.5]
contourDisplay = Show(contour, view)
contourDisplay.ColorArrayName = ["POINTS", ""]
contourDisplay.SetRepresentationType("Wireframe")
contourDisplay.AmbientColor = [1.0, 1.0, 1.0]
contourDisplay.DiffuseColor = [1.0, 1.0, 1.0]
contourDisplay.LineWidth = 4.5

# Velocity arrows (state.pvsm glyph settings, common scale factor)
glyph = Glyph(Input=reader, GlyphType="Arrow")
glyph.OrientationArray = ["CELLS", "U"]
glyph.ScaleArray = ["CELLS", "U"]
glyph.VectorScaleMode = "Scale by Magnitude"
glyph.ScaleFactor = args.scale
glyph.GlyphMode = "Uniform Spatial Distribution (Bounds Based)"
glyph.MaximumNumberOfSamplePoints = 5000
glyph.Seed = 10339
glyph.GlyphType.TipResolution = 6
glyph.GlyphType.TipRadius = 0.1
glyph.GlyphType.TipLength = 0.35
glyph.GlyphType.ShaftResolution = 6
glyph.GlyphType.ShaftRadius = 0.03
glyphDisplay = Show(glyph, view)
glyphDisplay.ColorArrayName = ["POINTS", ""]
glyphDisplay.AmbientColor = [1.0, 1.0, 1.0]
glyphDisplay.DiffuseColor = [1.0, 1.0, 1.0]
# The cell centres lie mid-way through the one-cell-thick slab, behind its front
# face; move the arrows towards the camera as state.pvsm does
glyphDisplay.Translation = [0.0, -1.0, 0.0]

if args.reference:
    # A 1 m/s arrow, in front of the mesh, in the still water at lower left
    z0 = args.zmin + 0.06*(args.zmax - args.zmin)
    ref = PointSource(Center=[0.03, -0.006, z0], NumberOfPoints=1, Radius=0.0)
    refU = Calculator(Input=ref)
    refU.AttributeType = "Point Data"
    refU.ResultArrayName = "Uref"
    refU.Function = "iHat"
    refGlyph = Glyph(Input=refU, GlyphType="Arrow")
    refGlyph.OrientationArray = ["POINTS", "Uref"]
    refGlyph.ScaleArray = ["POINTS", "Uref"]
    refGlyph.VectorScaleMode = "Scale by Magnitude"
    refGlyph.ScaleFactor = args.scale
    refGlyph.GlyphMode = "All Points"
    refGlyph.GlyphType.TipResolution = 6
    refGlyph.GlyphType.TipRadius = 0.1
    refGlyph.GlyphType.TipLength = 0.35
    refGlyph.GlyphType.ShaftResolution = 6
    refGlyph.GlyphType.ShaftRadius = 0.03
    refDisplay = Show(refGlyph, view)
    refDisplay.ColorArrayName = ["POINTS", ""]
    refDisplay.AmbientColor = [1.0, 1.0, 1.0]
    refDisplay.DiffuseColor = [1.0, 1.0, 1.0]
    refDisplay.Translation = [0.0, -1.0, 0.0]

    label = Text(Text="1 m/s")
    labelDisplay = Show(label, view)
    labelDisplay.WindowLocation = "Any Location"
    labelDisplay.Position = [(0.03 + args.scale + 0.012), 0.035]
    labelDisplay.FontSize = 72
    labelDisplay.FontFamily = "Times"
    labelDisplay.Color = [1.0, 1.0, 1.0]

# Parallel projection framing the full width and the band zmin < z < zmax
zc = 0.5*(args.zmin + args.zmax)
view.CameraParallelProjection = 1
view.CameraFocalPoint = [0.5, 0.0, zc]
view.CameraPosition = [0.5, -3.0, zc]
view.CameraViewUp = [0.0, 0.0, 1.0]
view.CameraParallelScale = 0.5*(args.zmax - args.zmin)

Render(view)
SaveScreenshot(args.out, view, ImageResolution=[3*args.width//2, 3*height//2])
print("wrote", args.out)
