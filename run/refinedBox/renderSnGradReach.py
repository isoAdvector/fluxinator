#!/usr/bin/env pvbatch
#
# Renders the paper figure panels reach_uncorrected.png and reach_corrected.png
# for the untilted box with a refined region (z_w = 0.49): cells adjacent to faces with a
# non-zero uncorrected and corrected snGrad(rho) in the initial field, coloured
# by ring level around the interface cells (0 amber, 1 green, 2 red with black
# outline), with the alpha = 0.5 contour. Requires the fields written by
# postProcess -time 0 -dict ../snGradRhoReach in scan/g1_zero_z0.4900.
#
# Usage (from this directory):
#   pvbatch --force-offscreen-rendering renderSnGradReach.py [--outDir paperFigures]

import argparse
import os

from paraview.simple import *

HERE = os.path.dirname(os.path.abspath(__file__))
CASE = "scan/g1_zero_z0.4900"

parser = argparse.ArgumentParser()
parser.add_argument("--outDir", default="paperFigures")
args = parser.parse_args()
outDir = os.path.join(HERE, args.outDir)
os.makedirs(outDir, exist_ok=True)

# (field, output file)
PANELS = [
    ("snGradReachUncorrected", "reach_uncorrected.png"),
    ("snGradReachCorrected", "reach_corrected.png"),
]

RING_COLOURS = {
    "-1": [1.0, 1.0, 1.0],
    "0": [0.95, 0.76, 0.31],
    "1": [0.60, 0.82, 0.71],
    "2": [0.84, 0.15, 0.16],
}

SIZE = [1400, 560]


def solidColour(display, colour):
    """Turn off scalar colouring and draw in a single colour."""
    display.ColorArrayName = ["POINTS", ""]
    display.AmbientColor = colour
    display.DiffuseColor = colour


def render(field, outFile):
    view = CreateView("RenderView")
    view.ViewSize = SIZE
    view.Background = [1, 1, 1]
    view.UseColorPaletteForBackground = 0
    view.OrientationAxesVisibility = 0
    # Hide the legend grid (axis labels along the image borders)
    view.LegendGrid.Visibility = 0
    view.LegendGrid.GridVisibility = 0
    view.LegendGrid.LegendLabelOpacity = 0

    foam = os.path.join(HERE, CASE, os.path.basename(CASE) + ".foam")
    reader = OpenFOAMReader(FileName=foam)
    reader.MeshRegions = ["internalMesh"]
    # Keep polyhedra whole so that their faces are drawn without extra edges
    # (the property name differs between ParaView versions)
    for name in ("Decomposepolyhedra", "DecomposePolyhedra"):
        if name in reader.ListProperties():
            setattr(reader, name, 0)
    reader.Createcelltopointfiltereddata = 1
    # The zero time, which holds the diagnostic fields, is only listed after
    # the files are reloaded with SkipZeroTime off
    reader.SkipZeroTime = 0
    ReloadFiles(reader)
    reader.UpdatePipelineInformation()
    reader.CellArrays = ["alpha.water", field]
    reader.UpdatePipeline(0.0)
    view.ViewTime = 0.0

    display = Show(reader, view)
    # No axes grid around the data (it labels the image borders)
    for owner in (display, view):
        for name in ("DataAxesGrid", "AxesGrid"):
            if name in owner.ListProperties():
                grid = getattr(owner, name)
                for prop in ("Visibility", "GridAxesVisibility"):
                    if prop in grid.ListProperties():
                        setattr(grid, prop, 0)
    display.SetRepresentationType("Surface With Edges")
    display.EdgeColor = [0.6, 0.6, 0.6]
    # Flat lighting so that the colours appear as specified
    display.Ambient = 1.0
    display.Diffuse = 0.0

    ColorBy(display, ("CELLS", field))
    lut = GetColorTransferFunction(field)
    lut.InterpretValuesAsCategories = 1
    annotations, colours = [], []
    for value, colour in RING_COLOURS.items():
        annotations += [value, value]
        colours += colour
    lut.Annotations = annotations
    lut.IndexedColors = colours
    display.SetScalarBarVisibility(view, False)

    ring2 = Threshold(Input=reader)
    ring2.Scalars = ["CELLS", field]
    ring2.LowerThreshold = 1.5
    ring2.UpperThreshold = 2.5
    ring2.ThresholdMethod = "Between"
    ring2Display = Show(ring2, view)
    ring2Display.SetRepresentationType("Wireframe")
    solidColour(ring2Display, [0, 0, 0])
    ring2Display.LineWidth = 5

    # The contour of the one cell thick mesh is seen edge-on, so draw its edges
    contour = Contour(Input=reader)
    contour.ContourBy = ["POINTS", "alpha.water"]
    contour.Isosurfaces = [0.5]
    contourDisplay = Show(contour, view)
    contourDisplay.SetRepresentationType("Wireframe")
    solidColour(contourDisplay, [0, 0, 0])
    contourDisplay.LineWidth = 3

    view.CameraPosition = [0.5, -3.0, 0.5]
    view.CameraFocalPoint = [0.5, 0.0, 0.5]
    view.CameraViewUp = [0.0, 0.0, 1.0]
    view.CameraParallelProjection = 1
    view.CameraParallelScale = 0.23

    out = os.path.join(outDir, outFile)
    SaveScreenshot(out, view, ImageResolution=SIZE, TransparentBackground=0)
    print("Wrote", out)

    for source in (contour, ring2, reader):
        Delete(source)
    Delete(view)


for panel in PANELS:
    render(*panel)
