# Copyright (C) 2026 Johan Roenby
# SPDX-License-Identifier: GPL-3.0-or-later
"""ParaView snapshots of the prescribed-flow cases (run/prescribedFlows/<case>) at t = 1 s:
volume fraction in the cells (white air, blue water) with the mesh lines, the prescribed
velocity as arrows at every 5th cell in each direction, and the interface alpha = 0.5 at
t = 0 (orange) and t = 1 s (black). Written fields: alpha.water from the run, U the
prescribed velocity.
Writes snap/<case>.png (input of make.py).
Usage: pvbatch snapshots.py"""
import os
from paraview.simple import *

HERE = os.path.dirname(os.path.abspath(__file__))
RUN = os.path.join(HERE, "..", "..", "..", "run", "prescribedFlows")
OUT = os.path.join(HERE, "snap")
os.makedirs(OUT, exist_ok=True)
# arrow length [m] per m/s
SCALE = {"shear": 0.2, "rotation": 0.09}
PX = 1600
ORANGE, BLACK, GREY = [0.92, 0.41, 0.20], [0.0, 0.0, 0.0], [0.15, 0.15, 0.15]

SUBSAMPLE = """
import numpy as np
from vtkmodules.vtkCommonCore import vtkPoints
from vtk.numpy_interface import dataset_adapter as dsa
inp = dsa.WrapDataObject(self.GetInput()); out = self.GetOutput()
x, z = inp.Points[:, 0], inp.Points[:, 2]
i, j = np.round(x*50 - 0.5).astype(int), np.round(z*50 - 0.5).astype(int)
keep = np.where((i % 5 == 2) & (j % 5 == 2))[0]
pts = vtkPoints()
for k in keep:
    pts.InsertNextPoint(x[k], -0.002, z[k])
out.SetPoints(pts)
dsa.WrapDataObject(out).PointData.append(inp.PointData['U'][keep], 'U')
"""


def plane(t, reader):
    """y = 0 slice of the case at time t, as one polydata."""
    s = Slice(Input=MergeBlocks(Input=ForceTime(Input=reader, ForcedTime=t)))
    s.SliceType = "Plane"
    s.SliceType.Origin = [0.5, 0.0, 0.5]
    s.SliceType.Normal = [0.0, 1.0, 0.0]
    s.Triangulatetheslice = 0
    return s


for case, scale in SCALE.items():
    r = OpenFOAMReader(FileName=os.path.join(RUN, case, "case.foam"), MeshRegions=["internalMesh"],
                       CellArrays=["alpha.water", "U"])
    # the reader skips the 0 directory by default, and builds its time list when the file
    # name is set: switching SkipZeroTime off takes effect only after a rescan, which a
    # change of the case type forces
    r.SkipZeroTime = 0
    r.CaseType = "Decomposed Case"
    r.UpdatePipelineInformation()
    r.CaseType = "Reconstructed Case"
    r.UpdatePipelineInformation()
    assert 0.0 in r.TimestepValues, r.TimestepValues
    view = CreateRenderView(ViewSize=[PX, PX])
    view.OrientationAxesVisibility = 0
    view.UseColorPaletteForBackground = 0
    view.Background = [1, 1, 1]
    view.AxesGrid.Visibility = 0

    p1 = plane(1.0, r)
    d = Show(p1, view)
    ColorBy(d, ("CELLS", "alpha.water"))
    lut = GetColorTransferFunction("alpha.water")
    lut.AutomaticRescaleRangeMode = "Never"
    lut.RGBPoints = [0.0, 1.0, 1.0, 1.0, 1.0, 0.45, 0.66, 0.90]
    lut.ColorSpace = "RGB"
    d.SetScalarBarVisibility(view, False)
    d.SetRepresentationType("Surface With Edges")
    d.EdgeColor = [0.80, 0.80, 0.80]
    d.Ambient, d.Diffuse = 1.0, 0.0

    for t, col in ((0.0, ORANGE), (1.0, BLACK)):
        c = Contour(Input=CellDatatoPointData(Input=plane(t, r)), ContourBy=["POINTS", "alpha.water"],
                    Isosurfaces=[0.5])
        cd = Show(c, view)
        cd.ColorArrayName = [None, ""]
        cd.AmbientColor = cd.DiffuseColor = col
        cd.LineWidth = 7.0

    pf = ProgrammableFilter(Input=CellCenters(Input=p1))
    pf.OutputDataSetType = "vtkPolyData"
    pf.Script = SUBSAMPLE
    g = Glyph(Input=pf, GlyphType="Arrow", OrientationArray=["POINTS", "U"],
              ScaleArray=["POINTS", "U"], ScaleFactor=scale, GlyphMode="All Points")
    g.VectorScaleMode = "Scale by Magnitude"
    g.GlyphType.TipResolution = 16
    g.GlyphType.ShaftResolution = 8
    g.GlyphType.TipLength = 0.3
    g.GlyphType.TipRadius = 0.12
    g.GlyphType.ShaftRadius = 0.035
    gd = Show(g, view)
    gd.ColorArrayName = [None, ""]
    gd.AmbientColor = gd.DiffuseColor = GREY
    gd.Ambient, gd.Diffuse = 1.0, 0.0

    Render(view)
    view.CameraPosition = [0.5, -3.0, 0.5]
    view.CameraFocalPoint = [0.5, 0.0, 0.5]
    view.CameraViewUp = [0.0, 0.0, 1.0]
    view.CameraParallelProjection = 1
    view.CameraParallelScale = 0.5
    Render(view)
    SaveScreenshot(os.path.join(OUT, case + ".png"), view, ImageResolution=[PX, PX])
    print("written", case)
    for s in list(GetSources().values()):
        Delete(s)
    Delete(view)
