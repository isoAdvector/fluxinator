"""Common paths, plotting style and output routing of the Part II figures and tables.

Every figure or table of Part II with LaTeX label fig:X or tab:X has its own
directory fig_X/ or tab_X/ here, with a make.py that writes the file the paper
shows into that directory. Scripts that make several outputs live in this
directory and are run by the make.py of each of their outputs; save() and
write() route each output by name to its directory (FIGDIR, TABDIR). Data that
several figures or tables share, and intermediate data, go to data/."""
import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.font_manager  # noqa: E402,F401

HERE = os.path.dirname(os.path.abspath(__file__))                   # papers/part2
G = os.path.abspath(os.path.join(HERE, "..", ".."))                 # repository root
RUN = os.path.join(G, "run")
WAVE = os.path.join(RUN, "progressiveWave")                         # Part I and final-solver studies
CONV = os.path.join(WAVE, "convection")                             # convection variants of Part II
RIG = os.path.join(RUN, "prescribedFlows")
DISC = os.path.join(RUN, "periodicDisc")
TGV = os.path.join(RUN, "taylorGreenVortex")
sys.path.insert(0, os.path.join(RUN, "common"))
import kinematics as K  # noqa: E402,F401

DATA = os.path.join(HERE, "data")
# figure name -> (directory, file name); data file name -> directory
FIGDIR = {"problem": ("fig_problem", "problem"), "mechanismLag": ("fig_mechanismLag", "mechanismLag"),
          "lagTGV": ("fig_lagTGV", "lagTGV"), "lagWave": ("fig_lagWave", "lagWave"),
          "rigCasesFlows": ("fig_rigCasesFlows", "rigCasesFlows"), "rigErrors": ("fig_rigErrors", "rigErrors"),
          "fullSolverUerr": ("fig_fullSolver", "fullSolver"), "mechanismFace": ("fig_mechanismFace", "mechanismFace"),
          "momFlux": ("fig_momFlux", "momFlux"), "waveField": ("fig_waveField", "waveField"),
          "waveCrest": ("fig_waveCrest", "waveCrest"), "lwCheck": ("fig_lwCheck", "lwCheck"),
          "selCase": ("fig_selCase", "selCase"), "selHeight": ("fig_selHeight", "selHeight"),
          "selCost": ("fig_selCost", "selCost"), "finalOverview": ("fig_finalOverview", "finalOverview")}
TABDIR = {"wave.txt": "tab_wave", "fullSolver.txt": "tab_fullSolver", "rigConvergence.txt": "tab_rigConvergence",
          "gain.txt": "fig_mechanismFace", "momFlux.txt": "fig_momFlux", "lwCheck.txt": "fig_lwCheck",
          "finalSolver.txt": "fig_finalOverview", "finalSolver.csv": "fig_finalOverview"}


def outpath(name):
    """Where a data file of this name is written (and read back)."""
    return os.path.join(HERE, TABDIR[name], name) if name in TABDIR else os.path.join(DATA, name)
# Latin Modern Roman (the OpenType Computer Modern of the paper text) and Computer Modern math
LMDIR = "/usr/share/texmf/fonts/opentype/public/lm"
for _f in ("lmroman10-regular", "lmroman10-bold", "lmroman10-italic", "lmroman10-bolditalic"):
    if os.path.exists(os.path.join(LMDIR, _f + ".otf")):
        matplotlib.font_manager.fontManager.addfont(os.path.join(LMDIR, _f + ".otf"))
STYLE = {"font.family": "serif", "font.serif": ["Latin Modern Roman", "DejaVu Serif"],
         "mathtext.fontset": "cm", "font.size": 9, "axes.linewidth": 0.7,
         "lines.linewidth": 1.2, "lines.markersize": 3.5, "legend.fontsize": 8,
         "legend.frameon": False, "axes.grid": True, "grid.linewidth": 0.4,
         "grid.alpha": 0.5}
plt.rcParams.update(STYLE)
# Figure widths in inches: FULL is the text width of the paper (455 pt), so figures are not scaled
FULL, HALF = 6.3, 3.1
# Variant colours, fixed per variant across all figures
COL = {"implicit": "#eb6834", "Crank-Nicolson": "#a35c9e", "explicit": "#2a78d6",
       "explicit + LW": "#86b6ef", "explicit + interfaceUpwind": "#1baf7a",
       "explicit + LW + interfaceUpwind": "#0f6e4a", "exact": "#1f1f1e"}


LOWER = ("inter", "fluxinator", "limitedLinear", "linear", "nOuter", "ddtCorr", "a1", "a2", "a3", "b2", "b3", "c2", "c3", "c5", "e3")


def capitalise(s):
    """Capital first letter of a figure text, keeping variables (single letters such as t, z),
    math, solver and keyword names (interFoam, limitedLinearV, ...) and panel prefixes '(a) '."""
    pre = ""
    if len(s) > 4 and s[0] == "(" and s[2] == ")" and s[3] == " ":
        pre, s = s[:4], s[4:]
    word = s.split(" ")[0].split("\n")[0]
    if not s or not s[0].islower() or len(word.rstrip(",.:;")) == 1 or s.startswith(LOWER):
        return pre + s
    return pre + s[0].upper() + s[1:]


def capitalise_texts(fig):
    for t in fig.findobj(matplotlib.text.Text):
        s = t.get_text()
        if s:
            t.set_text(capitalise(s))


def save(fig, name):
    capitalise_texts(fig)
    d, f = FIGDIR[name]
    os.makedirs(os.path.join(HERE, d), exist_ok=True)
    fig.savefig(os.path.join(HERE, d, f + ".pdf"), dpi=800, bbox_inches="tight", pad_inches=0.02)
    fig.savefig(os.path.join(HERE, d, f + ".png"), dpi=150, bbox_inches="tight", pad_inches=0.02)


def write(name, lines):
    p = outpath(name)
    os.makedirs(os.path.dirname(p), exist_ok=True)
    open(p, "w").write("\n".join(lines) + "\n")
    print("\n".join(lines))


def run_script(script, *args):
    """Run a shared script of this directory as __main__ (used by the make.py files)."""
    import runpy
    sys.argv = [script] + list(args)
    sys.path.insert(0, HERE)
    runpy.run_path(os.path.join(HERE, script), run_name="__main__")
