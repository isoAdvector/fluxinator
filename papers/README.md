# Figures and tables of the papers

One directory per figure or table, named after its LaTeX label: `fig_X/` for
`fig:X` and `tab_X/` for `tab:X`. Each holds

- `make.py` (or `make.sh`) that makes the figure or table from the runs in
  `../run`, with its inputs listed at the top. Run it with `python3 make.py`;
  ParaView figures (`fig_reach`, `fig_waveSetup`, the snapshots of
  `fig_rigCasesFlows`) with `pvbatch --force-offscreen-rendering make.py`.
- the file the paper shows: `<X>.pdf` (and `.png`) for a figure, `<X>.txt`
  with the numbers for a table.

The caption of each figure and table in the papers names its directory.

Scripts used by several figures or tables are in `part1/` and `part2/` next to
the directories (for example `part2/wave.py`, which makes `tab_wave` and
`fig_waveCrest`); each output's `make.py` runs them. `part1/common.py` and
`part2/common.py` hold the paths, the plotting style and, for Part II, the
routing of each output to its directory. `part2/data/` holds data that several
figures share or that is cited in the text (`browser.json`: metrics of all
fluxinator and Part I wave runs, `python3 browserData.py`; `rig.txt`,
`lag.txt`, `paperSelection.txt`, the crest-model files).

`part2/fig_lagWave` is no longer in the Part II paper (the paragraph on the
spurious water current at depth was removed on 2026-11-05). It is kept because
`part2/lag.py` writes it together with `fig_lagTGV`.
