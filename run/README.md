# Test cases

| Case | Content | Papers |
|---|---|---|
| `tiltedBox` | flat interface at rest in a tilted box | Part I Sec. 3 |
| `skewedBox` | flat interface at rest on a skewed (non-orthogonal) mesh | Part I Sec. 4 |
| `refinedBox` | flat interface at rest in a box with a refined region | Part I Sec. 4 |
| `standingWave` | standing wave across a refinement boundary | Part I Sec. 5 |
| `progressiveWave` | stream function wave in a periodic domain of one wavelength | Part I Secs. 2, 6; Part II Secs. 1, 3, 6, 7 |
| `taylorGreenVortex` | translating Taylor-Green vortex, single phase (the time lag) | Part II Sec. 3 |
| `prescribedFlows` | water disc in air, prescribed velocity, momentum equation only | Part II Secs. 4, 5 |
| `periodicDisc` | water disc in a periodic box, full solver, uniform acceleration or no force | Part II Secs. 4, 5 |
| `quickTest` | fluxinator on the accelerated disc and one wave period, about a minute | |
| `common` | shared Python modules and helper scripts | |

Every case has a README. A study directory holds the runs of one parameter
scan or set of variants; `<study>_ddtCorrOff` repeats a Part I study with
PIMPLE `ddtCorr no` (Part I, Appendix A; `common/runDdtCorrOffStudy.sh`).
The figures and tables of the papers are made in `../papers`.

Runs are not tracked by git, only their set-ups and the summaries in `data/`;
the run data behind the papers is published separately.
