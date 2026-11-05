# Translating Taylor-Green vortex (Part II, Section 3)

Single phase, doubly periodic, 64 x 64 cells on [0, 2 pi]^2:
u = U0 + sin(x - U0 t) cos z, w = -cos(x - U0 t) sin z, U0 = 1 (file `U0` in each
case; 0 in `implicitNoTranslation_dt0.01`). The exact x-mean of u stays U0; a time
lag theta dt between the mass flux and the convected velocity gives the spurious
mean flow -theta dt U0 t cos(2z)/2 (file `theta`: 1 implicit, 1/2 Crank-Nicolson,
0 explicit).

Each case is `baseCase` with its own `system/` settings; its `Allrun` writes the
initial velocity with `../writeU.py` and runs interIsoGoam. `writeData.py <case> ...`
prints the measured and predicted mean-flow amplitude (tgv.txt, tgv.png). The paper
figure is made in `papers/part2/fig_lagTGV`.
