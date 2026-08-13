# LRDs as Super-Eddington Fountain Flows — wind models (reproducibility package)

Data, generation scripts, and figure notebook for the 1-D super-Eddington
wind models (Owocki, Townsend & Quataert 2017 formulation) used in the paper.
The solver itself, `tiring_solver`, is included as a git submodule pinned to
the exact commit that produced these results.

## Getting the code

```bash
git clone --recurse-submodules <this-repo-url>
# or, after a plain clone:
git submodule update --init
```

This fetches `winds/tiring_solver/` from https://github.com/nkaaz/owocki_winds
at commit `8bf5754` (the pin is stored in this repo; see `.gitmodules`).
`winds/tiring_solver/README.md` documents the solver's method and range of
validity, and `winds/tiring_solver/validate.py` re-verifies solutions against
the governing equations.

## Layout

    winds/       solver (submodule), generation scripts, and all wind data
    analysis/    figure notebook(s)
    figures/     rendered paper figures

## Model setup

Spherically symmetric, radiation-driven wind from a central mass
M_BH = 10⁶ M☉ with κ_es = 0.34 cm²/g. The dimensionless problem
(Owocki et al. 2017, Appendix A) is two coupled ODEs for the scaled kinetic
energy w(x) = v²/v_esc² and scaled radiation pressure q(x) on
x = 1 − R/r, governed by the base Eddington ratio Γ₀ = L₀/L_Edd and the
photon-tiring number m = Ṁ G M_BH/(R L₀), with m_max = 1 − 1/Γ₀.

The heating (inner) radius follows the fountain-flow prescription
(Ṁ_out = Ṁ, h = 0.3):

    m_sph = 0.2/Γ₀   ⇒   R = Ṁ G M_BH / (0.2 L_Edd)

which ties m and Γ₀ together, so the wind family is parameterized by the
single tiring fraction f = m/m_max:

    Γ₀ = 1 + 0.2/f,   m = 0.2/Γ₀ = f·m_max.

Because R ∝ Ṁ, the dimensionless solution is **independent of Ṁ**
(self-similar): one solution per f serves every mass-loss rate, with Ṁ
entering only through the physical rescaling (R ∝ Ṁ, v_esc ∝ Ṁ^(−1/2)).

## Contents (three layers)

### 1. `winds/runs_f/f{F}.npz` — dimensionless solutions (`winds/solve_f_grid.py`)

Arrays `x, w, q, p` (4000-point uniform x-grid) plus scalars. Boundary
condition w_sg = 0 (idealized zero-sound-speed base). Every solution is
verified against the governing equations; `ode_residual` is the stored
normalized residual (see `tiring_solver/README.md`).

| f = m/m_max | Γ₀ | m | w₁ | v∞/v_esc | q₀ | ODE residual |
|---|---|---|---|---|---|---|
| 0.60 | 1.333333 | 0.150000 | 1.045e-01 | 0.323 | 4.361e-01 | 3.1e-07 |
| 0.75 | 1.266667 | 0.157895 | 4.636e-02 | 0.215 | 5.715e-01 | 1.5e-07 |
| 0.85 | 1.235294 | 0.161905 | 1.704e-02 | 0.131 | 8.713e-01 | 4.1e-07 |
| 0.90 | 1.222222 | 0.163636 | 7.224e-03 | 0.085 | 1.566e+00 | 4.2e-09 |
| 0.94 | 1.212766 | 0.164912 | 2.513e-03 | 0.050 | 3.781e+00 | 1.3e-06 |

f = 0.94 is the solver's validated ceiling: the w₁-continuation stalls at
f ≈ 0.947 regardless of mesh size (approaching the photon-tiring limit
f → 1 needs a matched asymptotic treatment, not more resolution).

### 2. `winds/runs_f/models/f{F}_Mdot{X}.npz` — physical profiles (`winds/generate_model_1d_f.py`)

Dimensionalized profiles for f ∈ {0.60, 0.90, 0.94} ×
Ṁ ∈ {2.5, 5, 10, 20} M☉/yr: CGS arrays `r, v, rho, tau` (Thomson depth)
plus all scalars (R, v_esc, L₀, τ★, τ_base, …).

### 3. `winds/sirocco_imports_f/outflow_f{F}_Mdot{X}.import.txt` — Sirocco inputs

The same models as 1-D spherical imports for the Sirocco radiative-transfer
code: columns (i, r [cm], v_r [cm/s], ρ [g/cm³], T_e [K]) on a 100-cell log
grid from r_in (τ_es = 100 surface, or R where the base is more transparent)
to 10× the τ_es = 1 photosphere, with one inner ghost cell at T_core.

## Reproducing everything from scratch

```bash
pip install -r requirements.txt

# Layer 1: dimensionless solutions (minutes; f=0.94 is the slow one)
cd winds
python solve_f_grid.py 0.60 0.75 0.85 0.90 0.94

# Layers 2+3: physical models + Sirocco imports
python generate_model_1d_f.py

# Paper figure (figures/FigPhotonTiredAndMBProfilesTiring.{png,pdf})
cd ../analysis
jupyter nbconvert --to notebook --execute --inplace FigPhotonTiredAndMBProfiles.ipynb
```

`solve_f_grid.py` with no arguments solves its full default grid, which
includes f = 0.95 — expect that one to fail with a RuntimeError explaining
the stall; this is the documented solver limit, not a bug.

Note on bit-level reproducibility: `scipy.integrate.solve_bvp` mesh
refinement can vary slightly across SciPy versions, so regenerated `.npz`
files may differ in the last digits. The physical results are insensitive to
this; `requirements.txt` records the exact versions used
(Python 3 + numpy 2.4.2, scipy 1.16.2, matplotlib 3.10.6).

## Figures

`analysis/FigPhotonTiredAndMBProfiles.ipynb` reads `winds/runs_f/models/` and produces the
stacked velocity/density figure (photon-tired winds at f = 0.94, 0.90, 0.60
with τ_es = 1 photosphere markers, plus analytic marginally-unbound
profiles). Pre-rendered copies are in `figures/`.
