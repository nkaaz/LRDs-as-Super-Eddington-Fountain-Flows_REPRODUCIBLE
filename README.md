# LRDs as Super-Eddington Fountain Flows — wind models (reproducibility package)

Everything needed to reproduce the paper's figures: the wind-model solver
(pinned git submodule) and generation scripts, the patches that define the
exact Sirocco radiative-transfer code, all Sirocco run configurations, and one
notebook per figure. The Sirocco *outputs* the notebooks read are archived on
Zenodo (DOI below).

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

    winds/             solver (submodule), generation scripts, and all wind data
    runs/              Sirocco run configurations (.pf + .slurm record per run; truncated and halpha_hires reruns as subdirectories)
    sirocco_patches/   the six patches that turn stock Sirocco (5aef5f17) into the code used here
    analysis/          figure notebooks (one per paper figure; see manifest below)
    figures/           rendered paper figures (PDF)
    ENVIRONMENT.md     how to build the patched Sirocco, run the models, and post-process them

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
Ṁ ∈ {2.5, 5, 10, 15, 20} M☉/yr: CGS arrays `r, v, rho, tau` (Thomson depth)
plus all scalars (R, v_esc, L₀, τ★, τ_base, …).

### 3. `winds/sirocco_imports_f/` — Sirocco inputs (regenerated on demand)

The same models as 1-D spherical imports for the Sirocco radiative-transfer
code: columns (i, r [cm], v_r [cm/s], ρ [g/cm³], T_e [K]) on a log grid
(100 cells; 300 for the production f = 0.94, Ṁ = 15 import — set `N_CELLS`
in `generate_model_1d_f.py`; see `ENVIRONMENT.md`, "Radial resolution") from
r_in (τ_es = 100 surface, or R where the base is more transparent) to 10× the
τ_es = 1 photosphere, with one inner ghost cell at T_core. The analogous
marginally-unbound imports come from `winds/generate_model_1d_mu.py`. The
copies actually run live in `runs/<name>/outflow.import.txt`; only the
300-cell production import (`outflow_f0.94_Mdot15.0_n300.import.txt`, from
`python generate_model_1d_f.py --ncells 300`) is also kept here. Both
generator scripts print the Sirocco `Central_object.{radius,temp,
dilution_factor}` and `Wind.radmax` values implied by each import (the
diluted-interior rule T_color = T_core·τ_in^(1/4), W = 1/τ_in with
τ_in = min(100, τ_base)).

## Reproducing everything from scratch

```bash
pip install -r winds/requirements.txt   # winds scripts. NB: analysis/requirements.txt
                                        # pins DIFFERENT versions — use a separate venv
                                        # (or any recent stack) for the notebooks.

# Layer 1: dimensionless solutions (minutes; f=0.94 is the slow one)
cd winds
python solve_f_grid.py 0.60 0.75 0.85 0.90 0.94

# Layers 2+3: physical models + Sirocco imports
python generate_model_1d_f.py

# Paper figure (figures/FigPhotonTiredAndMBProfiles.{png,pdf})
cd ../analysis
jupyter nbconvert --to notebook --execute --inplace FigPhotonTiredAndMBProfiles.ipynb
```

`solve_f_grid.py` with no arguments solves its full default grid, which
includes f = 0.95 — expect that one to print `FAILED:` (the stall) and
continue; this is the documented solver limit, not a bug.

Note on bit-level reproducibility: `scipy.integrate.solve_bvp` mesh
refinement can vary slightly across SciPy versions, so regenerated `.npz`
files may differ in the last digits. The physical results are insensitive to
this; `requirements.txt` records the exact versions used
(Python 3 + numpy 2.4.2, scipy 1.16.2, matplotlib 3.10.6).

## Figures

`analysis/FigPhotonTiredAndMBProfiles.ipynb` reads `winds/runs_f/models/` and produces the
stacked velocity/density figure (photon-tired winds at f = 0.94, 0.90, 0.60
with τ_es = 1 photosphere markers, plus analytic marginally-unbound
profiles) as `figures/FigPhotonTiredAndMBProfiles.{png,pdf}`.
Pre-rendered **PDF** copies of all figures are in `figures/` (PNGs are not
tracked).

---

# Radiative-transfer layer (Sirocco): `runs/`, `analysis/`, `figures/`

Everything below this line concerns the Monte-Carlo radiative-transfer models
that produce the paper's spectra and wind diagnostics. The code is **not stock
Sirocco** — see `ENVIRONMENT.md` for building the patched version
(`sirocco_patches/`), running the models, and generating the post-processing
tables the notebooks read.

## Wind families and naming

Two wind families, four mass-loss rates each (Ṁ = 2.5, 5, 10, 15 M☉/yr):

| physical family | `runs/` dirs | notebook keys | data dir (under `$SIROCCO_REPRO_DATA/repro/`) |
|---|---|---|---|
| marginally unbound (v = v_esc, "dilute_redd") | `mu_m*` | `mu_m*`, `MU_m*`, or bare `m*` | same names |
| photon-tired, f = 0.94 (Γ₀ = 1.2128) | `ptf_m*` | `pt_m*` / `PT_m*` | same names |

The photon-tired Ṁ = 15 run is `ptf_m15_n300`: it uses a 300-cell radial grid
because at the standard 100 cells its hydrogen-recombination front does not
converge (a grid-resolution limit cycle; see `ENVIRONMENT.md`, "Radial
resolution"). All other runs use 100 cells, verified steady.

`*_trunc_speconly` dirs are spectrum-only reruns on truncated winds (the
production source for all emergent spectra) and `halpha_hires/` subdirs are
narrow-band Hα reruns — both documented in `ENVIRONMENT.md`. `runs/mu_m15/`
also carries `outflow.restart.slurm`, the `-r` spectrum restart that completed
its 20 spectrum cycles after a walltime timeout.

## Data location: two environment variables

The nine Sirocco notebooks read post-processing tables and spectra that are
**not in this repo** (they will be archived on Zenodo; DOI to be added):

```bash
export SIROCCO_REPRO_DATA=/path/to/data      # dir containing repro/<run-name>/
```

Each Sirocco notebook's first code cell has a `# ==== config ====` block with
the run-name map (`RUN_DIRS`, or `PAIRS` in the truncation-comparison appendix
notebook) — that map is the authoritative record of which run feeds which
panel. `FigPhotonTiredAndMBProfiles` has no such block: it reads only
`winds/runs_f/models/` inside the repo.

## Figure manifest

| notebook | figure file(s) | families | data files read per run |
|---|---|---|---|
| `Fig_LRD_gallery.ipynb` | `Fig_LRD_gallery` | mu + ptf (8 panels) | `outflow.log_spec` / `spec.log_spec` |
| `Fig_LRD_halpha.ipynb` | `Fig_LRD_halpha` | mu + ptf (8 panels) | `log_spec` + `halpha_hires/outflow.log_spec` |
| `FigHeatCool.ipynb` | `FigHeatCool` | mu + ptf | `outflow.heat.txt`, `outflow.xspec.all.txt` |
| `FigHeatCool_TeXi.ipynb` | `FigHeatCool_TeXi` | mu + ptf | `outflow.heat.txt`, `outflow.master.txt` |
| `FigHaHbNetEmission.ipynb` | `FigHaHbNetEmission` | mu + ptf | `outflow.master.txt`, `outflow.H_1.levden.txt` |
| `FigTemp_n2_grid.ipynb` | `FigTemp_n2_grid` | mu + ptf | `master`, `H_1.levden`, `xspec.all` |
| `FigGamma1.ipynb` | `FigGamma1` | mu only (by design) | `master`, `H_1.levden`, `xspec.all` |
| `FigOpticalDepth.ipynb` | `FigOpticalDepth` | mu only (by design) | `master`, `H_1.levden`, `xspec.all` |
| `FigPhotonTiredAndMBProfiles.ipynb` | `FigPhotonTiredAndMBProfiles` | (analytic; `winds/` only — runs without Sirocco data) | `winds/runs_f/models/*.npz` |
| `FigAppTruncatedSpectraComparison.ipynb` | `FigAppTruncatedSpectraComparison` | mu_m15 + ptf_m15, truncated vs full domain | `log_spec` (+ `halpha_hires`) of both treatments |

**Truncation policy:** `Fig_LRD_gallery` reads the truncated
(converged-cells-only) reruns for mu_m15, ptf_m10 and ptf_m15;
`Fig_LRD_halpha` reads the truncated products for the two m15 runs and the
full-domain run for ptf_m10 (whose Hα is truncation-insensitive at the ≤4%
level). Flux blueward
of the Balmer edge depends on the unconverged outer wind (full-domain and
truncated treatments disagree there by orders of magnitude) and is not a
robust prediction of these models. `FigAppTruncatedSpectraComparison.ipynb`
quantifies this: redward of the Balmer edge the two treatments agree to
better than a factor ~2; blueward they diverge.
