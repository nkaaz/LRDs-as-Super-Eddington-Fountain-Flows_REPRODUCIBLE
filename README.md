# LRDs as Super-Eddington Fountain Flows — reproducibility package

This repository accompanies the paper (included as `paper_draft.pdf`) and
contains everything needed to reproduce its figures: the wind-model solver and
generation scripts, the exact modifications to the Sirocco radiative-transfer
code, the complete input files for every simulation, and one Jupyter notebook
per figure. The simulation *outputs* the notebooks read are archived
separately on Zenodo (DOI below) because they are too large for git.

## What's in this repository

    winds/             wind-model solver (git submodule), generation scripts, and all wind-model data
    runs/              Sirocco simulation configurations (input files + job-script record, one folder per run)
    sirocco_patches/   our six modifications to the public Sirocco code, as numbered patch files
    analysis/          figure notebooks (one per paper figure; see the manifest at the bottom)
    figures/           the rendered paper figures (PDF)
    ENVIRONMENT.md     step-by-step guide to building the modified Sirocco and re-running the simulations
    paper_draft.pdf    the manuscript this package accompanies

## Three ways to use this package

| You want to… | You need | Time | Follow |
|---|---|---|---|
| Reproduce every figure from the archived simulation outputs | a laptop with Python | ~30 min | **Part 1** (below) |
| Regenerate the analytic wind models themselves | a laptop with Python | minutes | **Part 2** |
| Re-run the radiative-transfer simulations from scratch | a Linux cluster with MPI + GSL | days of computing | **Part 3** and `ENVIRONMENT.md` |

Most readers want Part 1. It does **not** require installing or building any
simulation software.

---

## Part 1 — Reproduce the figures (laptop, ~30 minutes)

**Step 1. Download this repository.**

```bash
git clone --recurse-submodules <this-repo-url>
cd <repo-folder>
```

If you cloned without `--recurse-submodules` (or downloaded a ZIP), the folder
`winds/tiring_solver/` will be empty. That's fine for this part — it is only
needed in Part 2. To fetch it later: `git submodule update --init`.

**Step 2. Download and unpack the simulation outputs.**

From the Zenodo archive (DOI to be added here), download
`sirocco_outputs.tar.gz` (~84 MB). Move the downloaded file into your
current folder (browsers usually save to `~/Downloads`), then unpack it:

```bash
mkdir sirocco_data
tar -xzf sirocco_outputs.tar.gz -C sirocco_data
```

This creates `sirocco_data/repro/` with one folder per simulation, plus
`README_DATA.md` (describing every file) and `MANIFEST.sha256`. To verify the
download is intact:

```bash
cd sirocco_data && sha256sum -c MANIFEST.sha256 && cd ..
```

**Step 3. Tell the notebooks where the data is.**

```bash
export SIROCCO_REPRO_DATA=/full/path/to/sirocco_data
```

This sets an *environment variable* — a named value that programs started
from this terminal can read. It only lasts for the current terminal session,
so run it in the **same terminal** you will start Jupyter from (or add the
line to your `~/.bashrc`). If you forget this step, the notebooks fall back
to the authors' cluster path and fail with an error like
`FileNotFoundError: /scratch/gpfs/nk6352/...` — that error means "set
`SIROCCO_REPRO_DATA`".

**Step 4. Install the Python packages.**

```bash
pip install -r analysis/requirements.txt
```

These are the exact versions used for the paper; any reasonably recent
numpy/scipy/matplotlib/jupyter stack also works.

**Step 5. Run the notebooks.**

Either open them interactively (`jupyter lab analysis/`) and run all cells in
each, or render everything from the command line:

```bash
cd analysis
for nb in Fig*.ipynb; do jupyter nbconvert --to notebook --execute --inplace "$nb"; done
```

Each notebook writes its figure into `figures/` as a PDF (and a PNG, which is
not tracked by git). The PDFs already committed in `figures/` were produced by
exactly this procedure, so you can compare your output against them directly.
The first code cell of each notebook has a `# ==== config ====` block naming
the simulation folder(s) it reads — that block is the authoritative record of
which run feeds which panel. Which notebook makes which paper figure is
listed in the [figure manifest](#figure-manifest) below.

---

## Part 2 — Regenerate the wind models (laptop, minutes)

The wind models are solutions of the Owocki et al. (2017) photon-tired wind
equations; the solver lives in the `winds/tiring_solver/` submodule (fetch it
with `git submodule update --init` if the folder is empty; its own README
documents the method and range of validity, and `validate.py` re-verifies
every solution against the governing equations).

```bash
pip install -r winds/requirements.txt   # NB: pins DIFFERENT versions than
                                        # analysis/requirements.txt — use a
                                        # separate virtual environment (or any
                                        # recent stack) for the notebooks.

# Dimensionless solutions (minutes; f=0.94 is the slow one)
cd winds
python solve_f_grid.py 0.60 0.75 0.85 0.90 0.94

# Physical wind profiles + Sirocco input files
python generate_model_1d_f.py

# The analytic-profiles paper figure
cd ../analysis
jupyter nbconvert --to notebook --execute --inplace FigPhotonTiredAndMBProfiles.ipynb
```

Two expected oddities: (1) `solve_f_grid.py` run with *no* arguments also
tries f = 0.95 and prints `FAILED:` for it — that is the documented solver
limit (the continuation stalls at f ≈ 0.947 regardless of resolution), not a
bug. (2) Regenerated `.npz` files can differ from the committed ones in the
last decimal digits, because `scipy.integrate.solve_bvp` mesh refinement
varies slightly across SciPy versions; the physical results are insensitive
to this.

The marginally-unbound wind inputs come from `winds/generate_model_1d_mu.py`
the same way. Both generator scripts print the Sirocco
`Central_object.{radius,temp,dilution_factor}` and `Wind.radmax` values
implied by each model (the diluted-interior rule
T_color = T_core·τ_in^(1/4), W = 1/τ_in with τ_in = min(100, τ_base)).

---

## Part 3 — Re-run the radiative-transfer simulations (cluster, days)

The spectra and wind diagnostics come from **Sirocco** (Monte-Carlo radiative
transfer, https://github.com/sirocco-rt/sirocco) — not a stock release, but a
pinned version plus the six modifications in `sirocco_patches/`. Each
`runs/<name>/` folder contains a complete input set, and the included job
scripts record the exact invocations (ranks, photons, cycles, wall-times).
Expect each production model to need several hundred CPU cores for hours to
days.

**`ENVIRONMENT.md` walks through all of it command by command — including how
to apply the patches, assuming no prior experience with git — plus the
post-processing that turns raw outputs into the tables the notebooks read.**
Once your own runs are post-processed, point `SIROCCO_REPRO_DATA` at them and
Part 1 works unchanged.

---

# Reference

## The model in brief

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

## Wind-model data (three layers, all under `winds/`)

**1. `winds/runs_f/f{F}.npz` — dimensionless solutions** (from
`solve_f_grid.py`): arrays `x, w, q, p` (4000-point uniform x-grid) plus
scalars. Boundary condition w_sg = 0 (idealized zero-sound-speed base);
`ode_residual` is the stored normalized residual of each verified solution.

| f = m/m_max | Γ₀ | m | w₁ | v∞/v_esc | q₀ | ODE residual |
|---|---|---|---|---|---|---|
| 0.60 | 1.333333 | 0.150000 | 1.045e-01 | 0.323 | 4.361e-01 | 3.1e-07 |
| 0.75 | 1.266667 | 0.157895 | 4.636e-02 | 0.215 | 5.715e-01 | 1.5e-07 |
| 0.85 | 1.235294 | 0.161905 | 1.704e-02 | 0.131 | 8.713e-01 | 4.1e-07 |
| 0.90 | 1.222222 | 0.163636 | 7.224e-03 | 0.085 | 1.566e+00 | 4.2e-09 |
| 0.94 | 1.212766 | 0.164912 | 2.513e-03 | 0.050 | 3.781e+00 | 1.3e-06 |

f = 0.94 is the solver's validated ceiling (see Part 2).

**2. `winds/runs_f/models/f{F}_Mdot{X}.npz` — physical profiles** (from
`generate_model_1d_f.py`): CGS arrays `r, v, rho, tau` (Thomson depth) plus
all scalars (R, v_esc, L₀, τ★, τ_base, …) for f ∈ {0.60, 0.90, 0.94} ×
Ṁ ∈ {2.5, 5, 10, 15} M☉/yr.

**3. `winds/sirocco_imports_f/` — Sirocco input tables** (regenerated on
demand): the same models as 1-D spherical imports, columns
(i, r [cm], v_r [cm/s], ρ [g/cm³], T_e [K]) on a log grid — 100 cells, or
300 for the production f = 0.94, Ṁ = 15 import (`--ncells 300`; see
`ENVIRONMENT.md`, "Radial resolution") — from r_in (τ_es = 100 surface, or R
where the base is more transparent) to 10× the τ_es = 1 photosphere, with one
inner ghost cell at T_core. The copies actually run live in
`runs/<name>/outflow.import.txt`; only the 300-cell production import
(`outflow_f0.94_Mdot15.0_n300.import.txt`) is also kept here.

## Simulation families and naming

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
production source for the two m15 emergent spectra) and `halpha_hires/`
subdirs are narrow-band Hα reruns — both documented in `ENVIRONMENT.md`.
`runs/mu_m15/` also carries `outflow.restart.slurm`, the spectrum restart
that completed its 20 spectrum cycles after a wall-time timeout.

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
| `FigPhotonTiredAndMBProfiles.ipynb` | `FigPhotonTiredAndMBProfiles` | (analytic; `winds/` only — runs without the data archive) | `winds/runs_f/models/*.npz` |
| `FigAppTruncatedSpectraComparison.ipynb` | `FigAppTruncatedSpectraComparison` | mu_m15 + ptf_m15, truncated vs full domain | `log_spec` (+ `halpha_hires`) of both treatments |

(`spec.log_spec` is the same spectrum product as `outflow.log_spec`; the
`*_trunc_speconly` reruns simply name their output files `spec.*` instead of
`outflow.*`.)

**Truncation policy:** the two spectra figures read the truncated
(converged-cells-only) reruns for the two m15 runs; all other runs, including
ptf_m10, are read full-domain. Flux blueward of the Balmer edge depends on
the unconverged outer wind (full-domain and truncated treatments disagree
there by orders of magnitude) and is not a robust prediction of these models.
`FigAppTruncatedSpectraComparison.ipynb` quantifies this: redward of the
Balmer edge the two treatments agree to better than a factor ~2; blueward
they diverge.
