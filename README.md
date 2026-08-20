# LRDs as Super-Eddington Fountain Flows — reproducibility package

This repository accompanies the paper (`paper.pdf`) and
contains everything needed to reproduce its figures: the wind-model solver and
generation scripts, the modifications to the Sirocco radiative-transfer
code, the input files for every simulation, and one Jupyter notebook
per figure. The simulation *outputs* the notebooks read are archived
separately on Zenodo (DOI below).

## What's in this repository

```
winds/             wind-model solver (included as a git submodule), generation scripts, and all wind-model data
runs/              Sirocco simulation configurations (input files and example job scripts, one folder per run)
sirocco_patches/   our six modifications to the public Sirocco code as numbered patch files
analysis/          figure notebooks
figures/           the rendered paper figures
ENVIRONMENT.md     step-by-step guide to building the modified Sirocco and re-running the simulations
paper_draft.pdf    the manuscript this package accompanies
```



## Three ways to use this repo


| You want to…                                                | Ran on                   | Time              | Follow                          |
| ----------------------------------------------------------- | ------------------------ | ----------------- | ------------------------------- |
| Reproduce every figure from the archived simulation outputs | a laptop                 | ~30 min           | **Part 1**                      |
| Regenerate the analytic wind models                         | a laptop                 | minutes           | **Part 2**                      |
| Re-run the radiative-transfer simulations from scratch      | a cluster with MPI + GSL | days of computing | **Part 3** and `ENVIRONMENT.md` |


---



## Part 1 — Reproduce the figures (laptop, ~30 minutes)

**Step 1. Download this repository.**

```bash
git clone --recurse-submodules <this-repo-url>
cd <repo-folder>
```

Without `--recurse-submodules` the `winds/tiring_solver/` submodule is empty. `That only matters for Part 2`.

**Step 2. Download and unpack the simulation outputs.**

Download `sirocco_outputs.tar.gz` (~84 MB) from the Zenodo archive (DOI to be
added here) into this folder, then unpack it:

```bash
mkdir sirocco_data
tar -xzf sirocco_outputs.tar.gz -C sirocco_data
```

This creates `sirocco_data/repro/` with one folder per simulation, plus
`README_DATA.md` (describing every file) and `MANIFEST.sha256`. To verify the
download:

```bash
cd sirocco_data && sha256sum -c MANIFEST.sha256 && cd ..
```

**Step 3. Tell the notebooks where the data is.**

```bash
export SIROCCO_REPRO_DATA=/full/path/to/sirocco_data
```

**Step 4. Check your Python packages.**

The notebooks need numpy, scipy, pandas, matplotlib and Jupyter. You likely have these already; if not, run:

```bash
pip install -r analysis/requirements.txt   # optional
```

**Step 5. Run the notebooks.**

These can be ran interactively or to produce all the figures from the command line, run:

```bash
cd analysis
for nb in Fig*.ipynb; do jupyter nbconvert --to notebook --execute --inplace "$nb"; done
```

Each notebook writes its figure into `figures/` as a PDF and as a PNG. 

---

## Part 2 — Wind models

The wind models are solutions of the Owocki et al. (2017) photon-tired wind
equations. The solver is the `winds/tiring_solver/` submodule (run`git submodule update --init` if it is empty); it has its own readme, and `validate.py` re-verifies every solution against the governing
equations.

```bash
# pip install -r winds/requirements.txt   # optional, probably unnecessary

# Dimensionless solutions (minutes; f=0.94 is the slow one)
cd winds
python solve_f_grid.py 0.60 0.75 0.85 0.90 0.94

# Physical wind profiles + Sirocco input files
python generate_model_1d_f.py

# The analytic-profiles paper figure
cd ../analysis
jupyter nbconvert --to notebook --execute --inplace FigPhotonTiredAndMBProfiles.ipynb
```

Note: `solve_f_grid.py` run with *no* arguments will also try f = 0.95 and prints `FAILED:` the solver fails at this limit. 

The marginally-unbound wind inputs are analytic and come from `winds/generate_model_1d_mu.py`
the same way. Both generator scripts print the Sirocco
`Central_object.{radius,temp,dilution_factor}` and `Wind.radmax` values
implied by each model (the diluted-interior rule
T_color = T_core·τ_in^(1/4), W = 1/τ_in with τ_in = min(100, τ_base)).

---



## Part 3 — Re-run the radiative-transfer simulations

The spectra and wind diagnostics come from **Sirocco** (Monte-Carlo radiative
transfer, [https://github.com/sirocco-rt/sirocco](https://github.com/sirocco-rt/sirocco)) using a specific version plus

 six modifications in `sirocco_patches/`. Each `runs/<name>/` folder contains all the necessary input and example job scripts. You will need to modify the job scripts for your cluster. 

`ENVIRONMENT.md` **walks through Sirocco-specific instructions command by command, plus the post-processing that produces tables for the notebooks to read.** Once your own runs are post-processed, point `SIROCCO_REPRO_DATA` at them and
Part 1 works unchanged.

---



# Reference

## Wind data(under `winds/`)

**1.** `winds/runs_f/f{F}.npz` **are dimensionless solutions** (from `solve_f_grid.py`): arrays of `x, w, q, p` (see Owocki+2017 for definitions) plus scalars.`ode_residual` is the stored normalized residual of each verified solution. The runs used in the paper are,


| f = m/m_max | Γ₀       | m        | w₁        | v∞/v_esc | q₀        | ODE residual |
| ----------- | -------- | -------- | --------- | -------- | --------- | ------------ |
| 0.60        | 1.333333 | 0.150000 | 1.045e-01 | 0.323    | 4.361e-01 | 3.1e-07      |
| 0.75        | 1.266667 | 0.157895 | 4.636e-02 | 0.215    | 5.715e-01 | 1.5e-07      |
| 0.85        | 1.235294 | 0.161905 | 1.704e-02 | 0.131    | 8.713e-01 | 4.1e-07      |
| 0.90        | 1.222222 | 0.163636 | 7.224e-03 | 0.085    | 1.566e+00 | 4.2e-09      |
| 0.94        | 1.212766 | 0.164912 | 2.513e-03 | 0.050    | 3.781e+00 | 1.3e-06      |


**2.** `winds/runs_f/models/f{F}_Mdot{X}.npz` **are physical profiles** (from `generate_model_1d_f.py`) in CGS. The arrays are `r, v, rho, tau` (Thomson depth) plus all scalars (R, v_esc, L₀, τ★, τ_base, …) for f ∈ {0.60, 0.90, 0.94} × Ṁ ∈ {2.5, 5, 10, 15} M☉/yr.

**3.** `winds/sirocco_imports_f/` **are the Sirocco input tables**. These are columns (i, r [cm], v_r [cm/s], ρ [g/cm³], T_e [K]) on a log grid — 100 cells, or 300 for the production f = 0.94, Ṁ = 15 import (`--ncells 300`; see
`ENVIRONMENT.md`, "Radial resolution") — from r_in (τ_es = 100 surface, or R where the base is more transparent) to 10× the τ_es = 1 photosphere, with one inner ghost cell at T_core. Most copies are in `runs/<name>/outflow.import.txt`; the 300-cell production import (`outflow_f0.94_Mdot15.0_n300.import.txt`) is provided here instead. 

## Simulations

The simulations are divided into two wind families with four mass-loss rates each (Ṁ = 2.5, 5, 10, 15 M☉/yr):


| physical family                               | `runs/` dirs | notebook keys | data dir (under `$SIROCCO_REPRO_DATA/repro/`) |
| --------------------------------------------- | ------------ | ------------- | --------------------------------------------- |
| marginally unbound (v = v_esc, "dilute_redd") | `mu_m*`      | `mu_m*`       | same names                                    |
| photon-tired, f = 0.94 (Γ₀ = 1.2128)          | `ptf_m*`     | `pt_m*`       | same names                                    |


`*_trunc_speconly` dirs are spectrum-only reruns on truncated winds (the  
production source for the two m15 emergent spectra) and `halpha_hires/`  
subdirs are narrow-band Hα reruns — both documented in `ENVIRONMENT.md`.  
`runs/mu_m15/` also carries `outflow.restart.slurm`, a job script for a spectrum restart that completed its 20 spectrum cycles after a wall-time timeout.

