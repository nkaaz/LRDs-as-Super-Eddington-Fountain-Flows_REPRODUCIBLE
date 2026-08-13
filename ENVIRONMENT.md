# Radiative-transfer environment (Sirocco)

The synthetic spectra and per-cell wind diagnostics in this paper were produced
with **Sirocco** (Monte-Carlo radiative transfer with macro-atom NLTE;
https://github.com/sirocco-rt/sirocco) — **not a stock release**, but upstream
commit `5aef5f17` plus six small local patches shipped verbatim in
`sirocco_patches/`. Stock Sirocco will **reject the `.pf` files in `runs/`**
(they use a parameter it does not have) and, worse, would **silently
reproduce a photon-loss bug** in the optically-thick runs (patch 0004).
Build the patched code.

## Building the exact code

```bash
git clone https://github.com/sirocco-rt/sirocco.git
cd sirocco
git checkout 5aef5f17
git am /path/to/this-repo/sirocco_patches/*.patch
./configure          # needs MPI + GSL
cd source
make sirocco windsave2table modify_wind
```

Note: `make` runs a python-based indent step after linking; if the `python`
command is missing it errors *after* the binaries are already built — that
failure can be ignored.

Cluster modules used for the paper runs (Princeton Stellar, Intel nodes):
`gcc-toolset/13 openmpi/gcc/4.1.6 gsl/2.6`. Any MPI + GSL stack should work;
the `.slurm` scripts in `runs/` hardcode these modules and the author's paths,
so treat them as a record of the invocation (ranks, photons, cycles,
wall-times), not as portable submission scripts.

## The six patches

| # | Patch | Why it is required |
|---|---|---|
| 0001 | `Central_object.dilution_factor` | New `.pf` parameter: the central source emits a blackbody *spectrum* at `Central_object.temp` (the color temperature) with luminosity multiplied by W, so L = W·4πR²σT⁴. Every run here injects a diluted interior radiation field (W = 1/τ_es(r_in)); without this parameter no `.pf` in `runs/` parses. |
| 0002 | levden tables in `windsave2table` | Writes `rootname.<Elem>_<istate>.levden.txt` (NLTE level populations). Four notebooks (`FigGamma1`, `FigHaHbNetEmission`, `FigOpticalDepth`, `FigTemp_n2_grid`) read `outflow.H_1.levden.txt` — the n=2 populations behind the Balmer-break analysis. |
| 0003 | `NBINS_IN_CELL_SPEC` 1000→3000 | Per-cell J_ν resolution of `windsave2table -xall` (the `outflow.xspec.all.txt` files the notebooks read have 3001 rows). |
| 0004 | `MAXSCAT` 2000→200000 | **Critical.** These winds have τ_es(base) up to ~100; photons legitimately scatter tens of thousands of times. With the stock cap of 2000, ~80% of in-flight photons are silently destroyed and the photon-tired spectra come out ~20× too faint. |
| 0005 | `xdata/h20_hetop_standard80_z0.1.dat` | The 10%-metallicity atomic dataset every `.pf` references (`Atomic_data data/h20_hetop_standard80_z0.1.dat`). Recipe (also in the commit message): copy `h20_hetop_standard80` masterfile chain and scale every `Element` line with z ≥ 3 by −1.0 dex (abundances are log-number on the H=12 scale). |
| 0006 | `modify_wind -rcut` | Builds the truncated wind-saves used by the spectrum-only reruns (below). |

## Running a model

Each `runs/<name>/` contains the complete input set: `outflow.pf`,
`outflow.import.txt` (the 1-D spherical wind model; see `winds/` for how these
are generated), and the `.slurm` record. Run dirs need `data` and `zdata`
symlinks to the Sirocco `xdata`/`zdata` directories. Production settings:
2×10⁶ photons/cycle, `matrix_pow` ionization, `wide` photon sampling,
macro-atom H/He + simple-atom metals, 384 MPI ranks (4 nodes), 14 ionization
+ 20 (mu) or 40 (ptf) spectrum cycles. Use `-p` (logarithmic photon ramp-up)
for fresh runs only — never for `System_type=previous` restarts.

## Post-processing products the notebooks read

From a completed run directory:

```bash
windsave2table -xall outflow
```

produces `outflow.master.txt`, `outflow.heat.txt`,
`outflow.H_1.levden.txt` (patch 0002), and `outflow.xspec.all.txt`
(per-cell J_ν, patch 0003). Emergent spectra are `outflow.log_spec`
(log-uniform frequency grid — use this for optical work; the linear
`outflow.spec` under-resolves lines) with the 45° observer in column
`A45P0.50`.

## Truncated spectrum-only reruns (`*_trunc_speconly`)

In the highest-Ṁ runs the wind *beyond the hydrogen recombination front*
never converges (the gas there decouples from the radiation field and t_e
collapses); it is causally irrelevant to the interior solution but
contaminates the emergent blue/UV. The paper therefore takes those spectra
from spectrum-only reruns on a truncated wind:

```bash
# in the parent run dir, after completion:
cp outflow.wind_save full.wind_save
modify_wind -rcut <RCUT_CM> -out_root cut full     # floors ρ, n_e beyond RCUT
# then run spec.pf: System_type=previous, Wind.old_windfile=cut,
# Ionization_cycles=0 (see runs/ptf_m*_trunc_speconly/)
```

RCUT is 1.35× the recombination front radius of the converged parent
(hand-checked against the per-cell convergence pattern):

| run | recomb front | RCUT | cells kept |
|---|---|---|---|
| ptf_m10 | 3.05e17 cm (cell 47) | 4.111e17 cm | 0–49 of 99 |
| ptf_m20 | 4.4e17 cm (cell 36) | 5.958e17 cm | 0–40 of 99 |

Validation: the truncation changes the Balmer-red and NIR spectrum by ≤3%
(Hα ≤4%); everything blueward of the Balmer edge depends on the unconverged
outer wind in *both* treatments (parent and truncated disagree by 1–4 dex
there) and is not a robust prediction of these models — the gallery figure
renders that region grey.

## Hα line profiles (`halpha_hires/`)

`Fig_LRD_halpha` prefers a narrow-band restart when present:
`runs/<name>/halpha_hires/outflow.pf` = the run's `.pf` with
`System_type=previous` (reading a copy of the parent `outflow.wind_save`),
`Ionization_cycles=0`, and the spectral window narrowed to 6500–6700 Å with
10⁴ bins — all 2×10⁶ photons/cycle are then generated in-band, giving
~1 km/s sampling of the line.

## Data availability

The Sirocco *outputs* (wind-saves, spectra, per-cell tables) are not in this
git repository; they will be archived on Zenodo (DOI to be added here).
The notebooks locate them via two environment variables (see README):
`SIROCCO_REPRO_DATA` (data root; the run directories live under
`<root>/repro/<name>` with the same names as `runs/<name>`) and
`SIROCCO_REPRO_RUNSET` (`repro` = the production runs documented here,
default; `legacy` = the paper-freeze directory layout on the authors'
cluster, retained for provenance).
