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
| 0005 | `xdata/h20_hetop_standard80_z0.1.dat` | The 10%-metallicity atomic dataset every `.pf` references (`Atomic_data data/h20_hetop_standard80_z0.1.dat`). **The dataset itself is shipped in this patch** — no manual step; the commit message documents how it was made (every `Element` line with z ≥ 3 scaled −1.0 dex; abundances are log-number on the H=12 scale). |
| 0006 | `modify_wind -rcut` | Builds the truncated wind-saves used by the spectrum-only reruns (below). |

## Running a model

Each `runs/<name>/` contains the complete input set: `outflow.pf`,
`outflow.import.txt` (the 1-D spherical wind model; see `winds/` for how these
are generated), and the `.slurm` record. Run dirs need `data` and `zdata`
symlinks to the Sirocco `xdata`/`zdata` directories. Production settings:
2×10⁶ photons/cycle, `matrix_pow` ionization, photon sampling `user_bands`
(5 bands, 0.05–100 eV, boundaries 1.5/3.4/7/13.6 eV) for the mu family and
`wide` for the ptf family, macro-atom H/He + simple-atom metals, 14 ionization + 20 (mu) or 40 (ptf)
spectrum cycles; 384 MPI ranks at 4 GB each (4 nodes) for 100-cell runs, or
192 ranks at 8 GB for the 300-cell run (matrix-mode memory scales with cell
count: ~1.8 GB/rank at 100 cells, ~5.4 at 300). Use `-p` (logarithmic photon ramp-up)
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
# in the <name>_trunc_speconly dir:
cp ../<name>/outflow.wind_save full.wind_save
modify_wind -rcut <RCUT_CM> -out_root cut full     # floors ρ, n_e beyond RCUT -> cut.wind_save
# then: sbatch spec.slurm  (spec.pf: System_type=previous, Wind.old_windfile=cut,
# Ionization_cycles=0, Spectrum_cycles=400 — the truncated winds are photon-starved
# in the optical, so they get 400 spectrum cycles vs the parents' 20/40)
```

For a trunc dir's `halpha_hires/` rerun, copy `cut.wind_save` into that
subdirectory first; for a parent's `halpha_hires/`, copy the parent's
`outflow.wind_save` instead (its `.pf` reads `Wind.old_windfile outflow`).

RCUT is set at the outer edge of the contiguously converged core:
run `windsave2table` on the parent, read the `converge` column of
`outflow.converge.txt` (0 = converged), and place the cut at the first cell
of the sustained unconverged tail (≥5 consecutive failed cells); isolated
interior failures (runs of <5 cells) are tolerated. For ptf_m10 this boundary
coincides with 1.35× the hydrogen recombination front:

| run | RCUT | criterion | cells kept |
|---|---|---|---|
| mu_m15 | 8.509e18 cm | converged-core end (sustained tail from cell 76) | 0–75 of 99 |
| ptf_m15_n300 | 5.245e17 cm | converged-core end (sustained tail from cell 137); keeps the whole recombination front (r = 3.59e17) | 0–136 of 299 |
| ptf_m10 | 4.111e17 cm | 1.35× the recombination front (cell 47, r = 3.05e17) | 0–49 of 99 |

Truncation sensitivity (quantified in `FigAppTruncatedSpectraComparison`):
for mu_m15, whose cut sits essentially at its electron-scattering photosphere
(0.83 r_ph), the truncated and full-domain spectra agree closely redward of
the Balmer edge; for the photon-tired runs, whose cuts sit inside their
photospheres (0.26–0.37 r_ph), the removed unconverged tail absorbs up to
~half the optical, so the truncated (converged-cells-only) spectrum is the
production statement. Everything blueward of the Balmer edge depends on the
unconverged outer wind in *both* treatments (they disagree by 1–4 dex there)
and is not a robust prediction of these models.

## Radial resolution

The photon-tired Ṁ = 15 model is run on a **300-cell** radial grid
(`ptf_m15_n300`; import generated with `N_CELLS = 300` in
`winds/generate_model_1d_f.py`). At the standard 100 cells this run's
hydrogen-recombination front never converges: it executes a bounded limit
cycle (front radius sloshing between 3.0 and 3.9e17 cm, converged-cell count
oscillating, emergent Hα varying by ×4 between snapshots). The cause is
numerical — the front cells are internally optically thick to the Balmer
continuum that mediates the front's feedback (τ_BaC ≈ 34 per cell at 100
cells), and neither more photons (the front cells receive 4–14 × 10⁶
photons/cycle; estimator noise 0.03–0.05%) nor more ionization cycles changes
it. At 300 cells (τ_BaC ≈ 11 per cell) the front converges monotonically and
sits at the mean of the 100-cell oscillation. All other runs (both families)
are steady at 100 cells; their spectra are unaffected by this choice.

## Restarts used in production

`runs/mu_m15/outflow.restart.slurm` records the spectrum-only `-r` restart
that completed mu_m15's spectrum cycles 11–20 after a walltime timeout (in
`-r` mode `Spectrum_cycles` is the cumulative total). Truncated reruns and
`halpha_hires` runs use `System_type=previous` (cycle counts are NEW cycles
in that mode) and never the `-p` photon ramp-up.

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
The notebooks locate them via one environment variable (see README):
`SIROCCO_REPRO_DATA` (data root; the run directories live under
`<root>/repro/<name>` with the same names as `runs/<name>`).
