# Radiative-transfer environment (Sirocco)

The synthetic spectra and per-cell wind diagnostics in this paper were produced
with **Sirocco** (Monte-Carlo radiative transfer with macro-atom NLTE;
https://github.com/sirocco-rt/sirocco) — **not a stock release**, but upstream
commit `5aef5f17` plus six small local patches shipped verbatim in
`sirocco_patches/`. Stock Sirocco will **reject the `.pf` files in `runs/`**
(they use a parameter it does not have) and, worse, would **silently
reproduce a photon-loss bug** in the optically-thick runs (patch 0004).
Build the patched code.

## Building the exact code, step by step

The code used for the paper is the public Sirocco source at one specific
version, plus our six modifications. The modifications ship in this
repository's `sirocco_patches/` folder as numbered **patch files** — plain
text files that record changes to the code in git's standard format. You do
not need to know git to use them: the commands below download the public
code, rewind it to the exact version our patches expect, and apply the
patches automatically. This whole sequence (through the dry run in step 5)
has been verified end-to-end on a fresh download.

Prerequisites: a Linux machine with git, an MPI C compiler (e.g. OpenMPI's
`mpicc`), and the GSL library. On the Princeton Stellar cluster, where the
paper runs were done, these come from
`module load gcc-toolset/13 openmpi/gcc/4.1.6 gsl/2.6`; any MPI + GSL stack
should work.

**Step 1. Download the public Sirocco source and rewind it to our base
version.**

```bash
git clone https://github.com/sirocco-rt/sirocco.git
cd sirocco
git checkout 5aef5f17
```

The last command switches the source tree to the exact snapshot our patches
apply to. Git prints a notice about a "detached HEAD" — that is normal and
harmless.

**Step 2. Apply the six patches.**

```bash
git am /path/to/this-repo/sirocco_patches/*.patch
```

(Replace `/path/to/this-repo` with wherever you put this repository. The `*`
expands to the six files in numeric order, which is the required order.)
You should see six lines starting with `Applying:`, one per patch — plus a
few hundred "trailing whitespace" warnings, which come from an atomic-data
table in patch 0005 and are harmless. Two things that can go wrong:

- If git asks you to say who you are ("Please tell me who you are"), run
  `git config --global user.name "Your Name"` and
  `git config --global user.email you@example.com`, then repeat the
  `git am` command. (Applying patches creates commits, and commits need an
  author name.)
- If a patch fails to apply — which should only happen if step 1 was skipped
  or altered — run `git am --abort` to reset, and start over from
  `git checkout 5aef5f17`.

The source tree is now identical, file for file, to the code the paper used.

**Step 3. Build.**

```bash
./configure
cd source
make -k sirocco windsave2table modify_wind
```

`make` ends with an error (`/usr/bin/env: 'python': No such file or
directory` or similar) on systems with no `python` command: each program's
build ends with a cosmetic code-formatting step that needs it. The error is
harmless — it happens *after* each program is already built — and the `-k`
flag tells `make` to keep going to the remaining programs anyway (without it,
only the first one gets built). So the overall nonzero exit is expected.
Success means all three programs exist:

```bash
ls ../bin/sirocco ../bin/windsave2table ../bin/modify_wind
```

**Step 4. Check the patches took effect.**

```bash
grep MAXSCAT sirocco.h                      # should say 200000, not 2000
ls ../xdata/h20_hetop_standard80_z0.1.dat   # the paper's atomic dataset
```

**Step 5 (optional). Dry-run a production input file.**

From some scratch directory:

```bash
mkdir smoketest && cd smoketest
cp /path/to/this-repo/runs/mu_m2.5/outflow.pf .
cp /path/to/this-repo/runs/mu_m2.5/outflow.import.txt .
ln -s /path/to/sirocco/xdata data
ln -s /path/to/sirocco/zdata zdata
/path/to/sirocco/bin/sirocco -i outflow
```

`-i` reads the inputs and sets up the wind without running any photons
(seconds). A clean exit ending in `Error summary: dry run.` with only a
handful of benign atomic-data notes means the build works; stock (unpatched)
Sirocco would instead reject the input file at the
`Central_object.dilution_factor` line.

Note on the `.slurm` scripts in `runs/`: they hardcode the Stellar modules
and the author's paths, so treat them as a record of the invocation (ranks,
photons, cycles, wall-times), not as portable submission scripts.

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

Each `runs/<name>/` contains the complete input set: `outflow.pf` (the
parameter file), `outflow.import.txt` (the 1-D spherical wind model; see
`winds/` for how these are generated), and the `.slurm` record. Before
running, the run directory needs two symbolic links so Sirocco can find its
atomic data:

```bash
cd runs/<name>
ln -s /path/to/sirocco/xdata data
ln -s /path/to/sirocco/zdata zdata
```

Sirocco runs in parallel via MPI: you launch many copies of the program
("ranks") that share the photon work. A fresh run is started like this
(adapt the core count and paths; on a SLURM cluster the same command goes
inside a job script, launched with `srun` instead of `mpirun`):

```bash
cd runs/<name>
mpirun -np 384 /path/to/sirocco/bin/sirocco -p outflow.pf 2>&1 | tee outflow.log
```

The `-p` flag ramps the photon number up logarithmically over the early
cycles (much faster for these optically thick winds). Use it for fresh runs
only — never when restarting from a previous run (`System_type=previous`).

Production settings, recorded per run in the `.slurm` files:
2×10⁶ photons/cycle, `matrix_pow` ionization, photon sampling `user_bands`
(5 bands, 0.05–100 eV, boundaries 1.5/3.4/7/13.6 eV) for the mu family and
`wide` for the ptf family, macro-atom H/He + simple-atom metals, 14 ionization + 20 (mu) or 40 (ptf)
spectrum cycles; 384 MPI ranks at 4 GB each (4 nodes) for 100-cell runs, or
192 ranks at 8 GB for the 300-cell run (matrix-mode memory scales with cell
count: ~1.8 GB/rank at 100 cells, ~5.4 at 300).

A completed run leaves behind `outflow.wind_save` — Sirocco's binary
snapshot of the converged wind state. The restarts, truncated reruns, and
Hα reruns below all start from a copy of this file.

## Post-processing products the notebooks read

From a completed run directory:

```bash
/path/to/sirocco/bin/windsave2table -xall outflow
```

produces `outflow.master.txt`, `outflow.heat.txt`,
`outflow.H_1.levden.txt` (patch 0002), `outflow.xspec.all.txt`
(per-cell J_ν, patch 0003), and `outflow.converge.txt` (per-cell
convergence flags, used below to place the truncation radius). Emergent spectra are `outflow.log_spec`
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
# in the <name>_trunc_speconly dir (data/zdata symlinks needed here too):
cp ../<name>/outflow.wind_save full.wind_save
/path/to/sirocco/bin/modify_wind -rcut <RCUT_CM> -out_root cut full
#   (floors ρ, n_e beyond RCUT -> cut.wind_save)
mpirun -np 384 /path/to/sirocco/bin/sirocco spec.pf 2>&1 | tee spec.log
```

Note there is **no `-p`** on this run (it is a restart), and the parameter
file is `spec.pf`, not `outflow.pf`: it sets `System_type=previous`,
`Wind.old_windfile=cut`, `Ionization_cycles=0`, and `Spectrum_cycles=400`
(the truncated winds are photon-starved in the optical, so they get 400
spectrum cycles vs the parents' 20/40). Its outputs are therefore named
`spec.*` — the notebooks read `spec.log_spec` from these directories.
`spec.slurm` records the cluster invocation used for the paper.

For a trunc dir's `halpha_hires/` rerun, copy `cut.wind_save` into that
subdirectory first; for a parent's `halpha_hires/`, copy the parent's
`outflow.wind_save` instead (its `.pf` reads `Wind.old_windfile outflow`).

RCUT is set at the outer edge of the contiguously converged core:
run `windsave2table` on the parent (see "Post-processing" above), read the
`converge` column of `outflow.converge.txt` (0 = converged), and place the
cut at the first cell
of the sustained unconverged tail (≥5 consecutive failed cells); isolated
interior failures (runs of <5 cells) are tolerated:

| run | RCUT | criterion | cells kept |
|---|---|---|---|
| mu_m15 | 8.509e18 cm | converged-core end (sustained tail from cell 76) | 0–75 of 99 |
| ptf_m15_n300 | 5.245e17 cm | converged-core end (sustained tail from cell 137); keeps the whole recombination front (r = 3.59e17) | 0–136 of 299 |

Truncation sensitivity (quantified in `FigAppTruncatedSpectraComparison`):
for mu_m15, whose cut sits essentially at its electron-scattering photosphere
(0.83 r_ph), the truncated and full-domain spectra agree closely redward of
the Balmer edge; for ptf_m15, whose cut sits inside its photosphere
(0.26 r_ph), the removed unconverged tail absorbs up to ~half the optical,
so the truncated (converged-cells-only) spectrum is the production statement. Everything blueward of the Balmer edge depends on the
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
The notebooks locate them via the `SIROCCO_REPRO_DATA` environment variable —
see README Part 1 for the full figure-reproduction walkthrough. The run
directories live under `<data root>/repro/<name>` with the same names as
`runs/<name>`; if you re-ran the simulations yourself, arrange your output
directories the same way and point `SIROCCO_REPRO_DATA` at them.
