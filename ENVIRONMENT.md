# Radiative-transfer environment (Sirocco)

The synthetic spectra and wind diagnostics in this paper were produced
with **Sirocco** (https://github.com/sirocco-rt/sirocco). We used
commit `5aef5f17` plus six small local patches provided in
`sirocco_patches/`. 

## Instructions

The **Sirocco** modifications ship in this repo's `sirocco_patches/` folder as numbered **patch files**. The commands below download the public code, switch to the correct version, and apply the patches. 

Prerequisites: a Linux machine with git, an MPI C compiler (e.g. OpenMPI's
`mpicc`), and the GSL library. On the Princeton Stellar cluster these come from
`module load gcc-toolset/13 openmpi/gcc/4.1.6 gsl/2.6`.

**Step 1. Download the public Sirocco source and switch to the version used here.**

```bash
git clone https://github.com/sirocco-rt/sirocco.git
cd sirocco
git checkout 5aef5f17
```

**Step 2. Apply the six patches.**

```bash
git am /path/to/this-repo/sirocco_patches/*.patch
```

**Step 3. Build.**

```bash
./configure
cd source
make -k sirocco windsave2table modify_wind
```

Success means the three executables now exist:

```bash
ls ../bin/sirocco ../bin/windsave2table ../bin/modify_wind
```

## The six patches

| # | Patch | Why it is required |
|---|---|---|
| 0001 | `Central_object.dilution_factor` | New `.pf` parameter: the central source emits a blackbody spectrum at `Central_object.temp` (the color temperature) with luminosity multiplied by W, so L = W·4πR²σT⁴. Every run injects a diluted interior radiation field (W = 1/τ_es(r_in)). |
| 0002 | levden tables in `windsave2table` | Writes `rootname.<Elem>_<istate>.levden.txt` (NLTE level populations).|
| 0003 | `NBINS_IN_CELL_SPEC` 1000→3000 | Per-cell J_ν resolution of `windsave2table -xall` (the `outflow.xspec.all.txt` files the notebooks read have 3001 rows). |
| 0004 | `MAXSCAT` 2000→200000 | Prevents photons from being destroyed in our runs due to reaching a maximum scattering limit. |
| 0005 | `xdata/h20_hetop_standard80_z0.1.dat` | The 10%-metallicity atomic dataset (`Atomic_data data/h20_hetop_standard80_z0.1.dat`). |
| 0006 | `modify_wind -rcut` | Builds the truncated wind-saves used by the spectrum-only reruns (see below). |

## Running a model

Each `runs/<name>/` contains all necessary input: `outflow.pf` (the
parameter file), `outflow.import.txt` (the 1-D spherical wind model; see
`winds/` for how these are generated), and the `.slurm` script. Before
running, the run directory needs two symbolic links so Sirocco can find the
atomic data:

```bash
cd runs/<name>
ln -s /path/to/sirocco/xdata data
ln -s /path/to/sirocco/zdata zdata
```

Run Sirocco via MPI using a command like:

```bash
cd runs/<name>
mpirun -np 384 /path/to/sirocco/bin/sirocco -p outflow.pf 2>&1 | tee outflow.log
```

The `-p` flag ramps the photon number up over the early
cycles (this speeds up these optically thick wind runs). 

Some of our production settings are:
2×10⁶ photons/cycle, `matrix_pow` ionization, photon sampling `user_bands`
(5 bands, 0.05–100 eV, boundaries 1.5/3.4/7/13.6 eV) for the mu family and
`wide` for the ptf family, macro-atom H/He + simple-atom metals, 14 ionization + 20 (mu) or 40 (ptf)
spectrum cycles; 384 MPI ranks at 4 GB each (4 nodes) for 100-cell runs, or
192 ranks at 8 GB for the 300-cell run.

## Post-processing simulation output for the notebooks

From a completed run directory:

```bash
/path/to/sirocco/bin/windsave2table -xall outflow
```

produces `outflow.master.txt`, `outflow.heat.txt`,
`outflow.H_1.levden.txt`, `outflow.xspec.all.txt`, and `outflow.converge.txt`. Emergent spectra are `outflow.log_spec`.

## Truncated spectrum-only reruns (`*_trunc_speconly`)

In the highest-Ṁ runs the wind beyond the hydrogen recombination front
struggles to converge (the gas there decouples from the radiation field and t_e
collapses); this doesn't matter for the interior solution but
can contaminate the spectra. So, we also calculated spectra
from spectrum-only reruns on a truncated wind where the outer unconverged region is excised:

```bash
# in the <name>_trunc_speconly dir (data/zdata symlinks needed here too):
cp ../<name>/outflow.wind_save full.wind_save
/path/to/sirocco/bin/modify_wind -rcut <RCUT_CM> -out_root cut full
#   (floors ρ, n_e beyond RCUT -> cut.wind_save)
mpirun -np 384 /path/to/sirocco/bin/sirocco spec.pf 2>&1 | tee spec.log
```

Its outputs are named `spec.*`.

For a trunc dir's `halpha_hires/` rerun, we copy `cut.wind_save` into that
subdirectory first; for a parent's `halpha_hires/`, we copy the parent's
`outflow.wind_save` instead (its `.pf` reads `Wind.old_windfile outflow`).

RCUT is set at the outer edge of the contiguously converged core. This can be read from the
`converge` column of `outflow.converge.txt` (0 = converged). We place the
cut at the first cell of the unconverged region. We ignore the odd interior cell that is unconverged.
