#!/bin/bash
#SBATCH --job-name=ptf_m10_speconly
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=96
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --time=16:00:00
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=nk6352@princeton.edu

module purge
module load gcc-toolset/13 openmpi/gcc/4.1.6 gsl/2.6
export SIROCCO=/home/nk6352/sirocco/sirocco
export PATH=$SIROCCO/bin:$PATH

# Spectrum-only extraction from the converged inner wind of ptf_m10 (f=0.94).
# Outer unconverged neutral tail (cell centre r > 4.111e+17 cm = 1.35x recomb front,
# cells 50-98) made transparent via modify_wind -rcut -> cut.wind_save.
# System_type=previous reads cut.wind_save directly (no -p, no -r).
cd /scratch/gpfs/nk6352/sirocco_runs/repro/ptf_m10_trunc_speconly
srun sirocco spec.pf 2>&1 | tee spec.log
