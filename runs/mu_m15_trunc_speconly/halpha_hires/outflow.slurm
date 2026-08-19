#!/bin/bash
#SBATCH --job-name=muHA_m15t
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=96
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --time=8:00:00
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=nk6352@princeton.edu

module purge
module load gcc-toolset/13 openmpi/gcc/4.1.6 gsl/2.6
export SIROCCO=/home/nk6352/sirocco/sirocco
export PATH=$SIROCCO/bin:$PATH

# Spectrum-only extraction from the converged core of mu_m15. Sustained unconverged
# tail (cells 76/78-98; solid core ends at cell 75, interior blips 53-56/67-68 kept
# per family convention) made transparent via modify_wind -rcut 8.509e18 -> cut.wind_save.
# System_type=previous reads cut.wind_save directly (no -p, no -r). 20 spectrum cycles.
cd /scratch/gpfs/nk6352/sirocco_runs/repro/mu_m15_trunc_speconly/halpha_hires
srun sirocco outflow.pf 2>&1 | tee outflow.log
