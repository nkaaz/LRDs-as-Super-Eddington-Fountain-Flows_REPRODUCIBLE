#!/bin/bash
#SBATCH --job-name=ptfHA_m15n300t
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=48
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8G
#SBATCH --time=8:00:00
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=nk6352@princeton.edu

module purge
module load gcc-toolset/13 openmpi/gcc/4.1.6 gsl/2.6
export SIROCCO=/home/nk6352/sirocco/sirocco
export PATH=$SIROCCO/bin:$PATH

# Halpha-only (6500-6700 A, 10^4 bins) spectrum extraction from the truncated
# (converged-core) wind of ptf_m15_n300: System_type=previous reads cut.wind_save
# (copy it here from the trunc dir). 20 spectrum cycles; all photons generated
# in-band. No -p, no -r.
cd /scratch/gpfs/nk6352/sirocco_runs/repro/ptf_m15_n300_trunc_speconly/halpha_hires
srun sirocco outflow.pf 2>&1 | tee outflow.log
