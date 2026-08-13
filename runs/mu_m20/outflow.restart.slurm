#!/bin/bash
#SBATCH --job-name=mu_m20_r
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=96
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --time=48:00:00            # spectrum-only restart: 5/20 cycles took ~most of 18h -> 15 left
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=nk6352@princeton.edu
module purge
module load gcc-toolset/13 openmpi/gcc/4.1.6 gsl/2.6
export SIROCCO=/home/nk6352/sirocco/sirocco
export PATH=$SIROCCO/bin:$PATH
cd /scratch/gpfs/nk6352/sirocco_runs/repro/mu_m20
# -r = restart from wind_save+spec_save, continue spectrum cycles 6..20 (NOT -p). pf unchanged.
srun sirocco -r outflow.pf 2>&1 | tee outflow.restart.log
