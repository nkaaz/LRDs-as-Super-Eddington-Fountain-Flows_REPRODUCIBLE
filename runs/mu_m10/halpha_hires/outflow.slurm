#!/bin/bash
#SBATCH --job-name=muHA_m10
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=96
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --time=08:00:00
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=nk6352@princeton.edu

module purge
module load gcc-toolset/13 openmpi/gcc/4.1.6 gsl/2.6

export SIROCCO=/home/nk6352/sirocco/sirocco
export PATH=$SIROCCO/bin:$PATH

# Halpha-only spectrum restart off the parent's frozen ionization (System_type=previous).
# No -p (that's for fresh-run photon ramp; this is a restart).
cd /scratch/gpfs/nk6352/sirocco_runs/repro/mu_m10/halpha_hires

srun sirocco outflow.pf 2>&1 | tee outflow.log
