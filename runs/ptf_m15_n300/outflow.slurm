#!/bin/bash
#SBATCH --job-name=ptf_m15_n300
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=48
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8G
#SBATCH --time=32:00:00
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=nk6352@princeton.edu

module purge
module load gcc-toolset/13 openmpi/gcc/4.1.6 gsl/2.6
export SIROCCO=/home/nk6352/sirocco/sirocco
export PATH=$SIROCCO/bin:$PATH

# PRODUCTION photon-tired Mdot=15 run (300 radial cells). The 100-cell version of
# this model does not converge -- its recombination front executes a grid-resolution
# limit cycle (front cells internally thick to the Balmer continuum, tau_BaC ~ 34/cell);
# at 300 cells (tau_BaC ~ 11/cell) the front converges monotonically. See
# ENVIRONMENT.md, "Radial resolution". 48 ranks/node x 8G (matrix-mode memory scales
# with cell count). Fresh run -> -p ramp-up per family recipe.
cd /scratch/gpfs/nk6352/sirocco_runs/repro/ptf_m15_n300
srun sirocco -p outflow.pf 2>&1 | tee outflow.log
