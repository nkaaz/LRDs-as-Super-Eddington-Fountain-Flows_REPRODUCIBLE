#!/bin/bash
#SBATCH --job-name=ptfHA_m15n300p
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

# Spectrum-only extraction from the CONVERGED CORE of ptf_m15_n300 (the resolved-front
# 300-cell run). rcut=5.245e17 = start of the sustained unconverged tail (cell 137);
# keeps the entire steady recombination front (cell 120, r=3.59e17) and its converged
# wake; interior single-cell blips retained per convention. 400 spectrum cycles.
# 48 ranks/node x 8G for the 300-cell plasma arrays. No -p, no -r.
cd /scratch/gpfs/nk6352/sirocco_runs/repro/ptf_m15_n300/halpha_hires
srun sirocco outflow.pf 2>&1 | tee outflow.log
