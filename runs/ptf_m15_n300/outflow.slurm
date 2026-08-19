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

# RADIAL-RESOLUTION TEST: ptf_m15 at 300 cells (vs production 100). Same physics
# (f=0.94, Mdot=15; r_in/T_core/r_out identical) — tests whether the recombination-front
# limit cycle (bounded oscillation, conv ~35+-5, front sloshing cell 38<->42) is a grid
# artifact: front cells are internally thick to BaC (tau_BaC/cell ~34 at n100 -> ~11 here).
# Readout: does the oscillation amplitude shrink with resolution?
# 48 ranks/node x 8G: matrix mode needs ~5.4 GB/rank at 300 cells (1.8 GB @ 100).
# Fresh run -> -p ramp-up per family recipe.
cd /scratch/gpfs/nk6352/sirocco_runs/repro/ptf_m15_n300
srun sirocco -p outflow.pf 2>&1 | tee outflow.log
