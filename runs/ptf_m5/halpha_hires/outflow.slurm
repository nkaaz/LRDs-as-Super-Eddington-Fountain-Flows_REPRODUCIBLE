#!/bin/bash
#SBATCH --job-name=ptfHA_m5p
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

# Halpha-only (6500-6700 A) spectrum restart off the FULL-DOMAIN parent ionization of
# ptf_m5 (outflow.wind_save copied from the parent; System_type=previous, no -p, no -r).
# Untruncated counterpart to <trunc>/halpha_hires; Halpha is tail-insensitive (no n=2
# population in the cold neutral outer wind) so this and the trunc version should agree.
cd /scratch/gpfs/nk6352/sirocco_runs/repro/ptf_m5/halpha_hires
srun sirocco outflow.pf 2>&1 | tee outflow.log
