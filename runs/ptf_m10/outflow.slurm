#!/bin/bash
#SBATCH --job-name=ptf_m10
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=96
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --time=32:00:00            # f0.94 (G1.213) winds; same recipe as pt_*
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=nk6352@princeton.edu
module purge
module load gcc-toolset/13 openmpi/gcc/4.1.6 gsl/2.6
export SIROCCO=/home/nk6352/sirocco/sirocco
export PATH=$SIROCCO/bin:$PATH
cd /scratch/gpfs/nk6352/sirocco_runs/repro/ptf_m10
srun sirocco -p outflow.pf 2>&1 | tee outflow.log
