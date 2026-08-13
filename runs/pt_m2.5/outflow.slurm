#!/bin/bash
#SBATCH --job-name=pt_m2.5
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=96
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --time=48:00:00            # 10%Z + fixed MAXSCAT; single-shot ion+spec
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=nk6352@princeton.edu
module purge
module load gcc-toolset/13 openmpi/gcc/4.1.6 gsl/2.6
export SIROCCO=/home/nk6352/sirocco/sirocco
export PATH=$SIROCCO/bin:$PATH
cd /scratch/gpfs/nk6352/sirocco_runs/repro/pt_m2.5
srun sirocco -p outflow.pf 2>&1 | tee outflow.log
