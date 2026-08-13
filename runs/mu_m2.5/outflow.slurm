#!/bin/bash
#SBATCH --job-name=mu_m2.5            # marginally-unbound Mdot=2.5; matrix_pow + custombands; NBINS=3000
#SBATCH --nodes=4                  # node count (full nodes)
#SBATCH --ntasks-per-node=96       # full-node packing (96 cores/node on Stellar)
#SBATCH --cpus-per-task=1          # cpu-cores per task
#SBATCH --mem-per-cpu=4G           # memory per cpu-core (matrix mode ~1.8G/rank)
#SBATCH --time=48:00:00            # 10%Z + fixed MAXSCAT; single-shot ion+spec
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=nk6352@princeton.edu
module purge
module load gcc-toolset/13 openmpi/gcc/4.1.6 gsl/2.6
export SIROCCO=/home/nk6352/sirocco/sirocco
export PATH=$SIROCCO/bin:$PATH
cd /scratch/gpfs/nk6352/sirocco_runs/repro/mu_m2.5
srun sirocco -p outflow.pf 2>&1 | tee outflow.log
