#!/bin/bash
#SBATCH --ntasks=48
#SBATCH --time=96:0:0
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48
#SBATCH --mem-per-cpu=4096m

ulimit -s unlimited
set -e

module purge; module load bluebear
module load bear-apps/2021b
module load LAMMPS/23Jun2022-foss-2021b-kokkos

export OMP_PROC_BIND=false
export OMP_PLACES=threads
export OMP_NUM_THREADS=1
mpiexec -np 48 lmp -k on -sf kk -in NO5M_inf4.in
