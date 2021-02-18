#!/bin/bash
#SBATCH -p defq
#SBATCH -n 32
#SBATCH --ntasks-per-node=32
#SBATCH -J alpscthyb
#SBATCH -o stdout.%J
#SBATCH -e stderr.%J

export DCORE_MPIRUN_COMMAND="mpirun -np 32"
export DCORE_ALPSCTHYB_TIMELIMIT=60
export ROOTDIR=$HOME/repos/benchmarks

#'Sr2RuO4' does not work because it requires MPI in generating the TB model

for name in 'Hubbard_Atom' 'Dimer' 'Dimer_SOC' 'Trimer'; do
    cd $ROOTDIR/$name/scripts && ./alps_cthyb > output && cd ..
done
