#!/bin/bash
#SBATCH -p defq
#SBATCH -n 28
#SBATCH --ntasks-per-node=28
#SBATCH -J pomerol
#SBATCH -o stdout.%J
#SBATCH -e stderr.%J

export ROOTDIR=$HOME/repos/benchmarks
export COMPUTE_2PARTICLE=1

#for name in 'Dimer' ; do
#for name in 'Hubbard_Atom' ; do
for name in 'SIAM_Discrete_Bath' ; do
    cd $ROOTDIR/$name/scripts && mpirun -np 28 ./pomerol > output-pomerol && cd ..
done
