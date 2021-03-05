#!/bin/bash
#SBATCH -p defq
#SBATCH -n 28
#SBATCH --ntasks-per-node=28
#SBATCH -J alpscthyb
#SBATCH -o stdout.%J
#SBATCH -e stderr.%J

export DCORE_MPIRUN_COMMAND="mpirun -np 28"
export DCORE_ALPSCTHYB_TIMELIMIT=60
export ROOTDIR=$HOME/repos/benchmarks
export COMPUTE_2PARTICLE=1

for name in 'SIAM_Discrete_Bath' ; do
    cd $ROOTDIR/$name/scripts && ./alps_cthyb > output && cd ..
done

#'Sr2RuO4' does not work because it requires MPI in generating the TB model
#for name in 'Dimer' ; do
    #cd $ROOTDIR/$name/scripts && ./alps_cthyb > output && cd ..
#done

#for name in 'Hubbard_Atom' 'Dimer' 'Dimer_SOC' 'Trimer'; do
    #cd $ROOTDIR/$name/scripts && ./alps_cthyb > output && cd ..
#done
