#!/bin/sh
#PBS -l nodes=36:ppn=4

module add openmpi/gcc/64/1.10.1

cd $PBS_O_WORKDIR
num=` cat $PBS_NODEFILE | wc -l`
uniq $PBS_NODEFILE

mpiexec -n $num -machinefile $PBS_NODEFILE ./sub_convert_milc > quick.out

