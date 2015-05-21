#!/bin/sh
#PBS -e node-037:/home/krylov/bsv/output/
#PBS -o node-037:/home/krylov/bsv/output/
LD_LIBRARY_PATH=/usr/lib64/openmpi/lib
cd /home/krylov/bsv/output && /usr/lib64/openmpi/bin/mpiexec -mca btl openib,self -machinefile $PBS_NODEFILE $PBS_O_WORKDIR/debug
