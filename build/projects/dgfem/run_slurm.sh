#!/bin/bash

#SBATCH -n 1
#SBATCH --time=8:00:00
#SBATCH --mem-per-cpu=4000
#SBATCH --tmp=4000                       
#SBATCH --constraint=EPYC_7H12


#first command line argument is name of the run

##### options setup ####################

EXEC_NAME="projects.dgfem.dirichlet_$EXEC_NUM"

N_CELLS="4 8 16 32 64 128 256 512 1024 2048 4096"

echo "command: ./${EXEC_NAME} $RUN_NAME $INV $SIGMA $N_CELLS"

./${EXEC_NAME} $RUN_NAME $INV $SIGMA $N_CELLS

