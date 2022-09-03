#!/bin/bash

#name of the benchmarking run
export RUN_NAME=euler_1

#create necessary directory for outputs if it does not exist
mkdir measurements
mkdir measurements/${RUN_NAME}


N_CELLS="4 8 16 32 64 128"

BASE_COOMMAND="-n 1 -W 04:00 -R 'select[model=EPYC_7742]'"

EXEC_NAME="./projects.dgfem.quadratic_dirichlet"

C_INV=0.5
C_SIGMA=5

INV_MAX=2.0
SIGMA_MAX=100

INV_STEP_SIZE=0.2
SIGMA_STEP_SIZE=5

until [$C_INV -gt $INV_MAX]
do
    until [$C_SIGMA -gt $SIGMA_MAX]
    do

    bsub $BASE_COMMAND < $EXEC_NAME $RUN_NAME $C_INV $C_SIGMA $N_CELLS

    ((C_SIGMA += SIGMA_STEP_SIZE))
    done

((C_INV += INV_STEP_SIZE))
done






