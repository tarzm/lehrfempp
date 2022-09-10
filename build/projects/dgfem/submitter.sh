#!/bin/bash


#first command line argument is name of the run

##### options setup #############################
DIR_MEAS=measurements
DIR_RUN=${DIR_MEAS}/$1
DIR_OUT=outputs
DIR_OUT_RUN=outputs/$1

N_CELLS="4 8 16 32 64 128 256 512 1028"
MODEL_COMMAND="-R 'select[model=EPYC_7742]'"
EXEC_NAME=./projects.dgfem.quadratic_dirichlet

INV_MAX=2.0
SIGMA_MAX=100

INV_STEP_SIZE=0.2
SIGMA_STEP_SIZE=5
#####END options setup #############################

#create measurements dir if it does not exist
if [ ! -d "$DIR_MEAS" ]; then
    echo "dir does not exists"
    mkdir $DIR_MEAS
fi

#create run name directoy in measurements if it does not exist
#if it does, remove everything inside it
if [ -d "$DIR_RUN" ]; then
    rm -rf ${DIR_RUN}
fi

if [ ! -d "$DIR_RUN" ]; then
    mkdir $DIR_RUN
fi

#create measurements dir if it does not exist
if [ ! -d "$DIR_OUT" ]; then
    mkdir $DIR_OUT
fi

#create run name directoy in measurements if it does not exist
#if it does, remove everything inside it
if [ -d "$DIR_OUT_RUN" ]; then
    rm -rf ${DIR_OUT_RUN}
fi

if [ ! -d "$DIR_OUT_RUN" ]; then
    mkdir $DIR_OUT_RUN
fi


for C_INV in $(seq 1.0 0.3 5.0)
do
    for ((C_SIGMA = 30; C_SIGMA <= 170; C_SIGMA += 10))
    do

    bsub -n 1 -oo outputs/$1/${C_INV}_${C_SIGMA}.txt -W 04:00 -R 'select[model=EPYC_7742]' -J \'${C_INV}_${C_SIGMA}\'  ./quadratic $1 $C_INV $C_SIGMA $N_CELLS

    done

done






