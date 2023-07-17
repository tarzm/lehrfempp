#!/bin/bash

for EXEC_NUM in $(seq 3 5)
do

#first command line argument is name of the run
#second is c_inv
#third is c_sigma

##### options setup #############################

RUN_NAME="${1}_$EXEC_NUM"

DIR_MEAS=measurements
DIR_RUN=${DIR_MEAS}/$RUN_NAME
DIR_OUT=outputs
DIR_OUT_RUN=outputs/$RUN_NAME/output

INV=$2
SIGMA=$3

export RUN_NAME
export DIR_OUT_RUN
export EXEC_NUM
export INV
export SIGMA

#####END options setup #############################

#create measurements dir if it does not exist
if [ ! -d "$DIR_MEAS" ]; then
    echo "dir does not exists"
    mkdir $DIR_MEAS
fi

#SUBMIT
sbatch < run_slurm.sh
done


