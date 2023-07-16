#!/bin/bash


#first command line argument is name of the run

##### options setup #############################
DIR_MEAS=measurements
DIR_RUN=${DIR_MEAS}/$1
DIR_OUT=outputs
DIR_OUT_RUN=outputs/$1

N_CELLS="4 8 16 32 64 128 256 512 1024 2028 4096"

#SBATCH -n 1
#SBATCH --time=8:00:00
#SBATCH --mem-per-cpu=4000
#SBATCH --tmp=4000                       
#SBATCH --job-name=$1
#SBATCH --output=${DIR_OUT}
#SBATCH --error=error.err
#SBATCH --constraint=EPYC_7H12


EXEC_NAME="./projects.dgfem.quadratic_dirichlet"

INV=0.5
SIGMA=20

INV_MAX=2.0
SIGMA_MAX=100

INV_STEP_SIZE=0.2
SIGMA_STEP_SIZE=5

envl2mod
module load gcc/9.3.0
module load cmake/3.20.3

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


#SUBMIT
sbatch --wrap=${EXEC_NAME} $RUN_NAME $C_INV $C_SIGMA $N_CELLS