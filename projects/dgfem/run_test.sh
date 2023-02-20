#!/bin/bash

#first command line argument is name of the run

##### options setup #############################
DIR_MEAS=measurements
DIR_RUN=${DIR_MEAS}/$1
DIR_OUT=outputs
DIR_OUT_RUN=outputs/$1

#####END options setup #############################

#make the exectuable
make projects.dgfem.polytopic_from_hybrid

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

#create output dir if it does not exist
if [ ! -d "$DIR_OUT" ]; then
    mkdir $DIR_OUT
fi

#create output name directoy in measurements if it does not exist
#if it does, remove everything inside it
if [ -d "$DIR_OUT_RUN" ]; then
    rm -rf ${DIR_OUT_RUN}
fi

if [ ! -d "$DIR_OUT_RUN" ]; then
    mkdir $DIR_OUT_RUN
fi

#run test
./projects.dgfem.polytopic_from_hybrid $1