#!/bin/bash
#PBS -l select=1:ncpus=36:mem=16gb
#PBS -l walltime=00:20:00
#PBS -N simba_test

module purge
module load miniforge/3

eval "$(~/miniforge3/bin/conda shell.bash hook)"
conda activate simba_env

cd $HOME/auto_react_mech_construct

python main.py

