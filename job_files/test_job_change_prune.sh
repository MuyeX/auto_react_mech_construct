#!/bin/bash
#PBS -l select=1:ncpus=64:mem=32gb
#PBS -l walltime=02:00:00
#PBS -N simba_test_change_prune

module purge
module load miniforge/3

eval "$(~/miniforge3/bin/conda shell.bash hook)"
conda activate simba_env

cp -r $HOME/auto_react_mech_construct $TMPDIR/

cd $TMPDIR/auto_react_mech_construct

start=`date +%s`

python main.py

end=`date +%s`
runtime=$((end-start))
echo $runtime

# Copy the results back to the home directory
cp -r -u $TMPDIR/auto_react_mech_construct $HOME/