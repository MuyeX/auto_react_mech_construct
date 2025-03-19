#!/bin/bash
#PBS -l select=1:ncpus=32:mem=16gb
#PBS -l walltime=00:20:00
#PBS -N simba_test

module purge
module load miniforge/3

eval "$(~/miniforge3/bin/conda shell.bash hook)"
conda activate simba_env

cp -r $HOME/auto_react_mech_construct/ $TMPDIR/auto_react_mech_construct

cd $TMPDIR/auto_react_mech_construct

start=`date +%s`

python main.py config_hypo.json

end=`date +%s`
runtime=$((end-start))
echo $runtime

# Copy the results back to the home directory
cp -r -u $TMPDIR/auto_react_mech_construct $HOME/auto_react_mech_construct