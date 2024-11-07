#!/bin/tcsh
#SBATCH --partition=short-serial
#SBATCH --job-name=massplot_PPIs
#SBATCH --output=/home/users/eesjir/eesjir_group_workspace/transfer2_data_prep/scripts/slurmout/slurm-%j.out
#SBATCH --error=/home/users/eesjir/eesjir_group_workspace/transfer2_data_prep/scripts/slurmout/slurm-%j.err
#SBATCH --array=1-6200%200
#SBATCH --time=05:00:00

setenv PATH ${PATH}:/gws/smf/j04/ncas_radar/software/miniconda3_radar_group_20200519/envs/R_4_10_biorad_pyart_3_8/bin/

conda init tcsh

conda activate R_4_10_biorad_pyart_3_8

set scripts = /home/users/eesjir/eesjir_group_workspace/transfer2_data_prep/scripts/
set DIR = $1

#the x-band_PPI_massprocessor_indexgen.py must have been run first or the following will fail!

cd $DIR
set files = `ls *_XYZ_raw.nc`
@ nfiles = `ls *_XYZ_raw.nc | wc -l`

@ ST = 0
@ ST += $SLURM_ARRAY_TASK_ID

if($ST <= $nfiles) then
	python $scripts/x-band_PPI_massprocessor_NeelyPresences.py $DIR ${files[${ST}]}
endif