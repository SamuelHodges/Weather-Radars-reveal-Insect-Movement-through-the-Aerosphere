#!/bin/tcsh
#SBATCH --partition=short-serial
#SBATCH --job-name=x-band
#SBATCH --output=/home/users/eesjir/eesjir_group_workspace/transfer2_data_prep/scripts/slurmout/slurm-%j.out
#SBATCH --error=/home/users/eesjir/eesjir_group_workspace/transfer2_data_prep/scripts/slurmout/slurm-%j.err
#SBATCH --array=250-2000%200
#SBATCH --time=05:00:00
#SBATCH --mem=10G

setenv PATH ${PATH}:/gws/smf/j04/ncas_radar/software/miniconda3_radar_group_20200519/envs/R_4_10_biorad_pyart_3_8/bin/

conda init tcsh

conda activate R_4_10_biorad_pyart_3_8
cd /home/users/eesjir/eesjir_group_workspace

set targetdir=$1
cd targetdir


set files = `ls *SUR_v1.nc | grep -v raw`
@ nfiles = `ls *SUR_v1.nc | grep -v raw | wc -l`

set scripts = /home/users/eesjir/eesjir_group_workspace/transfer2_data_prep/scripts/

@ ST = 0
@ ST += $SLURM_ARRAY_TASK_ID

if($ST <= $nfiles) then
	python $scripts/Copypython_convert_xband_to_NeelyPres.py ${files[${ST}]} $targetdir
endif
