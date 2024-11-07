#!/bin/tcsh
#SBATCH --partition=short-serial
#SBATCH --job-name=x-band
#SBATCH --output=/home/users/eesjir/eesjir_group_workspace/transfer2_data_prep/scripts/slurmout/slurm-%j.out
#SBATCH --error=/home/users/eesjir/eesjir_group_workspace/transfer2_data_prep/scripts/slurmout/slurm-%j.err
#SBATCH --time=01:00:00

setenv PATH ${PATH}:/gws/smf/j04/ncas_radar/software/miniconda3_radar_group_20200519/envs/R_4_10_biorad_pyart_3_8/bin/

conda init tcsh

conda activate R_4_10_biorad_pyart_3_8
cd /home/users/eesjir/Chilbolton

set scripts = /home/users/eesjir/eesjir_group_workspace/transfer2_data_prep/scripts/
set date = $1
Rscript $scripts/copy_instantaneous_hourly_radar_data.r $date
