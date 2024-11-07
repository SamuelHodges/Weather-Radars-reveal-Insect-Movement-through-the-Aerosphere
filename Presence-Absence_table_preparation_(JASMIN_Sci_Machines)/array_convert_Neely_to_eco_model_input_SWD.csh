#!/bin/tcsh
#SBATCH --partition=short-serial
#SBATCH --job-name=patable
#SBATCH --output=/home/users/eesjir/eesjir_group_workspace/transfer2_data_prep/scripts/slurmout/slurm-%j.out
#SBATCH --error=/home/users/eesjir/eesjir_group_workspace/transfer2_data_prep/scripts/slurmout/slurm-%j.err
#SBATCH --array=0-31
#SBATCH --time=23:00:00
#SBATCH --mem=5G

setenv PATH ${PATH}:/gws/smf/j04/ncas_radar/software/miniconda3_radar_group_20200519/envs/R_4_10_biorad_pyart_3_8/bin/

conda init tcsh

conda activate R_4_10_biorad_pyart_3_8
cd /home/users/eesjir/eesjir_group_workspace

set scripts = /home/users/eesjir/eesjir_group_workspace/transfer2_data_prep/scripts/
# date should be numeric in yyyymm with a '01' appended, eg. 20170701 for July 2017
set date=$1
set hour1=$2
set hour2=$3
set targetdir=$4

cd $targetdir

@ ST = $date
@ ST += $SLURM_ARRAY_TASK_ID

Rscript $scripts/convert_lonlatalt_wradlib_to_eco_model_input_SWD_variant_type2.5_countsonly_Neely_cmd.r ${ST} $hour1 $hour2

