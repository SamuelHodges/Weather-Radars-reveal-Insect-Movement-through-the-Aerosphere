#!/bin/tcsh

@ start = $1
@ end = $2

set hour1 = $3
set hour2 = $4
set Dir = $5

set scripts = /home/users/eesjir/eesjir_group_workspace/transfer2_data_prep/scripts
cd $Dir

while ($start <= $end)
    sbatch $scripts/convert_to_eco_model_input_SWD.csh $start $hour1 $hour2 $Dir
    @ start ++
end
