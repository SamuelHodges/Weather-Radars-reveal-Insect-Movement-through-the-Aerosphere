#!/bin/tcsh

@ start = $1
@ end = $2

set Dir = $3

set scripts = /home/users/eesjir/eesjir_group_workspace/transfer2_data_prep/scripts
cd $Dir

while ($start <= $end)
    set Files = `ls *${start}*SUR_v1.nc | grep -v raw`

    foreach var ($Files)
        sbatch $scripts/convert_to_lonlatalt.csh $var $Dir
#        $scripts/convert_to_lonlatalt.csh $var $Dir
#        exit(1)
    end
    @ start ++
end
