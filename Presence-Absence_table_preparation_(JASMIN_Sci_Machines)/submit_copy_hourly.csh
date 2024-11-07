#!/bin/tcsh

@ start = $1
@ end = $2

while ($start <= $end)
  sbatch ./copy_hourly.csh $start
  @ start ++
end

