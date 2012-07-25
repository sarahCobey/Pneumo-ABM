#!/bin/bash
 
j=1

for treatment in `cat treatments.txt`; do
  for replicate in `cat simulations.txt`; do

    bsub -oo tr_${j}_sim_${replicate}_screen -eo tr_${j}_sim_${replicate}_err -q queue_name "echo \"$j $treatment $replicate\" | ./simulation.out"

  done
  let j=$j+1 
done


