#!/bin/bash
for run in `cat runlist.txt`
do
    echo "submitting job ${run}"
    sbatch process_image_${run}.sh
done
squeue -u $USER

