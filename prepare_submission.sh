#!/bin/bash
for run in `cat runlist.txt`
do
    echo "preparing ${run}"
    cp process_image_XXX.sh process_image_${run}.sh
    sed -i "s/XXX/${run}/g" process_image_${run}.sh
done
    
