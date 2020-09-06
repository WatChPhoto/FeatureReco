#!/bin/bash
# 
# Get list of unique bursts of images, and call averaging script
#

for i in `ls -1 B*.JPG | awk '{print substr($1,2,3)}' |uniq -d`;
do
    echo "Processing image ${i}"
    ./image_averaging.sh ${i}
    echo "Done image ${i}"
done

echo "Done."
