#!/bin/bash
#
# average images starting with B${1}*.jpg
# crop them to make images C${1}.jpg
# averaged image is ${1}.jpg

image_num=${1}

# First removed blurry photos by eye!
# resulting in usable photos listed in ListOfPhotos.txt


# Crop out the top row of drone watermark
for f in B${image_num}*.JPG;
do
    convert $f -crop 4000x2750+0+250 ${f/B/C};
done

# Align the images
align_image_stack -c 20 -v --corr=0.5 -m -a A${image_num} C${image_num}*.JPG
            
# Combine images taking the median pixel value
convert A${image_num}* -evaluate-sequence median ${image_num}.jpg

# Delete cropped images, and aligned images
rm A${image_num}*
rm C${image_num}*

