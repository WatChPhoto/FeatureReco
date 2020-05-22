# FeatureReco

Code to do feature recognition to locate bolts and PMT centers.  Implemented using OpenCV libraries.

You will need a copy of opencv libraries installed.  To build the code:

cd FeatureReco
cmake .
make


# Undistortion program

Simple program to take input JPG, apply undistortion, and write out undistoted image.

./UndistortImage /bigdisk/jamieson/TOW-Feb2020/BarrelSurveyFar/B0170163.JPG [output-file]

If output file is not specified, then it names output file undistorted<input-file-name>

# Feature finding code

Locate bolts in image by applying filters and hough transform.

./FindBoltLocations /bigdisk/jamieson/TOW-Feb2020/BarrelSurveyFar/B0170163.JPG 

Output file names are: 
  gausblur[input-file-name]  -- image with gaussian applied (if enabled) 
  bifilter[input-file-name]  -- image after bilateral filter applied (if enabled)
  sobel[input-file-name]     -- image after sobel filter applied (this is input to Hough)
  circles[input-file-name]   -- original image with circles found by hough added

To do list:
* add code to find features (Canny, Prewitt??)
* read in bolt locations from Patrick and Michael's output
* histogram shortest distance between each circle we find and closest bolt from Patrick's file
* optimize mean of above histogram by changing filter parametres and hough parameters.








