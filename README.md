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

Still needs lots of work.  So far it reads in an input file and applies a Gaussian blur, bilateral filter, and Sobel transform.  Outputs three images: one with with the Gaussian blur, one with bilateral filter, and third with bilateral filter follwed by Sobel.

./FindBoltLocations /bigdisk/jamieson/TOW-Feb2020/BarrelSurveyFar/B0170163.JPG 

Output file names are gausblur[input-file-name] and bifilter[input-file-name]

To do list:
* make a config file to read in parameters of bluring, and which filtering to apply
* add code to find features (Canny, Prewitt??)
* think about how to find the bolts.







