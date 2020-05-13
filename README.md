# FeatureReco

Code to do feature recognition to locate bolts and PMT centers.  Implemented using OpenCV libraries.

You will need a copy of opencv libraries installed.  To build the code:

cd FeatureReco
cmake .
make

The initial version only has some code to play with OpenCV to remove image distortion.  For example running:

./DisplayImage /bigdisk/jamieson/TOW-Feb2020/BarrelSurveyFar/B0170163.JPG

opens the input file, displays it, along with a black-and-white version that has the distortion removed.  The undistored image is written to output.jpg

Future improvements:
1) Output to filename based on input name
2) Start looking for features and output them to text files
3) Cleanup unused code


