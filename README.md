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

./FindBoltLocations [input-image-file-name-with-path] [bolt-location-textfile-with-path]
<pre>
Output file names are: 
*  gausblur[input-file-name]  -- image with gaussian applied (if enabled) 
*  bifilter[input-file-name]  -- image after bilateral filter applied (if enabled)
*  sobel[input-file-name]     -- image after sobel filter applied (this is input to Hough)
*  circles[input-file-name]   -- original image with circles found by hough added
*  FindBoltLocation.root      -- root histograms
</pre>
To do list:
Updates to output michel.jpg
* Make red and green circles 1 pixel wide                                             ---Completed
* Add a blue line from the matches between red and green circles (for distance < some threshold)
Updtes to histogram_ch0.root                                                          ---Completed. Added black line instead for visiblity reason.
* Rename output root file to FindBoltLocations.root                                   ---Completed
* make histogram of distance from each green circle (instead of from each red circle) ---Completed.
* update range of histogram... distance can only be positive 
Pick a run to use to do optimization on, and check optimization histogram is okay.
Optimize mean of above histogram by changing filter parametres and hough parameters
* write script to run many different parameter variations (eg. in bash, sed to edit Config.txt, awk to parse a string; or use python)
* each run done in a seperate folder, and save config file to that folder
* write a root macro to collect results from the parameter scans










