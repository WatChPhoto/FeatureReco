# FeatureReco

Code to do feature recognition to locate bolts and PMT centers.  Implemented using OpenCV libraries.

You will need a copy of opencv libraries installed.  To build the code:

```
> cd FeatureReco
> cmake .
> make
```

# Undistortion program

Simple program to take input JPG, apply undistortion, and write out undistoted image.

```
./UndistortImage /bigdisk/jamieson/TOW-Feb2020/BarrelSurveyFar/B0170163.JPG [output-file]
```

If output file is not specified, then it names output file undistorted<input-file-name>

# Feature finding code

Locate bolts in image by applying filters and hough transform.

```
./FindBoltLocations [input-image-file-name-with-path] [opt=bolt-location-textfile-with-path]
```

# Config.txt
Configuration file to change parameters and set verbosity flag.

<pre>
Output file names are: 
*  gausblur[input-file-name]  -- image with gaussian applied (if enabled) 
*  bifilter[input-file-name]  -- image after bilateral filter applied (if enabled)
*  sobel[input-file-name]     -- image after sobel filter applied (this is input to Hough)
*  blob_candidate             -- b&w image with bolts found using blob detection represented as white and rest black.
*  hough_candidate 	      -- b&w image with bolts found using hough transform represented as white and rest black.
*  circles[input-file-name]   -- original image with circles found by hough added
*  final.jpg 		      -- final image.
*  FindBoltLocation.root      -- root histograms
*  bolts[input-file-name].txt -- text file containing [pmtid(-1 for now) pmtx pmty pmtr boltid boltx bolty] 
</pre>

To do list:
* Consider ellipse hough instead of circle to find the PMT rings of bolts?
[Paper for hough ellipse][http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.1.8792&rep=rep1&type=pdf]
* Cleanup / modularize code more?




