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
*  Config_files		      -- Directory contaning config files for different images.
*houghellipse[file_number]    -- final ellipse found after pruning.
*houghellipse_before[file_number]      -- ellipses found before pruning.
</pre>

<table style ="width:100%;">
<tr>
<td>
<img src="./Config_files/239/houghellipse239.jpg" height="25%" width="100%"> 
</td>
<td>
<img src="./Config_files/379/houghellipse379.jpg" height="25%" width="100%"> 
</td>
</tr>
<tr>
<td>
  Side Image(houghellipse239.jpg)
</td>
<td>
  Corner Image(houghellipse379.jpg)
</td>
</tr>
</table>

To do list:
* Get it to work better for corner images
  * Find some way to filter blobs that are inside of the ellipses of bolts
  * Maybe could do a first hough transform on one of the layers (colors) or filetered image to find ellipses from dynode or edge of PMT    features that could help get rid of blobs inside the PMT ellipses
  * Could try picking from candidate ellipses by looking at properties of contained blobs

* For corner images, look for feature of edge of tank?  Could be used to separate PMTs on wall from those on bottom or top

* Start reading in some of the drone distance to wall and pointing data to get range of possible PMT ellipse shapes
    to tune the ellipse-hough parameters








