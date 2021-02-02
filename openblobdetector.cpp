#include "openblobdetector.hpp"

#include <opencv2/imgproc.hpp>
#include <opencv2/features2d.hpp>
#include <opencv2/core/core.hpp>

using namespace cv;

OpenBlobDetector::Params::Params()
{
  thresholdStep = 10;
  minThreshold = 50;
  maxThreshold = 220;
  minRepeatability = 2;
  minDistBetweenBlobs = 10;

  filterByColor = true;
  blobColor = 0;

  filterByArea = true;
  minArea = 25;
  maxArea = 5000;

  filterByCircularity = false;
  minCircularity = 0.8f;
  maxCircularity = std::numeric_limits<float>::max();

  filterByInertia = true;
  //minInertiaRatio = 0.6;
  minInertiaRatio = 0.1f;
  maxInertiaRatio = std::numeric_limits<float>::max();

  filterByConvexity = true;
  //minConvexity = 0.8;
  minConvexity = 0.95f;
  maxConvexity = std::numeric_limits<float>::max();
}

void OpenBlobDetector::Params::read(const cv::FileNode& fn )
{
  thresholdStep = fn["thresholdStep"];
  minThreshold = fn["minThreshold"];
  maxThreshold = fn["maxThreshold"];

  minRepeatability = (size_t)(int)fn["minRepeatability"];
  minDistBetweenBlobs = fn["minDistBetweenBlobs"];

  filterByColor = (int)fn["filterByColor"] != 0 ? true : false;
  blobColor = (uchar)(int)fn["blobColor"];

  filterByArea = (int)fn["filterByArea"] != 0 ? true : false;
  minArea = fn["minArea"];
  maxArea = fn["maxArea"];

  filterByCircularity = (int)fn["filterByCircularity"] != 0 ? true : false;
  minCircularity = fn["minCircularity"];
  maxCircularity = fn["maxCircularity"];

  filterByInertia = (int)fn["filterByInertia"] != 0 ? true : false;
  minInertiaRatio = fn["minInertiaRatio"];
  maxInertiaRatio = fn["maxInertiaRatio"];

  filterByConvexity = (int)fn["filterByConvexity"] != 0 ? true : false;
  minConvexity = fn["minConvexity"];
  maxConvexity = fn["maxConvexity"];
}

void OpenBlobDetector::Params::write(cv::FileStorage& fs) const
{
  fs << "thresholdStep" << thresholdStep;
  fs << "minThreshold" << minThreshold;
  fs << "maxThreshold" << maxThreshold;

  fs << "minRepeatability" << (int)minRepeatability;
  fs << "minDistBetweenBlobs" << minDistBetweenBlobs;

  fs << "filterByColor" << (int)filterByColor;
  fs << "blobColor" << (int)blobColor;

  fs << "filterByArea" << (int)filterByArea;
  fs << "minArea" << minArea;
  fs << "maxArea" << maxArea;

  fs << "filterByCircularity" << (int)filterByCircularity;
  fs << "minCircularity" << minCircularity;
  fs << "maxCircularity" << maxCircularity;

  fs << "filterByInertia" << (int)filterByInertia;
  fs << "minInertiaRatio" << minInertiaRatio;
  fs << "maxInertiaRatio" << maxInertiaRatio;

  fs << "filterByConvexity" << (int)filterByConvexity;
  fs << "minConvexity" << minConvexity;
  fs << "maxConvexity" << maxConvexity;
}


OpenBlobDetectorImpl::OpenBlobDetectorImpl(const OpenBlobDetector::Params &parameters) :
params(parameters)
{
}

void OpenBlobDetectorImpl::read( const cv::FileNode& fn )
{
  params.read(fn);
}

void OpenBlobDetectorImpl::write( cv::FileStorage& fs ) const
{
  writeFormat(fs);
  params.write(fs);
}

void OpenBlobDetectorImpl::findBlobs(InputArray _image, InputArray _binaryImage, std::vector<Center> &centers) const
{

  //CV_INSTRUMENT_REGION();

    Mat image = _image.getMat(), binaryImage = _binaryImage.getMat();
    CV_UNUSED(image);
    centers.clear();

    std::vector < std::vector<Point> > contours;
    findContours(binaryImage, contours, RETR_LIST, CHAIN_APPROX_NONE);

    for (size_t contourIdx = 0; contourIdx < contours.size(); contourIdx++)
    {
      Center center;

        center.confidence = 1;
        Moments moms = moments(contours[contourIdx]);
        if (params.filterByArea)
        {
            double area = moms.m00;
	    center.area = area;
            if (area < params.minArea || area >= params.maxArea)
                continue;
        }

        if (params.filterByCircularity)
        {
            double area = moms.m00;
            double perimeter = arcLength(contours[contourIdx], true);
            double ratio = 4 * CV_PI * area / (perimeter * perimeter);
	    center.circularity = ratio;
            if (ratio < params.minCircularity || ratio >= params.maxCircularity)
                continue;
        }

        if (params.filterByInertia)
        {
            double denominator = std::sqrt(std::pow(2 * moms.mu11, 2) + std::pow(moms.mu20 - moms.mu02, 2));
            const double eps = 1e-2;
            double ratio;
            if (denominator > eps)
            {
                double cosmin = (moms.mu20 - moms.mu02) / denominator;
                double sinmin = 2 * moms.mu11 / denominator;
                double cosmax = -cosmin;
                double sinmax = -sinmin;

                double imin = 0.5 * (moms.mu20 + moms.mu02) - 0.5 * (moms.mu20 - moms.mu02) * cosmin - moms.mu11 * sinmin;
                double imax = 0.5 * (moms.mu20 + moms.mu02) - 0.5 * (moms.mu20 - moms.mu02) * cosmax - moms.mu11 * sinmax;
                ratio = imin / imax;
            }
            else
            {
                ratio = 1;
            }
	    center.inertia = ratio;
            if (ratio < params.minInertiaRatio || ratio >= params.maxInertiaRatio)
                continue;

            center.confidence = ratio * ratio;
        }

        if (params.filterByConvexity)
        {
            std::vector < Point > hull;
            convexHull(contours[contourIdx], hull);
            double area = moms.m00;
            double hullArea = contourArea(hull);
            if (fabs(hullArea) < DBL_EPSILON)
                continue;
            double ratio = area / hullArea;
	    center.convexity = ratio;
            if (ratio < params.minConvexity || ratio >= params.maxConvexity)
                continue;
        }

        if(moms.m00 == 0.0)
            continue;
        center.location = Point2d(moms.m10 / moms.m00, moms.m01 / moms.m00);

        if (params.filterByColor)
        {
            if (binaryImage.at<uchar> (cvRound(center.location.y), cvRound(center.location.x)) != params.blobColor)
                continue;
        }

        //compute blob radius
        {
            std::vector<double> dists;
            for (size_t pointIdx = 0; pointIdx < contours[contourIdx].size(); pointIdx++)
            {
                Point2d pt = contours[contourIdx][pointIdx];
                dists.push_back(norm(center.location - pt));
            }
            std::sort(dists.begin(), dists.end());
            center.radius = (dists[(dists.size() - 1) / 2] + dists[dists.size() / 2]) / 2.;
        }

        centers.push_back(center);

    }

}

void OpenBlobDetectorImpl::detect(InputArray image, std::vector<cv::KeyPoint>& keypoints, InputArray mask)
{
  //CV_INSTRUMENT_REGION();

    keypoints.clear();
    blobinfo.clear();
    CV_Assert(params.minRepeatability != 0);
    Mat grayscaleImage;
    if (image.channels() == 3 || image.channels() == 4)
        cvtColor(image, grayscaleImage, COLOR_BGR2GRAY);
    else
        grayscaleImage = image.getMat();

    if (grayscaleImage.type() != CV_8UC1) {
        CV_Error(Error::StsUnsupportedFormat, "Blob detector only supports 8-bit images!");
    }

    std::vector < std::vector<Center> > centers;
    for (double thresh = params.minThreshold; thresh < params.maxThreshold; thresh += params.thresholdStep)
    {
        Mat binarizedImage;
        threshold(grayscaleImage, binarizedImage, thresh, 255, THRESH_BINARY);

        std::vector < Center > curCenters;
        findBlobs(grayscaleImage, binarizedImage, curCenters);
        std::vector < std::vector<Center> > newCenters;
        for (size_t i = 0; i < curCenters.size(); i++)
        {
            bool isNew = true;
            for (size_t j = 0; j < centers.size(); j++)
            {
                double dist = norm(centers[j][ centers[j].size() / 2 ].location - curCenters[i].location);
                isNew = dist >= params.minDistBetweenBlobs && dist >= centers[j][ centers[j].size() / 2 ].radius && dist >= curCenters[i].radius;
                if (!isNew)
                {
                    centers[j].push_back(curCenters[i]);

                    size_t k = centers[j].size() - 1;
                    while( k > 0 && curCenters[i].radius < centers[j][k-1].radius )
                    {
                        centers[j][k] = centers[j][k-1];
                        k--;
                    }
                    centers[j][k] = curCenters[i];

                    break;
                }
            }
            if (isNew)
                newCenters.push_back(std::vector<Center> (1, curCenters[i]));
        }
        std::copy(newCenters.begin(), newCenters.end(), std::back_inserter(centers));
    }

    for (size_t i = 0; i < centers.size(); i++)
    {
        if (centers[i].size() < params.minRepeatability)
            continue;
        Point2d sumPoint(0, 0);
        double normalizer = 0;
        for (size_t j = 0; j < centers[i].size(); j++)
        {
            sumPoint += centers[i][j].confidence * centers[i][j].location;
            normalizer += centers[i][j].confidence;
        }
        sumPoint *= (1. / normalizer);
        KeyPoint kpt(sumPoint, (float)(centers[i][centers[i].size() / 2].radius) * 2.0f);
        keypoints.push_back(kpt);
	blobinfo.push_back( centers[i][centers[i].size() / 2] );
    }

    if (!mask.empty())
    {
        KeyPointsFilter::runByPixelsMask(keypoints, mask.getMat());
    }
}

Ptr<OpenBlobDetector> OpenBlobDetector::create(const OpenBlobDetector::Params& params)
{
    return makePtr<OpenBlobDetectorImpl>(params);
}

String OpenBlobDetector::getDefaultName() const
{
    return (Feature2D::getDefaultName() + ".OpenBlobDetector");
}


