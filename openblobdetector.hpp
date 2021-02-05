#ifndef OPENBLOBDETECTOR_HPP
#define OPENBLOBDETECTOR_HPP

#include <iterator>
#include <limits>
#include <opencv2/features2d.hpp>

namespace cv
{

  class OpenBlobDetector : public Feature2D {
  public:

    struct Center
    {
      Point2d location;
      double radius;
      double confidence;
      double circularity;
      double inertia;
      double convexity;
      double area;
      double intensity;
    };
    
    // filled once for each keypoint returned by detect method
    std::vector< Center > blobinfo;
    
    
    struct Params
    {
      Params();
      float thresholdStep;
      float minThreshold;
      float maxThreshold;
      size_t minRepeatability;
      float minDistBetweenBlobs;
      
      bool filterByColor;
      uchar blobColor;
      
      bool filterByArea;
      float minArea, maxArea;
      
      bool filterByCircularity;
      float minCircularity, maxCircularity;
      
      bool filterByInertia;
      float minInertiaRatio, maxInertiaRatio;
      
      bool filterByConvexity;
      float minConvexity, maxConvexity;
      
      void read( const FileNode& fn );
      void write( FileStorage& fs ) const;
    };
    
    static Ptr<OpenBlobDetector>
      create(const OpenBlobDetector::Params &parameters = OpenBlobDetector::Params());
    virtual String getDefaultName() const CV_OVERRIDE;
  };
  
  
  class OpenBlobDetectorImpl : public OpenBlobDetector
  {
    
  public:

    explicit OpenBlobDetectorImpl(const OpenBlobDetector::Params &parameters = OpenBlobDetector::Params());
    virtual void detect( InputArray image, std::vector<KeyPoint>& keypoints, InputArray mask=noArray() ) CV_OVERRIDE;
    
    virtual void read( const FileNode& fn ) CV_OVERRIDE;
    virtual void write( FileStorage& fs ) const CV_OVERRIDE;
    
    
  protected:
    
    
    virtual void findBlobs(InputArray image, InputArray binaryImage, std::vector<Center> &centers) const;
    
    Params params;
    
  };
}

#endif
