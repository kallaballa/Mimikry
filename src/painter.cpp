#include <limits>
#include "painter.hpp"


namespace mimikry {
  Painter::Painter(const Mat& oImg, const Mat& pImg):
      oImg_(oImg.clone()), pImg_(pImg.clone()), result_(oImg.clone()), fitness_(0), genome_(1) {
  }

  Painter::Painter(): oImg_(), pImg_(), result_(1,1,1,1), fitness_(0), genome_(0) {
  }

  void Painter::paint() {
    for(Chromosome& c : genome_) {
      if(!c.isActive() && genome_.size() != 1)
        continue;

      Mat kernel = c.makeKernel();
      Point anchor( -1, 1);
      double delta = 0;
      int ddepth = -1;
      int kernel_size = kernel.cols;
      filter2D(result_.clone(), result_, ddepth , kernel, anchor, delta, BORDER_DEFAULT );
    }

    for(int row = 0; row < result_.rows; ++row) {
        uchar* pr = result_.ptr(row);
        uchar* pp = pImg_.ptr(row);

        for(int col = 0; col < result_.cols*3; ++col) {
          fitness_ += 1.0 - (abs((*pp++) - (*pr++)) / 255.0);
        }

    }
    fitness_ /= (pImg_.cols*pImg_.rows*3);

    bool uniform = true;
    bool accumulate = false;
    cv::Mat a1_hist, a2_hist;
    int dims = 3;
    const int sizes[] = {256,256,256};
    const int channels[] = {0,1,2};
    float rRange[] = {0,256};
    float gRange[] = {0,256};
    float bRange[] = {0,256};
    const float* ranges[] = {rRange,gRange,bRange};

    cv::calcHist(&result_, 1, channels, cv::Mat(), a1_hist, dims, sizes, ranges, uniform, accumulate );
    cv::calcHist(&pImg_, 1, channels, cv::Mat(), a2_hist, dims, sizes, ranges, uniform, accumulate );

    double diff = cv::compareHist(a1_hist, a2_hist, CV_COMP_INTERSECT);
    fitness_ += diff;

    std::cerr << "fitness: " << fitness_  << std::endl;
  }
} /* namespace mimikry */

