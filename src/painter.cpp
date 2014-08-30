#include <limits>
#include <opencv2/highgui/highgui.hpp>

#include "painter.hpp"

namespace mimikry {

Mat makeHistogram(const Mat& one) {
  bool uniform = true;
  bool accumulate = false;
  cv::Mat a1_hist;
  int dims = 3;
  const int sizes[] = { 256, 256, 256 };
  const int channels[] = { 0, 1, 2 };
  float rRange[] = { 0, 256 };
  float gRange[] = { 0, 256 };
  float bRange[] = { 0, 256 };
  const float* ranges[] = { rRange, gRange, bRange };

  cv::calcHist(&one, 1, channels, cv::Mat(), a1_hist, dims, sizes, ranges, uniform, accumulate);

  return a1_hist;
}

double compareHistorgram(const Mat& one, const Mat& two) {
  bool uniform = true;
  bool accumulate = false;
  cv::Mat a1_hist, a2_hist;
  int dims = 3;
  const int sizes[] = { 256, 256, 256 };
  const int channels[] = { 0, 1, 2 };
  float rRange[] = { 0, 256 };
  float gRange[] = { 0, 256 };
  float bRange[] = { 0, 256 };
  const float* ranges[] = { rRange, gRange, bRange };

  cv::calcHist(&one, 1, channels, cv::Mat(), a1_hist, dims, sizes, ranges, uniform, accumulate);
  cv::calcHist(&two, 1, channels, cv::Mat(), a2_hist, dims, sizes, ranges, uniform, accumulate);

  return cv::compareHist(a1_hist, a2_hist, CV_COMP_INTERSECT) / 67200;
}

double compareHistorgramGray(const Mat& one, const Mat& two) {
  bool uniform = true;
  bool accumulate = false;
  cv::Mat a1_hist, a2_hist;
  int dims = 1;
  const int sizes[] = { 1 };
  const int channels[] = { 0 };
  float gRange[] = { 0, 1 };
  const float* ranges[] = { gRange };

  cv::calcHist(&one, 1, channels, cv::Mat(), a1_hist, dims, sizes, ranges, uniform, accumulate);
  cv::calcHist(&two, 1, channels, cv::Mat(), a2_hist, dims, sizes, ranges, uniform, accumulate);

  return cv::compareHist(a1_hist, a2_hist, CV_COMP_INTERSECT);
}

double compareFloat(const Mat& one, const Mat& two) {
  double result = 0;
  for (int row = 0; row < one.rows; ++row) {
    for (int col = 0; col < one.cols; ++col) {
      result += 1.0 - fabs(one.at<float>(row, col) - two.at<float>(row, col));
    }
  }
  return (result / (one.rows * one.cols));
}

Mat makeFFT(const Mat& IRGB) {
  Mat I;
  cvtColor(IRGB, I, CV_RGB2GRAY);
  Mat padded;                            //expand input image to optimal size
  int m = getOptimalDFTSize(I.rows);
  int n = getOptimalDFTSize(I.cols); // on the border add zero values
  copyMakeBorder(I, padded, 0, m - I.rows, 0, n - I.cols, BORDER_CONSTANT, Scalar::all(0));

  Mat planes[] = { Mat_<float>(padded), Mat::zeros(padded.size(), CV_32F) };
  Mat complexI;
  merge(planes, 2, complexI);         // Add to the expanded another plane with zeros

  dft(complexI, complexI);            // this way the result may fit in the source matrix

  // compute the magnitude and switch to logarithmic scale
  // => log(1 + sqrt(Re(DFT(I))^2 + Im(DFT(I))^2))
  split(complexI, planes);                   // planes[0] = Re(DFT(I), planes[1] = Im(DFT(I))
  magnitude(planes[0], planes[1], planes[0]);                   // planes[0] = magnitude
  Mat magI = planes[0];

  magI += Scalar::all(1);                    // switch to logarithmic scale
  log(magI, magI);

  // crop the spectrum, if it has an odd number of rows or columns
  magI = magI(Rect(0, 0, magI.cols & -2, magI.rows & -2));

  // rearrange the quadrants of Fourier image  so that the origin is at the image center
  int cx = magI.cols / 2;
  int cy = magI.rows / 2;

  Mat q0(magI, Rect(0, 0, cx, cy));   // Top-Left - Create a ROI per quadrant
  Mat q1(magI, Rect(cx, 0, cx, cy));  // Top-Right
  Mat q2(magI, Rect(0, cy, cx, cy));  // Bottom-Left
  Mat q3(magI, Rect(cx, cy, cx, cy)); // Bottom-Right

  Mat tmp;                           // swap quadrants (Top-Left with Bottom-Right)
  q0.copyTo(tmp);
  q3.copyTo(q0);
  tmp.copyTo(q3);

  q1.copyTo(tmp);                    // swap quadrant (Top-Right with Bottom-Left)
  q2.copyTo(q1);
  tmp.copyTo(q2);

  normalize(magI, magI, 0, 1, CV_MINMAX);
  //imshow( "Display window", magI);
  //waitKey(0);
  return magI;
}

Painter::Painter(const Mat& oImg, const Mat& pImg) :
    oImg_(oImg.clone()), pImg_(pImg.clone()), pImgDFT_(makeFFT(pImg_.clone())), pImgHist_(), result_(oImg.clone()), fitness_(0), genome_(10) {
}

Painter::Painter() :
    oImg_(), pImg_(), pImgDFT_(), pImgHist_(), result_(1, 1, 1, 1), fitness_(0), genome_(0) {
}

void Painter::paint() {
  size_t dominant = 0;
  size_t numActive = genome_.countActiveChromosomes();

  if(numActive == 0)
    dominant = genome_.findDominantChromosome();

  for (size_t i = 0; i < genome_.size(); ++i) {
    Chromosome& c = genome_[i];

    if (!c.isActive() && !(numActive == 0 && i == dominant))
      continue;

    Mat kernel = c.makeKernel();
    Point anchor(-1, 1);
    double delta = 0;
    int ddepth = -1;
    filter2D(result_.clone(), result_, ddepth, kernel, anchor, delta, BORDER_DEFAULT);
  }
  /*
   for(int row = 0; row < result_.rows; ++row) {
   uchar* pr = result_.ptr(row);
   uchar* pp = pImg_.ptr(row);

   for(int col = 0; col < result_.cols*3; ++col) {
   fitness_ += 1.0 - (abs((*pp++) - (*pr++)) / 255.0);
   }

   }
   fitness_ /= (pImg_.cols*pImg_.rows*3);
   */

  fitness_ += (cv::compareHist(makeHistogram(result_), makeHistogram(pImg_), CV_COMP_INTERSECT) * compareFloat(makeFFT(result_), pImgDFT_));
}
} /* namespace mimikry */

