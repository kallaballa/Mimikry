#include <limits>
#include <opencv2/highgui/highgui.hpp>

#include "painter.hpp"

namespace mimikry {

int count_diff_pixels(cv::Mat in1, cv::Mat in2) {
    cv::Mat diff;
    cv::compare(in1, in2, diff, cv::CMP_NE);
    return cv::countNonZero(diff);
}

Mat make_histogram(const Mat& one) {
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

Mat make_histogram_gray(const Mat& oneRGB) {
  Mat one;
  cvtColor(oneRGB, one, CV_RGB2GRAY);

  bool uniform = true;
  bool accumulate = false;
  cv::Mat a1_hist, a2_hist;
  int dims = 1;
  const int sizes[] = { 255 };
  const int channels[] = { 0 };
  float gRange[] = { 0, 255 };
  const float* ranges[] = { gRange };

  cv::calcHist(&one, 1, channels, cv::Mat(), a1_hist, dims, sizes, ranges, uniform, accumulate);
  //imshow("Histogray", a1_hist);
  //waitKey(0);
  return a1_hist;
}

double compare_historgram(const Mat& one, const Mat& two) {
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

double compare_historgram_gray(const Mat& one, const Mat& two) {
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

double compare_float(const Mat& one, const Mat& two) {
  Scalar oMean;
  Scalar tMean;
  Scalar oSdv;
  Scalar tSdv;

  meanStdDev(one, oMean, oSdv);
  meanStdDev(two, tMean, tSdv);

  double mean = fabs(oMean.val[0] - tMean.val[0]);
  double sdv = fabs(oSdv.val[0] - tSdv.val[0]);

  return 1.0 - ((mean + sdv) / 2.0);
}

double compare_float_binarize(const Mat& one, const Mat& two) {
  Mat oNorm = one > 0.5;
  Mat tNorm = two > 0.5;

  Scalar oMean;
  Scalar tMean;
  Scalar oSdv;
  Scalar tSdv;

  double count = count_diff_pixels(oNorm, tNorm);
  return 1.0 - (count / (one.rows * one.cols));
}

Mat make_dft(const Mat& IRGB) {
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
/*
  int lowThreshold = 100;
  int ratio = 3;
  int kernel_size = 3;

  Mat detected_edges;
  Mat gray;
  Mat sharpened;
  magI.convertTo(gray, CV_8U);

  cv::GaussianBlur(gray, sharpened, cv::Size(0, 0), 3);
  cv::addWeighted(gray, 1.5, sharpened, -0.5, 0, sharpened);
  sharpened.convertTo(magI, CV_32F);*/

  return magI;
}

Painter::Painter(const Mat& oImg, const Mat& pImg) :
    oImg_(oImg),
    pImg_(pImg),
    pImgDFT_(make_dft(pImg)),
    pImgHist_(),
    result_(oImg.clone()),
    resultDFT_(),
    fitness_(0),
    pixelError_(0),
    histError_(0),
    fftError_(0),
    genome_(10) {
}

Painter::Painter() :
    oImg_(),
    pImg_(),
    pImgDFT_(),
    pImgHist_(),
    result_(1, 1, 1, 1),
    resultDFT_(),
    fitness_(0),
    pixelError_(0),
    histError_(0),
    fftError_(0),
    genome_(0) {
}

void Painter::paint() {
  size_t dominant = 0;
  size_t numActive = genome_.countActiveChromosomes();

  if (numActive == 0)
    dominant = genome_.findDominantChromosome();

  result_ = oImg_.clone();

  for (size_t i = 0; i < genome_.size(); ++i) {
    Chromosome& c = genome_[i];

    if (!c.isActive() && !(numActive == 0 && i == dominant))
      continue;

    Mat kernel = c.makeKernel();
    Point anchor(-1, -1);
    double delta = 0;
    int ddepth = -1;


    if(c.getOperation() == PASS) {
      filter2D(result_.clone(), result_, ddepth, kernel, anchor, delta, BORDER_DEFAULT);
    } else {
      Mat after;
      filter2D(oImg_.clone(), after, ddepth, kernel, anchor, delta, BORDER_DEFAULT);

      switch (c.getOperation()) {
      case PASS:
        CHECK(false);
        break;
      case SUBSTRACT:
        result_ = result_ - after;
        break;
      case ADD:
        result_ = result_ + after;
        break;
      default:
        CHECK(false);
        break;
      }
    }
  }

  double pixError = 0;


  Mat rNorm;
  Mat pNorm;
  normalize(result_, rNorm, 0, 255, CV_MINMAX);
  normalize(pImg_, pNorm, 0, 255, CV_MINMAX);

  Scalar rMean;
  Scalar pMean;
  Scalar rSdv;
  Scalar pSdv;

  meanStdDev(rNorm, rMean, rSdv);
  meanStdDev(pNorm, pMean, pSdv);

  double mean = (fabs(rMean.val[0] - pMean.val[0]) / 255.0
      + fabs(rMean.val[1] - pMean.val[1]) / 255.0
      + fabs(rMean.val[2] - pMean.val[2]) / 255.0)
          / 3.0;

  double sdv = (fabs(rSdv.val[0] - pSdv.val[0]) / 255.0
      + fabs(rSdv.val[1] - pSdv.val[1]) / 255.0
      + fabs(rSdv.val[2] - pSdv.val[2]) / 255.0)
          / 3.0;

  ///std::cerr << "mean: " << mean << "sdv: " << sdv << std::endl;

  pixError = 1-0 - ((mean + sdv) / 2.0);
  /*for (int row = 0; row < rNorm.rows; ++row) {
    uchar* pr = rNorm.ptr(row);
    uchar* pp = pNorm.ptr(row);
    for (int col = 0; col < rNorm.cols * 3; ++col) {
      pixError += 1.0 - (abs((*pp++) - (*pr++)) / 255.0);
    }
  }*/

  pixelError_ = pixError ;// / (result_.rows * result_.cols * 3);
 // histError_ = cv::compareHist(makeHistogramGray(result_), pImgHist_, CV_COMP_INTERSECT) / (result_.rows * result_.cols) ;
  resultDFT_ = make_dft(result_);
  fftError_ = compare_float(resultDFT_, pImgDFT_);
  histError_ = compare_float_binarize(resultDFT_, pImgDFT_);
  //(((double)(pImg_.rows * pImg_.cols) - count_diff_pixels(makeFFT(result_), pImgDFT_)) / ((double)pImg_.rows * pImg_.cols));
  fitness_ =  (histError_ * pixelError_ * fftError_);
  fitness_ = (fitness_ / pow(std::exp(-fitness_),2)) / (1 / pow(std::exp(-1),2));
}
} /* namespace mimikry */

