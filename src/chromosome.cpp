#include "chromosome.hpp"

namespace mimikry {

  double fRand2(double fMin, double fMax) {
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
  }

  void Chromosome::init(size_t maxKernelSize) {
    CHECK(maxKernelSize >= 3);

    this->maxKernelSize_ = maxKernelSize;

    /*
     * size is pow(maxKernelSize, 2) + 3 because
     * we have the active gene, the kernel size gene,
     * the amplify gene and the values of the kernel
     */
    for(size_t i = 0; i < pow(maxKernelSize, 2) + 5; ++i) {
      this->push_back(fRand2(-1,1));
    }
  }

  bool Chromosome::isActive() const {
    return (*this)[0] > 0;
  }

  size_t Chromosome::getKernelSize() const {
    CHECK(maxKernelSize_ != 0);

    size_t kernelSize = 3 + ((1 +(*this)[1]) / 2) * (maxKernelSize_ - 3);

    if(kernelSize < 3)
      kernelSize = 3;
    else if(kernelSize > 20)
      kernelSize = 20;

    return kernelSize;
  }

  Operation Chromosome::getOperation() const {
    int op = round((*this)[4]);
    if(op <= -1)
      return SUBSTRACT;
    else if(op == 0)
      return PASS;
    else if(op >= 1)
      return ADD;

    CHECK(false);
    return PASS;
  }

  Mat Chromosome::makeKernel() const {
    CHECK(!this->empty());
    size_t kernelSize = getKernelSize();
    //std::cerr << "ks: " << kernelSize << std::endl;
    Mat kernel(kernelSize, kernelSize, CV_64F);
    for(size_t x = 0; x < kernelSize; ++x) {
      for(size_t y = 0; y < kernelSize; ++y) {
        kernel.at<double>(x, y) = (*this)[2] * ((*this)[5 + (x * kernelSize) + y]) * 16;
        //kernel.at<double>(x, y) = ((*this)[3 + (x * kernelSize) + y]);
      }
    }
    float sigma = (*this)[3] * 4;
    if(sigma != 0) {
      Mat blurred;
      cv::GaussianBlur(kernel, blurred, cv::Size(0, 0), fabs(sigma));

      if(sigma < 0) {
        cv::addWeighted(kernel, 1.5, blurred, -0.5, 0, blurred);
      }

      kernel = blurred;
    }
    return kernel;
  }
}
