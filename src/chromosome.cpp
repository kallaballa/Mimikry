#include "chromosome.hpp"

namespace mimikry {
  double fRand2(double fMin, double fMax) {
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
  }

  void Chromosome::init(size_t size) {
    CHECK(size >= 10);
    CHECK(size % 2 == 0);

    for(size_t i = 0; i < size; ++i) {
      this->push_back(fRand2(-1,1));
    }
  }

  void Chromosome::check() {
    for(double& v : (*this)) {
      CHECK(v >= -1 && v <= 1);
    }
  }

  bool Chromosome::isActive() {
    return ((1 +(*this)[0]) / 2) * sqrt(this->size()) >= 3;
  }

  size_t Chromosome::getKernelSize() {
    size_t kernelSize = ((1 +(*this)[0]) / 2) * sqrt(this->size());
    std::cerr << "ks: " << kernelSize << std::endl;
    if(kernelSize < 3)
      kernelSize = 3;
    else if(kernelSize > 10)
      kernelSize = 10;

    return kernelSize;
  }

  Mat Chromosome::makeKernel() {
    CHECK(!this->empty());
    size_t kernelSize = getKernelSize();

    Mat kernel(kernelSize, kernelSize, CV_32F);
    for(int x = 0; x < kernelSize; ++x) {
      for(int y = 0; y < kernelSize; ++y) {
        kernel.at<float>(x, y) = ((*this)[1 + (x * kernelSize) + y] * 16);
      }
    }
    return kernel;
  }
}
