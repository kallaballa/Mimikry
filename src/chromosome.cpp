#include "chromosome.hpp"

namespace mimikry {
  double fRand2(double fMin, double fMax) {
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
  }

  void Chromosome::init(size_t maxKernelSize) {
    CHECK(maxKernelSize >= 3);
    CHECK(maxKernelSize >= 11);
    CHECK(maxKernelSize % 2 == 1);

    this->maxKernelSize_ = maxKernelSize;

    /*
     * size is pow(maxKernelSize, 2) + 2 because
     * we have the active gene, the kernel size gene
     * and the values of the kernel
     */
    for(size_t i = 0; i < pow(maxKernelSize, 2) + 2; ++i) {
      this->push_back(fRand2(-1,1));
    }
  }

  bool Chromosome::isActive() {
    return (*this)[0] > 0;
  }

  size_t Chromosome::getKernelSize() {
    CHECK(maxKernelSize_ != 0);

    size_t kernelSize = 3 + ((1 +(*this)[1]) / 2) * (maxKernelSize_ - 3);

    if(kernelSize < 3)
      kernelSize = 3;
    else if(kernelSize > 10)
      kernelSize = 10;

    return kernelSize;
  }

  Mat Chromosome::makeKernel() {
    CHECK(!this->empty());
    size_t kernelSize = getKernelSize();
    //std::cerr << "ks: " << kernelSize << std::endl;
    Mat kernel(kernelSize, kernelSize, CV_32F);
    for(size_t x = 0; x < kernelSize; ++x) {
      for(size_t y = 0; y < kernelSize; ++y) {
        //FIXME is 16 a good value?
        kernel.at<float>(x, y) = ((*this)[2 + (x * kernelSize) + y]);
      }
    }
    return kernel;
  }
}
