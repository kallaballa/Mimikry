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
     * size is pow(maxKernelSize, 2) + 6 because
     * we have the active gene, the kernel size gene,
     * the amplify gene and the values of the kernel
     */
    for(size_t i = 0; i < pow(maxKernelSize, 2) + getKernelOffset(); ++i) {
      this->push_back(fRand2(-1,1));
    }
  }

  size_t Chromosome::getKernelOffset() const {
    return 6;
  }

  bool Chromosome::isActive() const {
    return ((int)round((*this)[0])) % 2;
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

  double Chromosome::getAmplify() const {
    return (*this)[2];
  }

  double Chromosome::getSigma() const {
    return (*this)[3] * 4;
  }

  Operation Chromosome::getOperation() const {
    int op = abs(((int)round((*this)[4])) % 3);
    if(op == 0)
      return SUBSTRACT;
    else if(op == 1)
      return PASS;
    else if(op == 2)
      return ADD;

    CHECK(false);
    return PASS;
  }

  size_t Chromosome::getPrefilterIndex() const {
    return (abs(round((*this)[5] * 5))) % 4;
  }

  Mat Chromosome::makeKernel() const {
    CHECK(!this->empty());
    size_t kernelSize = getKernelSize();
    //std::cerr << "ks: " << kernelSize << std::endl;
    Mat kernel(kernelSize, kernelSize, CV_64F);
    for(size_t x = 0; x < kernelSize; ++x) {
      for(size_t y = 0; y < kernelSize; ++y) {
        kernel.at<double>(x, y) = getAmplify() * ((*this)[6 + (x * kernelSize) + y]) * 16;
        //kernel.at<double>(x, y) = ((*this)[3 + (x * kernelSize) + y]);
      }
    }
    float sigma = getSigma();
    size_t preFilter = getPrefilterIndex();
    Mat filtered;
    int s ;

    switch(preFilter) {
    case 0:
      //no prefilter
      break;
    case 1: //blur
      cv::GaussianBlur(kernel, filtered, cv::Size(0, 0), fabs(sigma));
      kernel = filtered;
      break;
    case 2: //unsharp mask
      cv::GaussianBlur(kernel, filtered, cv::Size(0, 0), fabs(sigma));
      cv::addWeighted(kernel, 1.5, filtered, -0.5, 0, filtered);
      kernel = filtered;
      break;
    case 3: //laplace
      s = 2 + abs(((int)round(sigma)));
      if((s % 2) == 0)
        ++s;
      Laplacian(kernel, filtered, -1, s);
      kernel = filtered;
      break;
    default:
      CHECK(false);
      break;
    }

    return kernel;
  }
}
