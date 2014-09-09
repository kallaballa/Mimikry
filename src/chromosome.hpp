#ifndef CHROMOSOME_H_
#define CHROMOSOME_H_

#include <vector>
#include <opencv2/imgproc/imgproc.hpp>
#include "error.hpp"

namespace mimikry {

using namespace cv;
using std::vector;
enum Operation {
  PASS,
  ADD,
  SUBSTRACT,
};

class Chromosome : public vector<double> {
  size_t maxKernelSize_ = 0;
public:
  Chromosome() {
  }

  virtual ~Chromosome() {
  }

  void init(size_t size);
  void check();
  bool isActive() const ;
  size_t getKernelOffset() const;
  size_t getKernelSize() const;
  double getAmplify() const;
  double getSigma() const;
  size_t getPrefilterIndex() const;
  Operation getOperation() const;
  Mat makeKernel() const;
};
} /* namespace mimikry */

#endif /* CHROMOSOME_H_ */
