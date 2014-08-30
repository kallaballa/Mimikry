#ifndef CHROMOSOME_H_
#define CHROMOSOME_H_

#include <vector>
#include <opencv2/imgproc/imgproc.hpp>
#include "error.hpp"

namespace mimikry {

using namespace cv;
using std::vector;

class Chromosome : public vector<double> {
public:
  Chromosome() {
  }

  virtual ~Chromosome() {
  }

  void init(size_t size);
  void check();
  bool isActive();
  size_t getKernelSize();
  Mat makeKernel();

  size_t maxKernelSize_ = 0;
};
} /* namespace mimikry */

#endif /* CHROMOSOME_H_ */
