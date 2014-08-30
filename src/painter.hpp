#ifndef PAINTER_HPP_
#define PAINTER_HPP_

#include <opencv2/imgproc/imgproc.hpp>
#include <string>
#include <vector>
#include "Genom.hpp"
#include "error.hpp"

namespace mimikry {

using std::string;
using std::vector;
typedef unsigned char sample_t;

class Painter {
public:
  Mat oImg_;
  Mat pImg_;
  Mat pImgDFT_;
  Mat pImgHist_;
  Mat result_;
  double fitness_;
  Genome genome_;

  Painter(const Mat& oImg, const Mat& pImg);
  Painter();

  void paint();

  Painter makeChild() const {
    Painter child(oImg_, pImg_);

    return child;
  }

  Painter clone() const {
    Painter child(oImg_, pImg_);

    //copy chromosomes
    for (size_t i = 0; i < genome_.size(); ++i) {
      for (size_t j = 0; j < genome_[i].size(); ++j) {
        child.genome_[i][j] = genome_[i][j] ;
      }
    }
    return child;
  }


  bool operator==(const Painter& other) const {
    for (size_t i = 0; i < genome_.size(); ++i) {
      for (size_t j = 0; j < genome_[i].size(); ++j) {
        if(other.genome_[i][j] != genome_[i][j])
          return false;
      }
    }
    return true;
  }

  bool operator!=(const Painter& other) const {
    return !this->operator ==(other);
  }

  bool operator<(const Painter& other) const {
    return (this->fitness_ < other.fitness_);
  }

};

typedef vector<Painter> Population;

} /* namespace mimikry */

#endif /* PAINTER_HPP_ */
