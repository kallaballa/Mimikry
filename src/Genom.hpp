/*
 * Genom.hpp
 *
 *  Created on: Aug 28, 2014
 *      Author: elchaschab
 */

#ifndef GENOM_HPP_
#define GENOM_HPP_

#include <vector>
#include "chromosome.hpp"

namespace mimikry {

using std::vector;

class Genome : public vector<Chromosome> {
public:

  Genome(size_t size) : vector<Chromosome>() {
    for(size_t i = 0; i < size; ++i) {
      Chromosome c;
      c.init(20);
      this->push_back(c);
    }
  }

  size_t countActiveChromosomes() {
    size_t cnt = 0;
    for(Chromosome& c : (*this)) {
      cnt += c.isActive() ? 1 : 0;
    }
    return cnt;
  }

  size_t getTotalKernelSize() {
    size_t total = 0;
    for(Chromosome& c : (*this)) {
      total += c.isActive() ? c.getKernelSize()  : 0;
    }
    return total;
  }

  size_t findDominantChromosome() {
    double max = -1;
    size_t dominant = 0;
    size_t i = 0;

    for(Chromosome& c : (*this)) {
      if(c[0] > max) {
        max = c[0];
        dominant = i;
      }
      ++i;
    }
    return dominant;
  }

  virtual ~Genome() {
  }

  void check() {
    for(Chromosome& c : (*this)) {
      c.check();
    }
  }
};


} /* namespace mimikry */

#endif /* GENOM_HPP_ */
