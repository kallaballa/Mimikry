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
      c.init(122);
      this->push_back(c);
    }
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
