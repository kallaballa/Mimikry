#include "genetic.hpp"
#include "error.hpp"
#include "painter.hpp"
#include <algorithm>
#include <iostream>
#include <set>

namespace mimikry {

int iRand(int x,int y) {
  return rand()%(y-x+1)+x;
}

double fRand(double fMin, double fMax) {
  double f = (double)rand() / RAND_MAX;
  return fMin + f * (fMax - fMin);
}

GeneticLayout make_default_genetic_layout() {
	return {
			0.2,// mutationRate
			0.7, // crossoverRate
			1,   // crossoverIterations
			0.3,   // maxPertubation
			4,   // numElite
			1,   // numEliteCopies
			false// usePerfDesc_
	};
}

using std::set;

GeneticPool::GeneticPool(GeneticLayout params) :
		layout_(params) {
	initialized_ = true;
}

GeneticPool::GeneticPool() :
		layout_() {
	initialized_ = false;
}

/*
 * mutates a chromosome by perturbing its weights by an amount not greater than params_.maxPerturbation_
 */
void GeneticPool::mutate(Genome& genome) {
	//traverse the chromosome and mutate each weight dependent
	//on the mutation rate

  for (Chromosome& c : genome) {
    for (double& w : c) {
      //do we perturb this weight?
      if (fRand(0, 1) < layout_.mutationRate) {
        //add or subtract a small value to the weight
        w += ((fRand(0, 1) - fRand(0, 1)) * layout_.maxPertubation);
      }
    }
  }
}

/*
 * Returns a Painter base on roulette wheel sampling
 */
Painter& GeneticPool::pickSpecimen(Population& pop) {
	double totalFitness = 0;
	for (size_t i = 0; i < pop.size(); ++i) {
	  CHECK(pop[i].fitness_ >= 0);
	  totalFitness += pop[i].fitness_;
	}

	//generate a random number between 0 & total fitness count
	double slice = (double) fRand(0, totalFitness);

	double fitnessSoFar = 0;
	for (size_t i = 0; i < pop.size(); ++i) {
		fitnessSoFar += pop[i].fitness_;

		//if the fitness so far > random number return the tank at this point
		if (fitnessSoFar >= slice) {
			return pop[i];
		}
	}

	CHECK(false);
	return pop[0]; 	//surpress warning
}

/*
 * With a chance defined by params_.crossoverRate_ perform a crossover of brain_.weights()
 */
std::pair<Painter, Painter> GeneticPool::crossover(Painter &mum, Painter &dad, size_t iterations) {
	Painter baby1 = mum.makeChild();
	Painter baby2 = mum.makeChild();

	for(size_t c = 0; c < mum.genome_.size(); ++c) {
		Chromosome& wMum = mum.genome_[c];
		Chromosome& wDad = dad.genome_[c];
		Chromosome& wBaby1 = baby1.genome_[c];
		Chromosome& wBaby2 = baby2.genome_[c];

		if ((fRand(0,1) > layout_.crossoverRate) || (mum == dad)) {
			for (size_t i = 0; i  < wMum.size(); ++i ) {
				wBaby1[i] = wMum[i];
				wBaby2[i] = wDad[i];
			}
		} else {
      size_t cp = iRand(0, wBaby1.size() - 1);

      for (size_t j = 0; j < cp; ++j) {
          wBaby1[j] = wMum[j];
          wBaby2[j] = wDad[j];
      }
      for (size_t j = cp; j < wBaby1.size(); ++j) {
          wBaby2[j] = wMum[j];
          wBaby1[j] = wDad[j];
      }
		}
	}
	return {baby1, baby2};
}

/*
 * copy numCopies copies of the n best specimen into the out population
 */
void GeneticPool::copyNBest(size_t n, const size_t numCopies,
		Population& in, Population& out) {
	//add the required amount of copies of the n most fittest to the supplied population
	while (n--) {
		for (size_t i = 0; i < numCopies; ++i) {
			Painter& t = in[(in.size() - 1) - n];
			Painter clone = t.clone();
			out.push_back(clone);
		}
	}
}

/*
 * generates the statistics for the given population
 */
void GeneticPool::calculateStatistics(Population& pop) {

}

/*
 * Use the genetic algorithm to construct a new population from the old
 */
Population GeneticPool::epoch(Population& old_pop) {
	if(!initialized_ || old_pop.size() == 1) {
		Population new_pop = old_pop;
		new_pop.clear();

		//sort the population (for scaling and elitism)
		sort(old_pop.begin(), old_pop.end());

		//calculate best, worst, average and total fitness
		calculateStatistics(old_pop);

		for(Painter& t : old_pop) {
			new_pop.push_back(t.makeChild());
		}

		return new_pop;
	}

	Population new_pop = old_pop;
	new_pop.clear();

	//sort the population (for scaling and elitism)
	sort(old_pop.begin(), old_pop.end());

	//calculate best, worst, average and total fitness
	calculateStatistics(old_pop);

  if (!(layout_.numEliteCopies_ * (layout_.numElite_ % 2))) {
    copyNBest(layout_.numElite_, layout_.numEliteCopies_, old_pop, new_pop);
  }

	//now we enter the GA loop
	//repeat until a new population is generated
	while (new_pop.size() < old_pop.size()) {
		//grab two chromosones
		Painter& mum = pickSpecimen(old_pop);
		Painter& dad = pickSpecimen(old_pop);

		//create some offspring via crossover
		std::pair<Painter, Painter> babies = crossover(mum, dad, layout_.crossoverIterations);

		//now we mutate
		mutate(babies.first.genome_);
		mutate(babies.second.genome_);

		//now copy into vecNewPop population
		new_pop.push_back(babies.first);
		if(new_pop.size() < old_pop.size())
			new_pop.push_back(babies.second);
	}
	CHECK(new_pop.size() == old_pop.size());
	return new_pop;
}
}
