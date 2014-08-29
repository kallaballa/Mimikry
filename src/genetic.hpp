#ifndef CGENALG_H
#define CGENALG_H

#include "painter.hpp"
#include "Genom.hpp"

using namespace std;

namespace mimikry {

struct BattleFieldLayout;

struct GeneticLayout {
	double mutationRate;
	double crossoverRate;
	size_t crossoverIterations;
	double maxPertubation;
	size_t numElite_;
	size_t numEliteCopies_;
	bool usePerfDesc_;
};

GeneticLayout make_default_genetic_layout();

//-----------------------------------------------------------------------
//
//	the genetic algorithm class
//-----------------------------------------------------------------------
class GeneticPool {
private:
	bool initialized_ = false;
	void mutate(Genome& genome);
	Painter& pickSpecimen(Population& pop);
	std::pair<Painter, Painter> crossover(Painter &mum, Painter &dad, size_t iterations);
	void copyNBest(size_t n, const size_t numCopies, Population& in, Population& out);
	void calculateStatistics(Population& pop);
public:
	GeneticLayout layout_;
	GeneticPool(GeneticLayout params);
	GeneticPool();

	//this runs the GA for one generation.
	virtual Population epoch(Population& old_pop);
};
}

#endif

