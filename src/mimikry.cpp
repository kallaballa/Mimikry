#include "error.hpp"
#include "painter.hpp"
#include "genetic.hpp"
#include <iostream>
#include <string>
#include <vector>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>

namespace mimikry {
  using std::string;
  using std::vector;
  typedef unsigned char sample_t;
  using namespace cv;

  void run(const string& original, const string& processed, size_t population_size, size_t iterations) {
    Mat oImg = imread( original.c_str());
    Mat pImg = imread( processed.c_str());
    vector<Painter> population;

    std::cerr << "create population" << std::endl;
    for(size_t i = 0; i < population_size; ++i) {
      std::cerr << "i: " << i << std::endl;
      population.push_back(Painter(oImg, pImg));
    }

    GeneticLayout gl = make_default_genetic_layout();
    GeneticPool pool(gl);
    std::cerr << "run" << std::endl;

    for(size_t i = 0; i < iterations; ++i) {
      std::cerr << "i: " << i << std::endl;
      for(size_t j = 0; j < population.size(); ++j) {
        Painter& p = population[j];
        p.paint();

        imwrite(("result/" + std::to_string(i) + "_" + std::to_string(j) + ".png"),p.result_);
      }

      population = pool.epoch(population);
    }

    for(size_t j = 0; j < population.size(); ++j) {
 //     population[j].result_.save_png(("result/" + std::to_string(j) + ".png").c_str());
    }
  }
} /* namespace mimikry */


inline void default_error_delegate(const std::string& msg) {
  std::cerr << "### Error: " << msg << std::endl;
  mimikry::print_stacktrace();
  exit(1);
}

int main(int argc, char** argv) {
  srand (time(NULL));
  if(argc != 3) {
    std::cerr << "Usage: mimikry <original image> <processed image>" << std::endl;
    exit(1);
  }

  mimikry::ErrorHandler::init(default_error_delegate);
  mimikry::run(std::string(argv[1]), std::string(argv[2]), 100, 10000000);
  return 0;
}
