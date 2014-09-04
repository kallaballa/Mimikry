#include "error.hpp"
#include "painter.hpp"
#include "genetic.hpp"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <iomanip>
#include <sstream>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>

namespace mimikry {
  using std::string;
  using std::vector;
  using std::stringstream;
  typedef unsigned char sample_t;
  using namespace cv;

  struct GlobalStatistics {
    size_t streak_;
    double lastBestFitness_;
  };

  GlobalStatistics global_stats = {0,0.0};

  struct PopulationStatistics {
    double bestFitness_;
    double avgFitness_;
    double avgPixelError_;
    double avgHistError_;
    double avgKernelSize_;
  };

  void updateGlobalStats(PopulationStatistics popStats) {
    if(global_stats.lastBestFitness_ == popStats.bestFitness_)
      global_stats.streak_++;
    else
      global_stats.streak_ = 0;

    global_stats.lastBestFitness_ = popStats.bestFitness_;
  }

  PopulationStatistics calcPopulationStats(const Population& pop) {
    vector<PopulationStatistics> stats;
    double bf = 0;
    double totalFitness = 0;
    double totalPixErr = 0;
    double totalHistErr = 0;
    double totalKernelSize = 0;
    double activeChromos = 0;

    for(const Painter& p : pop) {
      bf=std::max(p.fitness_, bf);
      totalFitness += p.fitness_;
      totalPixErr += p.pixelError_;
      totalHistErr += p.histError_;

      for(const Chromosome& c : p.genome_) {
        if(c.isActive()) {
          totalKernelSize += c.getKernelSize();
          ++activeChromos;
        }
      }
    }

    size_t size = pop.size();
    return {bf, totalFitness / size, totalPixErr / size, totalHistErr / size, totalKernelSize / activeChromos };
  }

  Mat combine_kernels(Mat& one, Mat& two) {
    size_t combinedKernelSize = one.cols + two.cols - 1;
    size_t off = (combinedKernelSize - one.cols) / 2;
    Mat combined(combinedKernelSize,combinedKernelSize, CV_64F, double(0));

    for(size_t x = 0; x < one.rows; ++x) {
      for(size_t y = 0; y < one.cols; ++y) {
        combined.at<double>(x + off, y + off) = one.at<double>(x, y);
      }
    }

    Point anchor(-1, -1);
    double delta = 0;
    int ddepth = -1;
    filter2D(combined.clone(), combined, ddepth, two, anchor, delta, BORDER_DEFAULT);
    return combined;
  }

  Chromosome make_chromo(const Mat& m, const int size) {
    Chromosome c;
    c.init(size);
    for(size_t x = 0; x < m.rows; ++x) {
      for(size_t y = 0; y < m.cols; ++y) {
        c[5 + (x * m.cols) + y] = m.at<double>(x, y);
      }
    }

    //FIXME properly calculate the weight
    c[1] = -1.0;
    c[0] = 0.1;
    while(c.getKernelSize() < m.rows) {
      c[1] += 0.01;
    }

    CHECK(c.getKernelSize() == m.rows);

    return c;
  }

  void pack_genome(Painter& p) {
    vector<Mat> activeKernels;

    //collective active kernels
    for(Chromosome& c : p.genome_) {
      if(c.isActive()) {
        c[0] = 0.001; //FIXME
        activeKernels.push_back(c.makeKernel());
      }
      else
        c[0] = -0.001; //FIXME
    }

    return; //FIXME

    //nothing else to do if there is only one active/dominant chromosome
    if(p.genome_.countActiveChromosomes() <= 1) {
      std::cerr << "#### not enough active chromos ####" << std::endl;
      return;
    }

    //sort by kernel size
    sort(activeKernels.begin(), activeKernels.end(), [](const Mat& one, const Mat& two) {
      return one.cols > two.cols;
    });

    //record current fitness so we can check correctness of the resulting kernel
    double fitnessBefore = p.fitness_;
    size_t i = 1;

    std::cerr << "active before: " << activeKernels.size() << std::endl;

    while(i < activeKernels.size() && activeKernels.size() > 1) {
      CHECK(i > 0);
      Mat& one = activeKernels[i -1];
      Mat& two = activeKernels[i];

      std::cerr << "one:" << one.cols << std::endl;
      std::cerr << "two:" << two.cols << std::endl;
      std::cerr << "pow(one + two -1):" << pow(one.cols + two.cols - 1,2) << std::endl;
      std::cerr << "pfs" << p.genome_.front().size() << std::endl;

      //check if the resulting kernel would fit into a chromosome
      if(pow(one.cols + two.cols - 1,2) < p.genome_.front().size()) {
        Mat combined = combine_kernels(one,two);
        CHECK(combined.cols == one.cols + two.cols - 1);
        activeKernels[i - 1] = combined;
        activeKernels.erase((activeKernels.begin() + i));
        CHECK(activeKernels[i - 1].cols == combined.cols);
        i=1;
      } else {
        ++i;
      }
    }

    std::cerr << "active after: " << activeKernels.size() << std::endl;

    Genome newGenome(p.genome_.size() - activeKernels.size());
    for(Chromosome& c : newGenome) {
      c[0] = -0.1;
    }

    for(const Mat& k : activeKernels) {
      newGenome.push_back(make_chromo(k, 20));
    }
    CHECK(p.genome_.size() == newGenome.size());

    imwrite(("result/before.png"),p.result_);
    p.genome_ = newGenome;
    p.paint();
    imwrite(("result/after.png"),p.result_);

    std::cerr << fitnessBefore << "/"<<  p.fitness_ << " " << p.genome_.getTotalKernelSize() << "/" << p.genome_.countActiveChromosomes() <<  std::endl;
  //  CHECK(fabs(fitnessBefore - p.fitness_) < 0.01);
  }

  Population make_population(const Mat& oImg, const Mat& pImg, const size_t size) {
    Population pop;
    for(size_t i = 0; i < size; ++i) {
      std::cerr << '\r' << i;
      pop.push_back(Painter(oImg, pImg));
    }
    std::cerr << std::endl;

    return pop;
  }

  void run_population(Population& pop) {
    for(size_t j = 0; j < pop.size(); ++j) {
      std::cerr << '\r' << j;
      Painter& p = pop[j];
      p.paint();
    }
    std::cerr << std::endl;
  }

  void run(const string& original, const string& processed, size_t prePopulationSize, size_t populationSize, size_t iterations) {
    Mat oImg = imread( original.c_str());
    Mat pImg = imread( processed.c_str());

    vector<Painter> population;
    GeneticLayout gl = make_default_genetic_layout();
    GeneticPool pool(gl);

    std::cerr << "create pre population" << std::endl;
    population = make_population(oImg, pImg, prePopulationSize);

    std::cerr << "run pre population" << std::endl;
    run_population(population);

    std::sort(population.rbegin(), population.rend());
    std::cerr << "best prepop fitness: " << population[0].fitness_ << std::endl;

    population.resize(populationSize);

    std::cerr << "run" << std::endl;

    for(size_t i = 0; i < iterations; ++i) {
      std::cerr << "i: " << i << std::endl;
      run_population(population);

      sort(population.begin(), population.end());

      for(size_t j = 0; j < population.size(); ++j) {
        Painter& p = population[j];
        stringstream ssname;
        ssname << "result/" << setfill('0') << setw(10) << i << "_" << setfill('0') << setw(5) << j;
        imwrite(ssname.str() + "_result.png",p.result_);
        normalize(p.resultDFT_, p.resultDFT_, 0, 255, CV_MINMAX);
        p.resultDFT_.convertTo(p.resultDFT_, CV_8U);
        imwrite(ssname.str() + "_dft.png",p.resultDFT_);
        std::cerr << "fitness: " << (std::to_string(p.fitness_) + " ( " + std::to_string(p.histError_) + " * " + std::to_string(p.pixelError_) + " * " + std::to_string(p.fftError_) + " )\t" + std::to_string(p.genome_.countActiveChromosomes()) + "/" + std::to_string(p.genome_.getTotalKernelSize()) + "\t" + ssname.str() + "\t");

        for(const Chromosome& c: p.genome_) {
          if(c.isActive())
            std::cerr << " " << c.getOperation();
        }

        std::cerr << std::endl;
      }

      PopulationStatistics stats = calcPopulationStats(population);
      std::cout << stats.avgFitness_ << ';' << stats.bestFitness_ << ';' << stats.avgPixelError_ << ';' << stats.avgHistError_ << ';' << stats.avgKernelSize_ << ';' << std::endl;

      updateGlobalStats(stats);
      if(false && global_stats.streak_ > 100) {
        global_stats.streak_ = 0;
        std::cerr << "### repopulate ###" << std::endl;
        Population freshPop = make_population(oImg, pImg, populationSize - 1);
        Painter& best = population.back();
        pack_genome(best);

        freshPop.push_back(best);
        population = freshPop;
      } else {
        population = pool.epoch(population);
      }
    }
  }
} /* namespace mimikry */


inline void default_error_delegate(const std::string& msg) {
  std::cerr << "### Error: " << msg << std::endl;
  mimikry::print_stacktrace();
  exit(1);
}

int main(int argc, char** argv) {
  std::cerr << "opencv optimized: " << cv::useOptimized() << std::endl;
  srand (time(NULL));
  if(argc != 3) {
    std::cerr << "Usage: mimikry <original image> <processed image>" << std::endl;
    exit(1);
  }

  mimikry::ErrorHandler::init(default_error_delegate);
  mimikry::run(std::string(argv[1]), std::string(argv[2]), 1000, 100, 10000000);
  return 0;
}
