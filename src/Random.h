// Copyright @ Chun Shen 2018

#ifndef SRC_RANDOM_H_
#define SRC_RANDOM_H_

#include <memory>
#include <random>

namespace RandomUtil {

class Random {
  private:
    int seed;
    std::random_device ran_dev;
    std::unique_ptr<std::mt19937> ran_generator;
    std::uniform_real_distribution<double> rand_uniform_dist;
    std::uniform_int_distribution<int> rand_int_uniform_dist;

  public:
    Random(int seed_in, double min = 0.0, double max = 1.0);
    double rand_uniform() { return (rand_uniform_dist(*ran_generator)); }
    int rand_int_uniform() { return (rand_int_uniform_dist(*ran_generator)); }
    int get_seed() const { return (seed); }
};

}  // namespace RandomUtil

#endif  // SRC_RANDOM_H_
