// Copyright @ Chun Shen 2018

#include "Random.h"

namespace RandomUtil {

Random::Random(int seed_in, double min, double max)
    : rand_uniform_dist(min, max) {
    seed = seed_in;
    if (seed == -1) {
        seed = std::random_device{}();
    }
    ran_generator = std::unique_ptr<std::mt19937>(new std::mt19937(seed));
}

Random::Random(const std::array<std::uint32_t, 3> &stream_key)
    : seed(static_cast<int>(stream_key.back())), rand_uniform_dist(0.0, 1.0) {
    std::seed_seq keys(stream_key.begin(), stream_key.end());
    ran_generator.reset(new std::mt19937(keys));
}

}  // namespace RandomUtil
