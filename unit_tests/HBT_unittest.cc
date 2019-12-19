#include "../src/HBT_correlation.h"
#include "gtest/gtest.h"

#include <string>

using std::string;

namespace {
// Tests the Glauber class

TEST(HBT_correlation, DefaultConstructor) {
    int randomSeed = 0;
    std::shared_ptr<RandomUtil::Random> ran_gen_ptr(
                                    new RandomUtil::Random(randomSeed));
    ParameterReader *paraRdr = new ParameterReader;
    paraRdr->readFromFile("parameters.dat");
    paraRdr->setVal("run_mode", 1);
    string path="test_gzip_reader";
    particleSamples particle_list(paraRdr, path, ran_gen_ptr);
    HBT_correlation test(paraRdr, path, ran_gen_ptr, &particle_list);
    EXPECT_EQ(0, 0);
}

TEST(HBT_correlation, calculate_flow_event_plane_angle) {
    int randomSeed = 0;
    std::shared_ptr<RandomUtil::Random> ran_gen_ptr(
                                    new RandomUtil::Random(randomSeed));
    ParameterReader *paraRdr = new ParameterReader;
    paraRdr->readFromFile("parameters.dat");
    paraRdr->setVal("run_mode", 1);
    string path="test_gzip_reader";
    particleSamples particle_list(paraRdr, path, ran_gen_ptr);
    particle_list.read_in_particle_samples();

    HBT_correlation test(paraRdr, path, ran_gen_ptr, &particle_list);
    test.calculate_flow_event_plane_angle(1);
    double psi_ref = test.get_psi_ref();
    EXPECT_NEAR(psi_ref, -1.6751628499713109, 1e-8);

    test.calculate_flow_event_plane_angle(2);
    psi_ref = test.get_psi_ref();
    EXPECT_NEAR(psi_ref, -0.85236384292521539, 1e-8);

    test.calculate_flow_event_plane_angle(3);
    psi_ref = test.get_psi_ref();
    EXPECT_NEAR(psi_ref, 0.06674083818478263, 1e-8);
}

TEST(HBT_correlation, calculate_HBT_correlation_function_inv) {
    int randomSeed = 0;
    std::shared_ptr<RandomUtil::Random> ran_gen_ptr(
                                    new RandomUtil::Random(randomSeed));
    ParameterReader *paraRdr = new ParameterReader;
    paraRdr->readFromFile("parameters.dat");
    paraRdr->setVal("run_mode", 1);
    paraRdr->setVal("number_of_oversample_events", 2);
    paraRdr->setVal("number_of_mixed_events", 1);
    paraRdr->setVal("invariant_radius_flag", 1);
    string path="test_gzip_reader";
    particleSamples particle_list(paraRdr, path, ran_gen_ptr);
    HBT_correlation test(paraRdr, path, ran_gen_ptr, &particle_list);
    test.calculate_HBT_correlation_function();
    EXPECT_EQ(0, 0);
}

TEST(HBT_correlation, calculate_HBT_correlation_function) {
    int randomSeed = 0;
    std::shared_ptr<RandomUtil::Random> ran_gen_ptr(
                                    new RandomUtil::Random(randomSeed));
    ParameterReader *paraRdr = new ParameterReader;
    paraRdr->readFromFile("parameters.dat");
    paraRdr->setVal("run_mode", 1);
    paraRdr->setVal("number_of_oversample_events", 2);
    paraRdr->setVal("number_of_mixed_events", 1);
    paraRdr->setVal("invariant_radius_flag", 0);
    paraRdr->setVal("azimuthal_flag", 0);
    string path="test_gzip_reader";
    particleSamples particle_list(paraRdr, path, ran_gen_ptr);
    HBT_correlation test(paraRdr, path, ran_gen_ptr, &particle_list);
    test.calculate_HBT_correlation_function();
    EXPECT_EQ(0, 0);
}

}  // namespace
