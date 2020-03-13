#include "../src/HBT_correlation.h"
#include "gtest/gtest.h"

#include <string>

using std::string;

namespace {
// Tests the Glauber class

TEST(HBT_correlation, DefaultConstructor) {
    ParameterReader paraRdr;
    paraRdr.readFromFile("parameters.dat");
    paraRdr.setVal("run_mode", 1);
    string path="test_gzip_reader";
    int randomSeed = 0;
    std::shared_ptr<RandomUtil::Random> ran_gen_ptr(
                                    new RandomUtil::Random(randomSeed));
    HBT_correlation test(paraRdr, path, ran_gen_ptr);
    EXPECT_EQ(0, 0);
}

TEST(HBT_correlation, calculate_flow_event_plane_angle) {
    ParameterReader paraRdr;
    paraRdr.readFromFile("parameters.dat");
    paraRdr.setVal("run_mode", 1);
    string path="test_gzip_reader";
    int randomSeed = 0;
    std::shared_ptr<RandomUtil::Random> ran_gen_ptr(
                                    new RandomUtil::Random(randomSeed));
    auto particle_list = std::make_shared<particleSamples> (paraRdr, path,
                                                            ran_gen_ptr);
    particle_list->read_in_particle_samples_and_filter();

    HBT_correlation test(paraRdr, path, ran_gen_ptr);
    test.set_particle_list(particle_list);
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
    ParameterReader paraRdr;
    paraRdr.readFromFile("parameters.dat");
    paraRdr.setVal("run_mode", 1);
    paraRdr.setVal("number_of_oversample_events", 2);
    paraRdr.setVal("number_of_mixed_events", 1);
    paraRdr.setVal("invariant_radius_flag", 1);
    string path="test_gzip_reader";
    int randomSeed = 0;
    std::shared_ptr<RandomUtil::Random> ran_gen_ptr(
                                    new RandomUtil::Random(randomSeed));
    auto particle_list = std::make_shared<particleSamples> (paraRdr, path,
                                                            ran_gen_ptr);
    particle_list->read_in_particle_samples_and_filter();
    particle_list->read_in_particle_samples_mixed_event_and_filter();
    HBT_correlation test(paraRdr, path, ran_gen_ptr);
    test.calculate_HBT_correlation_function(particle_list);
    EXPECT_EQ(0, 0);
}

TEST(HBT_correlation, calculate_HBT_correlation_function) {
    ParameterReader paraRdr;
    paraRdr.readFromFile("parameters.dat");
    paraRdr.setVal("run_mode", 1);
    paraRdr.setVal("number_of_oversample_events", 2);
    paraRdr.setVal("number_of_mixed_events", 1);
    paraRdr.setVal("invariant_radius_flag", 0);
    paraRdr.setVal("azimuthal_flag", 0);
    string path="test_gzip_reader";
    int randomSeed = 0;
    std::shared_ptr<RandomUtil::Random> ran_gen_ptr(
                                    new RandomUtil::Random(randomSeed));
    auto particle_list = std::make_shared<particleSamples> (paraRdr, path,
                                                            ran_gen_ptr);
    particle_list->read_in_particle_samples_and_filter();
    particle_list->read_in_particle_samples_mixed_event_and_filter();
    HBT_correlation test(paraRdr, path, ran_gen_ptr);
    test.calculate_HBT_correlation_function(particle_list);
    EXPECT_EQ(0, 0);
}

}  // namespace
