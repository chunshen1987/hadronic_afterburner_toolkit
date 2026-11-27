#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "../src/HBT_correlation.h"
#include "../src/parameters.h"
#include "binary_event_writer.h"
#include "gtest/gtest.h"

using std::string;

namespace {
// Tests the Glauber class

TEST(HBT_correlation, DefaultConstructor) {
    ParameterReader paraRdr;
    paraRdr.readFromFile("parameters.dat");
    paraRdr.setVal("run_mode", 1);
    string path = "test_gzip_reader";
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
    string path = "test_gzip_reader";
    int randomSeed = 0;
    std::shared_ptr<RandomUtil::Random> ran_gen_ptr(
        new RandomUtil::Random(randomSeed));
    auto particle_list =
        std::make_shared<particleSamples>(paraRdr, path, ran_gen_ptr);
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
    string path = "test_gzip_reader";
    int randomSeed = 0;
    std::shared_ptr<RandomUtil::Random> ran_gen_ptr(
        new RandomUtil::Random(randomSeed));
    auto particle_list =
        std::make_shared<particleSamples>(paraRdr, path, ran_gen_ptr);
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
    string path = "test_gzip_reader";
    int randomSeed = 0;
    std::shared_ptr<RandomUtil::Random> ran_gen_ptr(
        new RandomUtil::Random(randomSeed));
    auto particle_list =
        std::make_shared<particleSamples>(paraRdr, path, ran_gen_ptr);
    particle_list->read_in_particle_samples_and_filter();
    particle_list->read_in_particle_samples_mixed_event_and_filter();
    HBT_correlation test(paraRdr, path, ran_gen_ptr);
    test.calculate_HBT_correlation_function(particle_list);
    EXPECT_EQ(0, 0);
}

TEST(HBT_correlation, BinsOnePairAtItsOutSideLongMomentumAndPhase) {
    // One event with a pi+ pair at K_T = 0.3 GeV, azimuth 0.6 rad and
    // nonzero K_z, built with q_out = 0.04 GeV, q_side = -0.03 GeV and lab
    // q_z = 0.042 GeV. A third pi+ is too far away in momentum to pair.
    using test_fixtures::Track;
    const double m = 0.13957, KT = 0.3, phi_K = 0.6, K_z = 0.11;
    const double q_out = 0.04, q_side = -0.03, q_z = 0.042;
    const double c = std::cos(phi_K), s = std::sin(phi_K);
    const double K[3] = {KT * c, KT * s, K_z};
    const double q[3] = {q_out * c - q_side * s, q_out * s + q_side * c, q_z};
    const Track first {211,
                       m,
                       9.0,
                       3.0,
                       -2.0,
                       4.0,
                       K[0] + q[0] / 2,
                       K[1] + q[1] / 2,
                       K[2] + q[2] / 2};
    const Track second {211,
                        m,
                        8.0,
                        0.0,
                        0.0,
                        0.0,
                        K[0] - q[0] / 2,
                        K[1] - q[1] / 2,
                        K[2] - q[2] / 2};
    const Track far = test_fixtures::track_from_rapidity(211, m, 0.9, 0.1, 3.5);
    const std::string path = "test_hbt_single_pair";
    test_fixtures::write_binary_events(path, {{first, second, far}});

    ParameterReader params;
    params.readFromFile("parameters.dat");
    params.setVal("read_in_mode", 9);
    params.setVal("read_in_real_mixed_events", 0);
    params.setVal("particle_monval", 211);
    params.setVal("number_of_oversample_events", 1);
    params.setVal("invariant_radius_flag", 0);
    params.setVal("azimuthal_flag", 0);
    params.setVal("long_comoving_boost", 1);
    auto random = std::make_shared<RandomUtil::Random>(8);
    auto samples = std::make_shared<particleSamples>(params, path, random);
    samples->read_in_particle_samples_and_filter();
    ASSERT_EQ(samples->get_number_of_particles(0), 3);

    // Expected values from the stored (float-rounded) tracks, with
    // q = p1 - p2, K = (p1 + p2)/2, out along K_T, side = z-hat x out, long
    // in the frame where K_z = 0, and phase q.x / hbar c with metric (+---).
    const particle_info p1 = samples->get_particle(0, 0);
    const particle_info p2 = samples->get_particle(0, 1);
    const double qx = p1.px - p2.px, qy = p1.py - p2.py, qz = p1.pz - p2.pz;
    const double qE = p1.E - p2.E;
    const double Kx = 0.5 * (p1.px + p2.px), Ky = 0.5 * (p1.py + p2.py);
    const double Kz = 0.5 * (p1.pz + p2.pz), KE = 0.5 * (p1.E + p2.E);
    const double kT = std::hypot(Kx, Ky);
    const double expected_out = (qx * Kx + qy * Ky) / kT;
    const double expected_side = (qy * Kx - qx * Ky) / kT;
    const double beta = Kz / KE, gamma = 1.0 / std::sqrt(1.0 - beta * beta);
    const double expected_long = gamma * (qz - beta * qE);
    const double phase = (qE * (p1.t - p2.t) - qx * (p1.x - p2.x)
                          - qy * (p1.y - p2.y) - qz * (p1.z - p2.z))
                         / hbarC;
    ASSERT_NEAR(expected_out, q_out, 1e-6);
    ASSERT_NEAR(expected_side, q_side, 1e-6);
    ASSERT_GT(std::abs(expected_long - qz), 0.01);  // frames differ by > 1 bin

    HBT_correlation hbt(params, path, random);
    hbt.calculate_HBT_correlation_function(samples);
    hbt.output_HBTcorrelation();

    // Bins of width 0.01 GeV centred on -0.15, ..., 0.15 (31 points).
    const auto centre = [](double value) {
        const int index = static_cast<int>((value + 0.155) / 0.01);
        return -0.15 + 0.01 * index;
    };
    std::ifstream raw(
        (path + "/HBT_raw_correlation_function_KT_0.2_0.4.dat").c_str());
    ASSERT_TRUE(raw.good());
    std::string line;
    double total_same = 0.0, total_mixed = 0.0;
    int populated = 0;
    while (std::getline(raw, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream fields(line);
        double out, side, longitudinal, same, cos_sum, mixed;
        fields >> out >> side >> longitudinal >> same >> cos_sum >> mixed;
        total_same += same;
        total_mixed += mixed;
        if (same == 0.0) continue;
        populated++;
        EXPECT_NEAR(out, centre(expected_out), 1e-9);
        EXPECT_NEAR(side, centre(expected_side), 1e-9);
        EXPECT_NEAR(longitudinal, centre(expected_long), 1e-9);
        EXPECT_DOUBLE_EQ(same, 1.0);
        EXPECT_NEAR(cos_sum, std::cos(phase), 1e-12);
    }
    EXPECT_EQ(populated, 1);
    EXPECT_DOUBLE_EQ(total_same, 1.0);
    EXPECT_DOUBLE_EQ(total_mixed, 0.0);  // no other collision to mix with
}

}  // namespace
