
#include "../src/particleSamples.h"

#include <cmath>
#include <memory>
#include <string>

#include "../src/Random.h"
#include "gtest/gtest.h"
#include "zlib.h"

using std::string;

namespace {

TEST(ParticleDecay, ConservesFourMomentumForRestingAndMovingMothers) {
    particle_decay decayer(std::make_shared<RandomUtil::Random>(24680), 0);
    for (int pdg : {113, 223}) {
        for (double pz : {0.0, 0.6}) {
            particle_info mother{}, plus{}, minus{}, neutral{};
            mother.monval = pdg;
            mother.mass = decayer.get_particle_mass(pdg);
            mother.pz = pz;
            mother.E = std::hypot(mother.mass, mother.pz);
            plus.mass = decayer.get_particle_mass(211);
            minus.mass = decayer.get_particle_mass(-211);
            if (pdg == 113) {
                decayer.perform_two_body_decay(mother, plus, minus);
            } else {
                neutral.mass = decayer.get_particle_mass(111);
                decayer.perform_three_body_decay(mother, plus, minus, neutral);
            }

            EXPECT_NEAR(plus.E + minus.E + neutral.E, mother.E, 1e-12);
            EXPECT_NEAR(plus.px + minus.px + neutral.px, mother.px, 1e-12);
            EXPECT_NEAR(plus.py + minus.py + neutral.py, mother.py, 1e-12);
            EXPECT_NEAR(plus.pz + minus.pz + neutral.pz, mother.pz, 1e-12);
        }
    }
}

TEST(particleSamples, RapidityShiftPreservesHbtPhase) {
    ParameterReader params;
    params.readFromFile("parameters.dat");
    params.setVal("read_in_mode", 0);
    particleSamples samples(
        params, "test_reader_files", std::make_shared<RandomUtil::Random>(0));

    std::vector<particle_info> pair(2);
    for (auto &particle : pair) {
        particle.mass = 0.14;
        particle.px = 0.2;
        particle.E = std::sqrt(0.14 * 0.14 + 0.2 * 0.2 + 0.1 * 0.1);
    }
    pair[0].pz = 0.1;
    pair[1].pz = -0.1;
    pair[0].z = 1.0;
    std::vector<std::vector<particle_info> *> events{&pair};
    samples.boostParticles(&events, 0.7);
    const double phase =
        (pair[0].E - pair[1].E) * (pair[0].t - pair[1].t)
        - (pair[0].pz - pair[1].pz) * (pair[0].z - pair[1].z);
    EXPECT_NEAR(phase, -0.2, 1e-12);
}

TEST(particleSamples, DefaultConstructor) {
    ParameterReader paraRdr;
    paraRdr.readFromFile("parameters.dat");
    paraRdr.setVal("run_mode", 0);
    string path = "test_reader_files";
    int randomSeed = 0;
    std::shared_ptr<RandomUtil::Random> ran_gen_ptr(
        new RandomUtil::Random(randomSeed));
    particleSamples particle_list(paraRdr, path, ran_gen_ptr);

    EXPECT_EQ(0, 0);
}

TEST(particleSamples, filter_particles_from_events) {
    ParameterReader paraRdr;
    paraRdr.readFromFile("parameters.dat");
    paraRdr.setVal("run_mode", 0);
    paraRdr.setVal("read_in_mode", 0);
    paraRdr.setVal("particle_monval", 321);
    string path = "test_reader_files";
    int randomSeed = 0;
    std::shared_ptr<RandomUtil::Random> ran_gen_ptr(
        new RandomUtil::Random(randomSeed));
    particleSamples particle_list(paraRdr, path, ran_gen_ptr);
    particle_list.read_in_particle_samples();
    particle_list.filter_particles_from_events(321);
    int nev = particle_list.get_number_of_events();
    EXPECT_EQ(nev, 2);

    int n_particles = particle_list.get_number_of_particles(0);
    EXPECT_EQ(n_particles, 11);
    n_particles = particle_list.get_number_of_particles(1);
    EXPECT_EQ(n_particles, 14);

    paraRdr.setVal("particle_monval", 2212);
    particle_list.filter_particles_from_events(2212);
    nev = particle_list.get_number_of_events();
    EXPECT_EQ(nev, 2);

    n_particles = particle_list.get_number_of_particles(0);
    EXPECT_EQ(n_particles, 1);
    n_particles = particle_list.get_number_of_particles(1);
    EXPECT_EQ(n_particles, 3);
}

TEST(particleSamples, read_in_particle_samples_OSCAR) {
    ParameterReader paraRdr;
    paraRdr.readFromFile("parameters.dat");
    paraRdr.setVal("run_mode", 0);
    paraRdr.setVal("read_in_mode", 0);
    paraRdr.setVal("particle_monval", 211);
    string path = "test_reader_files";
    int randomSeed = 0;
    std::shared_ptr<RandomUtil::Random> ran_gen_ptr(
        new RandomUtil::Random(randomSeed));
    particleSamples particle_list(paraRdr, path, ran_gen_ptr);
    particle_list.read_in_particle_samples_and_filter();
    int nev = particle_list.get_number_of_events();
    EXPECT_EQ(nev, 2);

    int n_particles = particle_list.get_number_of_particles(0);
    EXPECT_EQ(n_particles, 31);
    n_particles = particle_list.get_number_of_particles(1);
    EXPECT_EQ(n_particles, 27);
}

TEST(particleSamples, read_in_particle_samples_OSCAR_mixed_event) {
    ParameterReader paraRdr;
    paraRdr.readFromFile("parameters.dat");
    paraRdr.setVal("run_mode", 1);
    paraRdr.setVal("read_in_mode", 0);
    paraRdr.setVal("particle_monval", 211);
    string path = "test_reader_files";
    int randomSeed = 0;
    std::shared_ptr<RandomUtil::Random> ran_gen_ptr(
        new RandomUtil::Random(randomSeed));
    particleSamples particle_list(paraRdr, path, ran_gen_ptr);
    particle_list.read_in_particle_samples_mixed_event_and_filter();
    int nev = particle_list.get_number_of_mixed_events();
    EXPECT_EQ(nev, 2);

    int n_particles = particle_list.get_number_of_particles_mixed_event(0);
    EXPECT_EQ(n_particles, 31);
    n_particles = particle_list.get_number_of_particles_mixed_event(1);
    EXPECT_EQ(n_particles, 27);
}

TEST(particleSamples, read_in_particle_samples_UrQMD) {
    ParameterReader paraRdr;
    paraRdr.readFromFile("parameters.dat");
    paraRdr.setVal("run_mode", 0);
    paraRdr.setVal("read_in_mode", 1);
    paraRdr.setVal("particle_monval", 211);
    string path = "test_reader_files";
    int randomSeed = 0;
    std::shared_ptr<RandomUtil::Random> ran_gen_ptr(
        new RandomUtil::Random(randomSeed));
    particleSamples particle_list(paraRdr, path, ran_gen_ptr);
    particle_list.read_in_particle_samples();
    particle_list.filter_particles_from_events(211);
    int nev = particle_list.get_number_of_events();
    EXPECT_EQ(nev, 2);

    int n_particles = particle_list.get_number_of_particles(0);
    EXPECT_EQ(n_particles, 116);
    n_particles = particle_list.get_number_of_particles(1);
    EXPECT_EQ(n_particles, 97);
}

TEST(particleSamples, read_in_particle_samples_UrQMD_mixed_event) {
    ParameterReader paraRdr;
    paraRdr.readFromFile("parameters.dat");
    paraRdr.setVal("run_mode", 1);
    paraRdr.setVal("read_in_mode", 1);
    paraRdr.setVal("particle_monval", 211);
    string path = "test_reader_files";
    int randomSeed = 0;
    std::shared_ptr<RandomUtil::Random> ran_gen_ptr(
        new RandomUtil::Random(randomSeed));
    particleSamples particle_list(paraRdr, path, ran_gen_ptr);
    particle_list.read_in_particle_samples_mixed_event_and_filter();
    int nev = particle_list.get_number_of_mixed_events();
    EXPECT_EQ(nev, 2);

    int n_particles = particle_list.get_number_of_particles_mixed_event(0);
    EXPECT_EQ(n_particles, 116);
    n_particles = particle_list.get_number_of_particles_mixed_event(1);
    EXPECT_EQ(n_particles, 97);
}

TEST(particleSamples, read_in_particle_samples_UrQMD_zipped) {
    ParameterReader paraRdr;
    paraRdr.readFromFile("parameters.dat");
    paraRdr.setVal("run_mode", 0);
    paraRdr.setVal("read_in_mode", 2);
    paraRdr.setVal("particle_monval", 211);
    string path = "test_gzip_reader";
    int randomSeed = 0;
    std::shared_ptr<RandomUtil::Random> ran_gen_ptr(
        new RandomUtil::Random(randomSeed));
    particleSamples particle_list(paraRdr, path, ran_gen_ptr);
    particle_list.read_in_particle_samples();
    particle_list.filter_particles_from_events(211);
    int nev = particle_list.get_number_of_events();
    EXPECT_EQ(nev, 2);

    int n_particles = particle_list.get_number_of_particles(0);
    EXPECT_EQ(n_particles, 116);
    n_particles = particle_list.get_number_of_particles(1);
    EXPECT_EQ(n_particles, 97);
}

TEST(particleSamples, read_in_particle_samples_UrQMD_mixed_event_zipped) {
    ParameterReader paraRdr;
    paraRdr.readFromFile("parameters.dat");
    paraRdr.setVal("run_mode", 1);
    paraRdr.setVal("read_in_mode", 2);
    paraRdr.setVal("particle_monval", 211);
    string path = "test_gzip_reader";
    int randomSeed = 0;
    std::shared_ptr<RandomUtil::Random> ran_gen_ptr(
        new RandomUtil::Random(randomSeed));
    particleSamples particle_list(paraRdr, path, ran_gen_ptr);
    particle_list.read_in_particle_samples_mixed_event_and_filter();
    int nev = particle_list.get_number_of_mixed_events();
    EXPECT_EQ(nev, 2);

    int n_particles = particle_list.get_number_of_particles_mixed_event(0);
    EXPECT_EQ(n_particles, 116);
    n_particles = particle_list.get_number_of_particles_mixed_event(1);
    EXPECT_EQ(n_particles, 97);
}

TEST(particleSamples, read_in_particle_samples_gzipped) {
    ParameterReader paraRdr;
    paraRdr.readFromFile("parameters.dat");
    paraRdr.setVal("run_mode", 0);
    paraRdr.setVal("read_in_mode", 10);
    paraRdr.setVal("particle_monval", 211);
    string path = "test_gzip_reader";
    int randomSeed = 0;
    std::shared_ptr<RandomUtil::Random> ran_gen_ptr(
        new RandomUtil::Random(randomSeed));
    particleSamples particle_list(paraRdr, path, ran_gen_ptr);
    particle_list.read_in_particle_samples_and_filter();
    int nev = particle_list.get_number_of_events();
    EXPECT_EQ(nev, 2);

    int n_particles = particle_list.get_number_of_particles(0);
    EXPECT_EQ(n_particles, 619);
    n_particles = particle_list.get_number_of_particles(1);
    EXPECT_EQ(n_particles, 635);
}

TEST(particleSamples, read_in_particle_samples_mixed_event_gzipped) {
    ParameterReader paraRdr;
    paraRdr.readFromFile("parameters.dat");
    paraRdr.setVal("run_mode", 1);
    paraRdr.setVal("read_in_mode", 10);
    paraRdr.setVal("particle_monval", 211);
    string path = "test_gzip_reader";
    int randomSeed = 0;
    std::shared_ptr<RandomUtil::Random> ran_gen_ptr(
        new RandomUtil::Random(randomSeed));
    particleSamples particle_list(paraRdr, path, ran_gen_ptr);
    particle_list.read_in_particle_samples_mixed_event_and_filter();
    int nev = particle_list.get_number_of_mixed_events();
    EXPECT_EQ(nev, 2);

    int n_particles = particle_list.get_number_of_particles_mixed_event(0);
    EXPECT_EQ(n_particles, 619);
    n_particles = particle_list.get_number_of_particles_mixed_event(1);
    EXPECT_EQ(n_particles, 635);
}

TEST(particleSamples, perform_resonance_feed_down) {
    ParameterReader paraRdr;
    paraRdr.readFromFile("parameters.dat");
    paraRdr.setVal("run_mode", 0);
    paraRdr.setVal("read_in_mode", 10);
    paraRdr.setVal("particle_monval", 211);
    paraRdr.setVal("resonance_feed_down_flag", 1);
    string path = "test_gzip_reader";
    int randomSeed = 0;
    std::shared_ptr<RandomUtil::Random> ran_gen_ptr(
        new RandomUtil::Random(randomSeed));
    particleSamples particle_list(paraRdr, path, ran_gen_ptr);
    particle_list.read_in_particle_samples_and_filter();
    int nev = particle_list.get_number_of_events();
    EXPECT_EQ(nev, 2);

    int n_particles = particle_list.get_number_of_particles(0);
    EXPECT_EQ(n_particles, 744);
    n_particles = particle_list.get_number_of_particles(1);
    EXPECT_EQ(n_particles, 780);
}

TEST(particleSamples, map_urqmd_to_pdg) {
    ParameterReader paraRdr;
    paraRdr.readFromFile("parameters.dat");
    paraRdr.setVal("run_mode", 0);
    string path = "test_reader_files";
    int randomSeed = 0;
    std::shared_ptr<RandomUtil::Random> ran_gen_ptr(
        new RandomUtil::Random(randomSeed));
    particleSamples particle_list(paraRdr, path, ran_gen_ptr);
    EXPECT_EQ(particle_list.get_pdg_id(1, 1), 2212);
    EXPECT_EQ(particle_list.get_pdg_id(101, 2), 211);
    EXPECT_EQ(particle_list.get_pdg_id(40, 2), 3222);
    EXPECT_EQ(particle_list.get_pdg_id(23, 0), 0);
}

}  // namespace
