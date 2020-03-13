
#include "../src/particleSamples.h"
#include "../src/Random.h"
#include "gtest/gtest.h"
#include "zlib.h"

#include <memory>
#include <string>

using std::string;

namespace {

TEST(particleSamples, DefaultConstructor) {
    ParameterReader paraRdr;
    paraRdr.readFromFile("parameters.dat");
    paraRdr.setVal("run_mode", 0);
    string path="test_reader_files";
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
    string path="test_reader_files";
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
    string path="test_reader_files";
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
    string path="test_reader_files";
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
    string path="test_reader_files";
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
    string path="test_reader_files";
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
    string path="test_gzip_reader";
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
    string path="test_gzip_reader";
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
    string path="test_gzip_reader";
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
    string path="test_gzip_reader";
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
    string path="test_gzip_reader";
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
    string path="test_reader_files";
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
