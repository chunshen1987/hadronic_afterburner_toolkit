#include "../src/particleSamples.h"
#include "gtest/gtest.h"

namespace {
// Tests the Glauber class

TEST(particleSamples, DefaultConstructor) {
    ParameterReader *paraRdr = new ParameterReader;
    paraRdr->readFromFile("parameters.dat");
    paraRdr->setVal("run_mode", 0);
    string path="test_files";
    particleSamples particle_list(paraRdr, path);
        
    EXPECT_EQ(0, 0);
}

TEST(particleSamples, read_in_particle_samples_OSCAR) {
    ParameterReader *paraRdr = new ParameterReader;
    paraRdr->readFromFile("parameters.dat");
    paraRdr->setVal("run_mode", 0);
    paraRdr->setVal("read_in_mode", 0);
    paraRdr->setVal("particle_monval", 211);
    string path="test_reader_files";
    particleSamples particle_list(paraRdr, path);
    particle_list.read_in_particle_samples();
    int nev = particle_list.get_number_of_events();
    EXPECT_EQ(nev, 2);

    int n_particles = particle_list.get_number_of_particles(0);
    EXPECT_EQ(n_particles, 31);
    n_particles = particle_list.get_number_of_particles(1);
    EXPECT_EQ(n_particles, 27);
}

TEST(particleSamples, read_in_particle_samples_OSCAR_mixed_event) {
    ParameterReader *paraRdr = new ParameterReader;
    paraRdr->readFromFile("parameters.dat");
    paraRdr->setVal("run_mode", 1);
    paraRdr->setVal("read_in_mode", 0);
    paraRdr->setVal("particle_monval", 211);
    string path="test_reader_files";
    particleSamples particle_list(paraRdr, path);
    particle_list.read_in_particle_samples_OSCAR_mixed_event();
    int nev = particle_list.get_number_of_mixed_events();
    EXPECT_EQ(nev, 2);

    int n_particles = particle_list.get_number_of_particles_mixed_event(0);
    EXPECT_EQ(n_particles, 31);
    n_particles = particle_list.get_number_of_particles_mixed_event(1);
    EXPECT_EQ(n_particles, 27);
}

TEST(particleSamples, read_in_particle_samples_UrQMD) {
    ParameterReader *paraRdr = new ParameterReader;
    paraRdr->readFromFile("parameters.dat");
    paraRdr->setVal("run_mode", 0);
    paraRdr->setVal("read_in_mode", 1);
    paraRdr->setVal("particle_monval", 211);
    string path="test_reader_files";
    particleSamples particle_list(paraRdr, path);
    particle_list.read_in_particle_samples();
    int nev = particle_list.get_number_of_events();
    EXPECT_EQ(nev, 2);

    int n_particles = particle_list.get_number_of_particles(0);
    EXPECT_EQ(n_particles, 116);
    n_particles = particle_list.get_number_of_particles(1);
    EXPECT_EQ(n_particles, 97);
}

TEST(particleSamples, read_in_particle_samples_UrQMD_mixed_event) {
    ParameterReader *paraRdr = new ParameterReader;
    paraRdr->readFromFile("parameters.dat");
    paraRdr->setVal("run_mode", 1);
    paraRdr->setVal("read_in_mode", 1);
    paraRdr->setVal("particle_monval", 211);
    string path="test_reader_files";
    particleSamples particle_list(paraRdr, path);
    particle_list.read_in_particle_samples_UrQMD_mixed_event();
    int nev = particle_list.get_number_of_mixed_events();
    EXPECT_EQ(nev, 2);

    int n_particles = particle_list.get_number_of_particles_mixed_event(0);
    EXPECT_EQ(n_particles, 116);
    n_particles = particle_list.get_number_of_particles_mixed_event(1);
    EXPECT_EQ(n_particles, 97);
}

TEST(particleSamples, read_in_particle_samples_UrQMD_zipped) {
    ParameterReader *paraRdr = new ParameterReader;
    paraRdr->readFromFile("parameters.dat");
    paraRdr->setVal("run_mode", 0);
    paraRdr->setVal("read_in_mode", 2);
    paraRdr->setVal("particle_monval", 211);
    string path="test_gzip_reader";
    particleSamples particle_list(paraRdr, path);
    particle_list.read_in_particle_samples();
    int nev = particle_list.get_number_of_events();
    EXPECT_EQ(nev, 2);

    int n_particles = particle_list.get_number_of_particles(0);
    EXPECT_EQ(n_particles, 116);
    n_particles = particle_list.get_number_of_particles(1);
    EXPECT_EQ(n_particles, 97);
}

TEST(particleSamples, read_in_particle_samples_UrQMD_mixed_event_zipped) {
    ParameterReader *paraRdr = new ParameterReader;
    paraRdr->readFromFile("parameters.dat");
    paraRdr->setVal("run_mode", 1);
    paraRdr->setVal("read_in_mode", 2);
    paraRdr->setVal("particle_monval", 211);
    string path="test_gzip_reader";
    particleSamples particle_list(paraRdr, path);
    particle_list.read_in_particle_samples_UrQMD_mixed_event_zipped();
    int nev = particle_list.get_number_of_mixed_events();
    EXPECT_EQ(nev, 2);

    int n_particles = particle_list.get_number_of_particles_mixed_event(0);
    EXPECT_EQ(n_particles, 116);
    n_particles = particle_list.get_number_of_particles_mixed_event(1);
    EXPECT_EQ(n_particles, 97);
}

//    EXPECT_NEAR(1.5, z_test, 1e-8);
//    EXPECT_DOUBLE_EQ(2.6571870894737679, t_test);
}  // namespace
