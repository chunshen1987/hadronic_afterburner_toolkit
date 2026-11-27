
#include "../src/particleSamples.h"

#include <sys/stat.h>

#include <cmath>
#include <fstream>
#include <memory>
#include <stdexcept>
#include <string>

#include "../src/Random.h"
#include "binary_event_writer.h"
#include "gtest/gtest.h"
#include "zlib.h"

#ifdef _OPENMP
#include <omp.h>
#endif

using std::string;

namespace {

TEST(ParticleDecay, ConservesFourMomentumForRestingAndMovingMothers) {
    particle_decay decayer(std::make_shared<RandomUtil::Random>(24680), 0);
    for (int pdg : {113, 223}) {
        for (double pz : {0.0, 0.6}) {
            particle_info mother {}, plus {}, minus {}, neutral {};
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

TEST(ParticleDecay, ConservesRoundedFourMomentumAtLargeBoost) {
    particle_decay decayer(std::make_shared<RandomUtil::Random>(24680), 0);
    for (int pdg : {113, 223, 225}) {
        int four_body = 0;
        for (double gamma : {3.0, 30.0, 300.0}) {
            particle_info mother {};
            mother.monval = pdg;
            mother.mass = decayer.get_particle_mass(pdg);
            const double momentum = mother.mass * std::sqrt(gamma * gamma - 1);
            mother.E = static_cast<float>(gamma * mother.mass);
            mother.px = static_cast<float>(0.36 * momentum);
            mother.py = static_cast<float>(-0.48 * momentum);
            mother.pz = static_cast<float>(0.80 * momentum);

            for (int sample = 0; sample < 30; sample++) {
                std::vector<particle_info> daughters;
                decayer.perform_decays(mother, daughters);
                if (daughters.size() == 4) four_body++;
                double energy = 0, px = 0, py = 0, pz = 0;
                for (const auto &daughter : daughters) {
                    energy += daughter.E;
                    px += daughter.px;
                    py += daughter.py;
                    pz += daughter.pz;
                    const double mass_squared =
                        daughter.E * daughter.E - daughter.px * daughter.px
                        - daughter.py * daughter.py - daughter.pz * daughter.pz;
                    EXPECT_NEAR(
                        mass_squared, daughter.mass * daughter.mass, 1e-10);
                }
                const double tolerance = 1e-12 * mother.E;
                EXPECT_NEAR(energy, mother.E, tolerance);
                EXPECT_NEAR(px, mother.px, tolerance);
                EXPECT_NEAR(py, mother.py, tolerance);
                EXPECT_NEAR(pz, mother.pz, tolerance);
            }
        }
        if (pdg == 225) {
            EXPECT_GT(four_body, 0);
        }
    }
}

TEST(ParticleDecay, SamplesFourBodyAndSubthresholdChannels) {
    particle_decay decayer(std::make_shared<RandomUtil::Random>(12345), 0);
    int four_body = 0;
    for (int pdg : {225, 10333}) {
        particle_info mother {};
        mother.monval = pdg;
        mother.mass = decayer.get_particle_mass(pdg);
        mother.pz = 0.4;
        mother.E = std::hypot(mother.mass, mother.pz);
        for (int sample = 0; sample < 100; sample++) {
            std::vector<particle_info> daughters;
            decayer.perform_decays(mother, daughters);
            EXPECT_TRUE(daughters.size() == 2 || daughters.size() == 4);
            if (daughters.size() == 4) four_body++;
            double energy = 0.0, px = 0.0, py = 0.0, pz = 0.0;
            int charge = 0, strangeness = 0;
            bool sampled_mass = false;
            for (const auto &daughter : daughters) {
                energy += daughter.E;
                px += daughter.px;
                py += daughter.py;
                pz += daughter.pz;
                charge += decayer.get_particle_charge(daughter.monval);
                strangeness +=
                    decayer.get_particle_strange_number(daughter.monval);
                if (daughter.mass
                    < decayer.get_particle_mass(daughter.monval)) {
                    sampled_mass = true;
                    EXPECT_EQ(decayer.check_particle_stable(daughter), 0);
                }
            }
            if (pdg == 10333) {
                EXPECT_TRUE(sampled_mass);
            }
            EXPECT_NEAR(energy, mother.E, 1e-12);
            EXPECT_NEAR(px, mother.px, 1e-12);
            EXPECT_NEAR(py, mother.py, 1e-12);
            EXPECT_NEAR(pz, mother.pz, 1e-12);
            EXPECT_EQ(charge, 0);
            EXPECT_EQ(strangeness, 0);
        }
    }
    EXPECT_GT(four_body, 0);
}

TEST(ParticleDecay, SkipsExactThresholdAndKeepsTerminalResonance) {
    particle_decay decayer(std::make_shared<RandomUtil::Random>(12345), 0);
    particle_info mother {};
    mother.monval = 225;
    mother.E = mother.mass = 4 * decayer.get_particle_mass(111);
    for (int sample = 0; sample < 30; sample++) {
        std::vector<particle_info> daughters;
        decayer.perform_decays(mother, daughters);
        EXPECT_EQ(daughters.size(), 2u);
    }

    mother.monval = 113;
    mother.E = mother.mass = 2 * decayer.get_particle_mass(211);
    std::vector<particle_info> daughters;
    decayer.perform_decays(mother, daughters);
    ASSERT_EQ(daughters.size(), 1u);
    EXPECT_EQ(daughters.front().monval, mother.monval);
    EXPECT_EQ(decayer.check_particle_stable(daughters.front()), 1);
}

TEST(ParticleDecay, HandlesExactThreeBodyThresholdAndRejectsInconsistentInput) {
    particle_decay decayer(std::make_shared<RandomUtil::Random>(12345), 0);
    particle_info mother {}, first {}, second {}, third {};
    mother.monval = 223;
    mother.E = mother.mass = 0.375;
    first.mass = second.mass = third.mass = 0.125;
    decayer.perform_three_body_decay(mother, first, second, third);
    for (const auto *daughter : {&first, &second, &third}) {
        EXPECT_DOUBLE_EQ(daughter->E, daughter->mass);
        EXPECT_DOUBLE_EQ(daughter->px, 0);
        EXPECT_DOUBLE_EQ(daughter->py, 0);
        EXPECT_DOUBLE_EQ(daughter->pz, 0);
    }

    mother.monval = 113;
    mother.mass = 0.30;
    mother.E = 0.25;
    first.mass = second.mass = 0.138;
    EXPECT_THROW(
        decayer.perform_two_body_decay(mother, first, second),
        std::runtime_error);

    const double threshold = first.mass + second.mass;
    mother.mass = std::nextafter(threshold, 1.0);
    mother.E = std::nextafter(threshold, 0.0);
    decayer.perform_two_body_decay(mother, first, second);
    EXPECT_NEAR(first.E + second.E, mother.mass, 1e-15);
    EXPECT_NEAR(
        first.E * first.E - first.px * first.px - first.py * first.py
            - first.pz * first.pz,
        first.mass * first.mass, 1e-15);
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
    std::vector<std::vector<particle_info> *> events {&pair};
    samples.boostParticles(&events, 0.7);
    const double phase = (pair[0].E - pair[1].E) * (pair[0].t - pair[1].t)
                         - (pair[0].pz - pair[1].pz) * (pair[0].z - pair[1].z);
    EXPECT_NEAR(phase, -0.2, 1e-12);
}

TEST(particleSamples, ResonanceDecayStreamsAreIndependentOfWorkerCount) {
    ParameterReader params;
    params.readFromFile("parameters.dat");
    params.setVal("read_in_mode", 0);
    params.setVal("resonance_feed_down_flag", 1);
    params.setVal("resonance_decay_parallel", 1);
    const auto run = [&params](int threads) {
#ifdef _OPENMP
        const int old_threads = omp_get_max_threads();
        omp_set_num_threads(threads);
#else
        (void)threads;
#endif
        particleSamples samples(
            params, "test_reader_files",
            std::make_shared<RandomUtil::Random>(4));
        particle_info mother {};
        mother.monval = 225;
        mother.E = mother.mass = 1.275;
        std::vector<std::vector<particle_info>> events(
            8, std::vector<particle_info>(32, mother));
        std::vector<std::vector<particle_info> *> event_ptrs;
        for (auto &event : events) event_ptrs.push_back(&event);
        samples.perform_resonance_feed_down(&event_ptrs);
#ifdef _OPENMP
        omp_set_num_threads(old_threads);
#endif
        return events;
    };

    const auto one = run(1);
    const auto four = run(4);
    EXPECT_NE(one[0][0].px, one[1][0].px);
    for (size_t i = 0; i < one.size(); i++) {
        ASSERT_EQ(one[i].size(), four[i].size());
        for (size_t j = 0; j < one[i].size(); j++) {
            const auto &a = one[i][j], &b = four[i][j];
            EXPECT_EQ(a.monval, b.monval);
            EXPECT_EQ(a.E, b.E);
            EXPECT_EQ(a.px, b.px);
            EXPECT_EQ(a.py, b.py);
            EXPECT_EQ(a.pz, b.pz);
            EXPECT_EQ(a.t, b.t);
        }
    }
}

TEST(particleSamples, SigmaFeedDownSharesOneDecayAcrossLists) {
    // Two read_in_mode = 9 events with two Sigma0 each.
    const string path = "test_sigma_feed_down";
    mkdir(path.c_str(), 0755);
    {
        std::ofstream out(path + "/particle_list.bin", std::ios::binary);
        const float mass = 1.192642f, pT = 0.8f;
        for (int iev = 0; iev < 2; iev++) {
            const int count = 2;
            out.write(reinterpret_cast<const char *>(&count), sizeof(int));
            for (int i = 0; i < count; i++) {
                const int pdg = 3212;
                const double y = 0.3 * (2 * iev + i) - 0.4;
                const double mT = std::sqrt(mass * mass + pT * pT);
                const float record[9] = {
                    mass, 0.f, 0.f,
                    0.f,  0.f, static_cast<float>(mT * std::cosh(y)),
                    pT,   0.f, static_cast<float>(mT * std::sinh(y))};
                out.write(reinterpret_cast<const char *>(&pdg), sizeof(int));
                out.write(
                    reinterpret_cast<const char *>(record), sizeof(record));
            }
        }
    }
    ParameterReader params;
    params.readFromFile("parameters.dat");
    params.setVal("read_in_mode", 9);
    params.setVal("read_in_real_mixed_events", 0);
    params.setVal("resonance_feed_down_flag", 0);
    params.setVal("net_particle_flag", 0);
    params.setVal("analyze_balance_function", 1);
    params.setVal("particle_monval", 3122);
    params.setVal("particle_alpha", 3122);
    params.setVal("particle_beta", -3122);
    params.setVal("resonance_weak_feed_down_Sigma_to_Lambda_flag", 1);
    particleSamples samples(
        params, path, std::make_shared<RandomUtil::Random>(5));
    samples.read_in_particle_samples_and_filter();
    samples.read_in_particle_samples_mixed_event_and_filter();

    // particle_beta = -3122 makes bbar a Lambda list as well.
    const auto *a = samples.get_balance_function_particle_list_a();
    const auto *bbar = samples.get_balance_function_particle_list_bbar();
    const auto *a_mixed =
        samples.get_balance_function_particle_list_a_mixed_event();
    ASSERT_EQ(samples.get_number_of_events(), 2);
    ASSERT_EQ(a->size(), 2u);
    for (int iev = 0; iev < 2; iev++) {
        ASSERT_EQ(samples.get_number_of_particles(iev), 2);
        ASSERT_EQ((*a)[iev]->size(), 2u);
        ASSERT_EQ((*bbar)[iev]->size(), 2u);
        ASSERT_EQ((*a_mixed)[iev]->size(), 2u);
        for (int i = 0; i < 2; i++) {
            const particle_info lambda = samples.get_particle(iev, i);
            EXPECT_EQ(lambda.monval, 3122);
            for (const auto *list : {a, bbar, a_mixed}) {
                const particle_info &other = (*(*list)[iev])[i];
                EXPECT_EQ(other.event_track_id, lambda.event_track_id);
                EXPECT_DOUBLE_EQ(other.E, lambda.E);
                EXPECT_DOUBLE_EQ(other.px, lambda.px);
                EXPECT_DOUBLE_EQ(other.py, lambda.py);
                EXPECT_DOUBLE_EQ(other.pz, lambda.pz);
            }
        }
    }
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

TEST(ParameterReader, ParsesAssignmentsAndLetsLaterOverridesWin) {
    ParameterReader params;
    params.phraseOneLine("  Particle_Monval =  211   # pi+");
    params.phraseOneLine("# a comment line");
    params.phraseOneLine("   ");
    params.phraseOneLine("rap_max=0.5");
    EXPECT_TRUE(params.exist("particle_monval"));
    EXPECT_TRUE(params.exist("PARTICLE_MONVAL"));  // names are lower-cased
    EXPECT_DOUBLE_EQ(params.getVal("particle_monval"), 211);
    EXPECT_DOUBLE_EQ(params.getVal("rap_max"), 0.5);

    // Command-line key=value overrides replace file values.
    char program[] = "hadronic_afterburner_tools.e";
    char monval[] = "particle_monval=321";
    char rapidity[] = "rap_max = 1.0#note";
    char *argv[] = {program, monval, rapidity};
    params.readFromArguments(3, argv);
    EXPECT_DOUBLE_EQ(params.getVal("particle_monval"), 321);
    EXPECT_DOUBLE_EQ(params.getVal("rap_max"), 1.0);

    EXPECT_FALSE(params.exist("not_a_parameter"));
    EXPECT_DOUBLE_EQ(params.getVal("not_a_parameter", -7.5), -7.5);
    params.setVal("Not_A_Parameter", 2);
    EXPECT_DOUBLE_EQ(params.getVal("not_a_parameter", -7.5), 2);
}

TEST(particleSamples, BinaryReaderMapsFieldsAndKeepsEmptyEvents) {
    using test_fixtures::Track;
    const string path = "test_binary_reader";
    // Event 0: one proton with distinct values in every field. Event 1 is
    // empty. Event 2: a pi+ pi- pair.
    const Track proton {2212, 0.938, 7.5, -1.25, 2.5, 0.75, 0.3, -0.4, 1.2};
    const Track pi_plus =
        test_fixtures::track_from_rapidity(211, 0.13957, 0.5, 0.2, 1.0);
    const Track pi_minus =
        test_fixtures::track_from_rapidity(-211, 0.13957, 0.7, -0.3, -2.0);
    test_fixtures::write_binary_events(
        path, {{proton}, {}, {pi_plus, pi_minus}});

    ParameterReader params;
    params.readFromFile("parameters.dat");
    params.setVal("read_in_mode", 9);
    params.setVal("read_in_real_mixed_events", 0);
    params.setVal("resonance_feed_down_flag", 0);
    params.setVal("particle_monval", 99999);  // keep every track
    particleSamples samples(
        params, path, std::make_shared<RandomUtil::Random>(1));
    samples.read_in_particle_samples_and_filter();

    ASSERT_EQ(samples.get_number_of_events(), 3);
    EXPECT_EQ(samples.get_number_of_particles(0), 1);
    EXPECT_EQ(samples.get_number_of_particles(1), 0);
    EXPECT_EQ(samples.get_number_of_particles(2), 2);

    const auto as_float = [](double value) {
        return static_cast<double>(static_cast<float>(value));
    };
    const particle_info p = samples.get_particle(0, 0);
    EXPECT_EQ(p.monval, 2212);
    EXPECT_EQ(p.event_track_id, 1u);
    EXPECT_DOUBLE_EQ(p.mass, as_float(proton.mass));
    EXPECT_DOUBLE_EQ(p.t, as_float(proton.t));
    EXPECT_DOUBLE_EQ(p.x, as_float(proton.x));
    EXPECT_DOUBLE_EQ(p.y, as_float(proton.y));
    EXPECT_DOUBLE_EQ(p.z, as_float(proton.z));
    EXPECT_DOUBLE_EQ(p.E, as_float(proton.E()));
    EXPECT_DOUBLE_EQ(p.px, as_float(proton.px));
    EXPECT_DOUBLE_EQ(p.py, as_float(proton.py));
    EXPECT_DOUBLE_EQ(p.pz, as_float(proton.pz));

    // Derived kinematics: pT, azimuth, y = asinh(pz/mT), eta = asinh(pz/pT).
    const double pT = std::hypot(p.px, p.py);
    EXPECT_DOUBLE_EQ(p.pT, pT);
    EXPECT_DOUBLE_EQ(p.phi_p, std::atan2(p.py, p.px));
    EXPECT_NEAR(p.rap_y, std::asinh(p.pz / std::hypot(pT, p.mass)), 1e-12);
    EXPECT_NEAR(p.rap_eta, std::asinh(p.pz / pT), 1e-12);

    const particle_info second = samples.get_particle(2, 1);
    EXPECT_EQ(second.monval, -211);
    EXPECT_EQ(second.event_track_id, 2u);
    EXPECT_NEAR(second.rap_y, -0.3, 1e-6);
    EXPECT_NEAR(second.pT, 0.7, 1e-6);
    EXPECT_NEAR(second.phi_p, -2.0, 1e-6);
}

TEST(particleSamples, AggregateSelectorsFollowQuantumNumbers) {
    // Quantum numbers come from EOS/pdg.dat (UrQMD table): K+ has S = +1,
    // Lambda S = -1, Xi- S = -2.
    ParameterReader params;
    params.readFromFile("parameters.dat");
    params.setVal("particle_monval", 211);
    particleSamples samples(
        params, "test_reader_files", std::make_shared<RandomUtil::Random>(0));

    struct Case {
        int selector, monval;
        bool picked;
    };
    const Case cases[] = {
        {99999, 211, true},   {99999, 111, true},  {99999, -3122, true},
        {9999, 211, true},    {9999, -211, true},  {9999, 2212, true},
        {9999, 111, false},   {9999, 3122, false}, {9998, 211, true},
        {9998, -211, false},  {9998, 111, false},  {-9998, -211, true},
        {-9998, 211, false},  {9997, 2212, true},  {9997, 3122, true},
        {9997, -2212, false}, {9997, 211, false},  {-9997, -2212, true},
        {-9997, 2212, false}, {9996, 321, true},   {9996, -321, false},
        {9996, 3122, false},  {-9996, 3122, true}, {-9996, 3312, true},
        {-9996, 321, false},  {321, 321, true},    {321, -321, false},
        {-211, -211, true},   {-211, 211, false},
    };
    for (const auto &c : cases) {
        EXPECT_EQ(samples.decide_to_pick_OSCAR(c.selector, c.monval), c.picked)
            << "selector " << c.selector << ", particle " << c.monval;
    }

    EXPECT_EQ(samples.decide_to_pick_charge(211), 1);
    EXPECT_EQ(samples.decide_to_pick_charge(-211), -1);
    EXPECT_EQ(samples.decide_to_pick_charge(111), 0);
    EXPECT_EQ(samples.decide_to_pick_charge(-2212), -1);

    // Antiparticle of an identified particle of interest (211).
    EXPECT_EQ(samples.decide_to_pick_anti_particles(-211), 1);
    EXPECT_EQ(samples.decide_to_pick_anti_particles(211), 0);
    EXPECT_EQ(samples.decide_to_pick_anti_particles(-321), 0);

    // For 9998 (positive hadrons) the antiparticles are negative hadrons.
    params.setVal("particle_monval", 9998);
    particleSamples positive(
        params, "test_reader_files", std::make_shared<RandomUtil::Random>(0));
    EXPECT_EQ(positive.decide_to_pick_anti_particles(-211), 1);
    EXPECT_EQ(positive.decide_to_pick_anti_particles(-2212), 1);
    EXPECT_EQ(positive.decide_to_pick_anti_particles(111), 0);
    EXPECT_EQ(positive.decide_to_pick_anti_particles(321), 0);
}

TEST(particleSamples, LongitudinalBoostShiftsRapiditiesAndKeepsInvariants) {
    ParameterReader params;
    params.readFromFile("parameters.dat");
    particleSamples samples(
        params, "test_reader_files", std::make_shared<RandomUtil::Random>(0));

    // A pi+ at y = 0.4 emitted at tau = 6 fm, eta_s = -0.25.
    const double mass = 0.13957, pT = 0.6, y = 0.4, phi = 2.2;
    const double tau = 6.0, eta_s = -0.25, shift = 0.7;
    const double mT = std::hypot(mass, pT);
    particle_info track {};
    track.monval = 211;
    track.mass = mass;
    track.px = pT * std::cos(phi);
    track.py = pT * std::sin(phi);
    track.pz = mT * std::sinh(y);
    track.E = mT * std::cosh(y);
    track.t = tau * std::cosh(eta_s);
    track.z = tau * std::sinh(eta_s);
    track.x = 1.5;
    track.y = -0.5;
    std::vector<particle_info> event {track};
    std::vector<std::vector<particle_info> *> events {&event};

    samples.boostParticles(&events, shift);
    const particle_info &b = event[0];
    EXPECT_NEAR(b.rap_y, y + shift, 1e-12);
    EXPECT_NEAR(
        0.5 * std::log((b.t + b.z) / (b.t - b.z)), eta_s + shift, 1e-12);
    EXPECT_NEAR(b.E * b.E - b.pz * b.pz, mT * mT, 1e-12);  // m_T invariant
    EXPECT_NEAR(b.t * b.t - b.z * b.z, tau * tau, 1e-11);  // tau invariant
    EXPECT_DOUBLE_EQ(b.pT, pT);
    EXPECT_NEAR(b.phi_p, phi, 1e-12);
    EXPECT_NEAR(b.rap_eta, std::asinh(b.pz / pT), 1e-12);
    EXPECT_DOUBLE_EQ(b.x, 1.5);
    EXPECT_DOUBLE_EQ(b.y, -0.5);

    // A zero shift leaves four-vectors unchanged but fills the kinematics.
    std::vector<particle_info> unboosted {track};
    std::vector<std::vector<particle_info> *> single {&unboosted};
    samples.boostParticles(&single, 0.0);
    EXPECT_DOUBLE_EQ(unboosted[0].pz, track.pz);
    EXPECT_DOUBLE_EQ(unboosted[0].t, track.t);
    EXPECT_NEAR(unboosted[0].rap_y, y, 1e-12);
}

}  // namespace
