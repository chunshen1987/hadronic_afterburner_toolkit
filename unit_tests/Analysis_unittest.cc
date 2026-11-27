
#include "Analysis.h"

#include <cmath>
#include <memory>
#include <string>
#include <vector>

#include "BalanceFunction.h"
#include "single_particleSpectra.h"
#include "gtest/gtest.h"

namespace {

TEST(BalanceFunction, PreservesDistinctTracksAtEqualRapidity) {
    ParameterReader params;
    params.setVal("particle_alpha", 211);
    params.setVal("particle_beta", 211);
    params.setVal("BpT_min", 0.2);
    params.setVal("BpT_max", 3.0);
    params.setVal("Bnpts", 21);
    params.setVal("Brap_max", 2.0);
    params.setVal("rap_type", 1);
    auto ran = std::make_shared<RandomUtil::Random>(1);
    BalanceFunction analysis(params, ".", ran);
    particle_info first{}, second{};
    first.monval = second.monval = 211;
    first.pT = second.pT = 0.3;
    first.event_track_id = 1;
    second.event_track_id = 2;
    std::vector<particle_info> a{first, second}, b{first, second};
    std::vector<std::vector<particle_info> *> selected_a{&a}, selected_b{&b};
    std::vector<std::vector<double>> same(21, std::vector<double>(20, 0));
    analysis.combine_and_bin_particle_pairs(same, &selected_a, &selected_b);
    EXPECT_DOUBLE_EQ(same[10][5], 2);
}

TEST(BalanceFunction, SameBatchMixingNeverPairsAnEventWithItself) {
    ParameterReader params;
    params.setVal("particle_alpha", 211);
    params.setVal("particle_beta", -211);
    params.setVal("BpT_min", 0.2);
    params.setVal("BpT_max", 3.0);
    params.setVal("Bnpts", 21);
    params.setVal("Brap_max", 2.0);
    params.setVal("rap_type", 1);
    auto ran = std::make_shared<RandomUtil::Random>(1);
    BalanceFunction analysis(params, ".", ran);
    // One like-sign track per event, so a reference drawn from the same
    // event would pair each track with itself at Delta y = 0.
    particle_info first{}, second{};
    first.monval = second.monval = 211;
    first.pT = second.pT = 0.3;
    first.event_track_id = second.event_track_id = 1;
    first.rap_y = 0.0;
    second.rap_y = 1.0;
    std::vector<particle_info> event_0{first}, event_1{second};
    std::vector<std::vector<particle_info> *> events{&event_0, &event_1};
    for (int repeat = 0; repeat < 20; repeat++) {
        std::vector<std::vector<double>> mixed(21, std::vector<double>(20, 0));
        analysis.combine_and_bin_mixed_particle_pairs(
            mixed, &events, &events, true);
        double self_bin = 0, minus_one = 0, plus_one = 0;
        for (int j = 0; j < 20; j++) {
            self_bin += mixed[10][j];
            minus_one += mixed[5][j];
            plus_one += mixed[15][j];
        }
        EXPECT_DOUBLE_EQ(self_bin, 0);
        EXPECT_DOUBLE_EQ(minus_one, 1);
        EXPECT_DOUBLE_EQ(plus_one, 1);
    }

    // A single-event batch has no independent reference.
    std::vector<std::vector<particle_info> *> single{&event_0};
    std::vector<std::vector<double>> mixed(21, std::vector<double>(20, 0));
    analysis.combine_and_bin_mixed_particle_pairs(
        mixed, &single, &single, true);
    for (const auto &row : mixed) {
        for (double entry : row) EXPECT_DOUBLE_EQ(entry, 0);
    }
}

TEST(Analysis, Constructor) {
    std::string path = "results";
    Analysis test_analysis(path);
    EXPECT_EQ(0, 0);
}

TEST(Analysis, ReadInParameters) {
    std::string path = "results";
    Analysis test_analysis(path);
    test_analysis.UpdateParameterDict("parameters.dat");
    EXPECT_EQ(0, 0);
}

TEST(Analysis, InitializeAnalysis) {
    std::string path = "test_gzip_reader";
    Analysis test_analysis(path);
    test_analysis.UpdateParameterDict("parameters.dat");
    test_analysis.InitializeAnalysis();
    EXPECT_EQ(0, 0);
}

TEST(Analysis, FlowAnalysis) {
    std::string path = "test_gzip_reader";
    Analysis test_analysis(path);
    test_analysis.UpdateParameterDict("parameters.dat");
    test_analysis.PerformAnalysis();
    EXPECT_EQ(0, 0);
}

TEST(Analysis, IntegratedTupleCountsIgnorePackedTransverseMomentum) {
    ParameterReader params;
    params.readFromFile("parameters.dat");
    params.setVal("particle_monval", 9999);
    params.setVal("compute_correlation", 1);
    params.setVal("rapidity_distribution", 0);
    singleParticleSpectra spectra(
        params, ".", std::make_shared<RandomUtil::Random>(42));

    // Four tracks at phi = 0 and pT = 0.25 GeV.
    std::vector<double> real(9, 4.0), imag(9, 0.0);
    imag[0] = 1.0;
    std::vector<double> pairs(9), pair_squares(9), triples(9), triple_squares(9);
    std::vector<double> cn4(5), cn4_squares(5), sc(6), sc_squares(6);

    spectra.calculate_two_particle_correlation_charge_base(
        real, imag, real, imag, 0, pairs, pair_squares);
    spectra.calculate_three_particle_correlation(
        real, imag, real, imag, real, imag, 0, triples, triple_squares);
    spectra.calculate_four_particle_correlation_Cn4(
        real, imag, cn4, cn4_squares);
    spectra.calculate_four_particle_correlation_SC(
        real, imag, real, imag, real, imag, real, imag, 0, sc, sc_squares);

    EXPECT_DOUBLE_EQ(pairs[0], 12.0);
    EXPECT_DOUBLE_EQ(triples[0], 24.0);
    EXPECT_DOUBLE_EQ(cn4[0], 24.0);
    EXPECT_DOUBLE_EQ(sc[0], 24.0);
    EXPECT_DOUBLE_EQ(imag[0], 1.0);
}

TEST(Analysis, RapidityDifferentialTripletsMatchEnumeration) {
    ParameterReader params;
    params.readFromFile("parameters.dat");
    params.setVal("particle_monval", 9999);
    params.setVal("compute_correlation", 1);
    params.setVal("n_rap", 3);
    params.setVal("rapidity_dis_min", -1);
    params.setVal("rapidity_dis_max", 1);
    singleParticleSpectra spectra(
        params, ".", std::make_shared<RandomUtil::Random>(42));
    struct Track {
        int bin;
        bool positive;
        double phi;
    };
    const std::vector<Track> tracks = {
        {0, true, -1.1}, {1, true, 0.2}, {1, true, 1.5},
        {1, false, 2.3}, {2, false, -0.6}
    };
    using Grid = std::vector<std::vector<double>>;
    const auto q_vectors = [&tracks](int charge) {
        Grid real(3, std::vector<double>(9)), imag = real;
        for (const auto &track : tracks) {
            if (charge > 0 && !track.positive) continue;
            if (charge < 0 && track.positive) continue;
            real[track.bin][0] += 1.0;
            imag[track.bin][0] += 0.25;
            for (int n = 1; n < 9; n++) {
                real[track.bin][n] += std::cos(n * track.phi);
                imag[track.bin][n] += std::sin(n * track.phi);
            }
        }
        return std::make_pair(real, imag);
    };
    auto inclusive = q_vectors(0);
    for (int charge : {0, 1, -1}) {
        auto first = q_vectors(charge == 0 ? 0 : 1);
        auto second = q_vectors(charge);
        for (int flag : {1, 2}) {
            Grid corr(9, std::vector<double>(3)), err = corr;
            spectra.calculate_three_particle_correlation_deltaeta(
                first.first, first.second, second.first, second.second,
                inclusive.first, inclusive.second, flag,
                charge < 0 ? 2 : 1, corr, err);

            std::vector<double> count(3), c123(3);
            for (const auto &a : tracks) {
                if (charge != 0 && !a.positive) continue;
                for (const auto &b : tracks) {
                    if (&a == &b || (charge > 0 && !b.positive)
                        || (charge < 0 && b.positive)) {
                        continue;
                    }
                    for (const auto &c : tracks) {
                        if (&a == &c || &b == &c) continue;
                        const int bin = a.bin - (flag == 1 ? b.bin : c.bin) + 1;
                        if (bin < 0 || bin >= 3) continue;
                        count[bin]++;
                        c123[bin] += std::cos(a.phi + 2 * b.phi - 3 * c.phi);
                    }
                }
            }
            for (int bin = 0; bin < 3; bin++) {
                EXPECT_NEAR(corr[0][bin], count[bin], 1e-12);
                EXPECT_NEAR(corr[2][bin], c123[bin], 1e-12);
            }
        }
    }
}

}  // namespace
