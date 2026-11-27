
#include "Analysis.h"

#include <cmath>
#include <string>

#include "single_particleSpectra.h"
#include "gtest/gtest.h"

namespace {

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
            for (int n = 0; n < 9; n++) {
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
