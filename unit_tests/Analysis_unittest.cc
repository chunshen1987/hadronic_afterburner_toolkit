
#include "Analysis.h"

#include <sys/stat.h>

#include <cmath>
#include <cstdio>
#include <fstream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include "BalanceFunction.h"
#include "gtest/gtest.h"
#include "single_particleSpectra.h"

namespace {

// Data rows of a text output file, skipping comment lines.
std::vector<std::vector<double>> read_rows(const std::string &filename) {
    std::ifstream in(filename.c_str());
    EXPECT_TRUE(in.good()) << filename;
    std::vector<std::vector<double>> rows;
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream fields(line);
        std::vector<double> row;
        double value;
        while (fields >> value) row.push_back(value);
        if (!row.empty()) rows.push_back(row);
    }
    return rows;
}

// One UrQMD text track: 9 leading columns, ityp, 2*I3, charge, lcl#, ncl,
// origin, then t, x, y, z, E, px, py, pz.
std::string urqmd_track(
    double mass, int ityp, int iso3, int charge, double t, double E, double px,
    double py, double pz) {
    const double columns[] = {t, 0.0, 0.0, 0.0, E, px, py, pz, mass};
    std::string line;
    char buffer[64];
    for (double value : columns) {
        std::snprintf(buffer, sizeof(buffer), "%16.8e", value);
        line += buffer;
    }
    std::snprintf(
        buffer, sizeof(buffer), "%11d%3d%3d%9d%5d%4d", ityp, iso3, charge, 0, 0,
        0);
    line += buffer;
    for (int i = 0; i < 8; i++) {
        std::snprintf(buffer, sizeof(buffer), "%16.8e", columns[i]);
        line += buffer;
    }
    return line;
}

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
    particle_info first {}, second {};
    first.monval = second.monval = 211;
    first.pT = second.pT = 0.3;
    first.event_track_id = second.event_track_id = 1;
    first.rap_y = 0.0;
    second.rap_y = 1.0;
    std::vector<particle_info> event_0 {first}, event_1 {second};
    std::vector<std::vector<particle_info> *> events {&event_0, &event_1};
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
    std::vector<std::vector<particle_info> *> single {&event_0};
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

TEST(Analysis, ReconstructedPhiAloneCarriesTheBranchingRatio) {
    // One UrQMD event with a stored phi and a K+ K- pair from a phi decay,
    // both phi at y = 0 and pT = 0.8 GeV. Only the rebuilt phi is rescaled.
    const std::string path = "test_phi_branching_ratio";
    mkdir(path.c_str(), 0755);
    {
        std::ifstream header_in("test_reader_files/particle_list.dat");
        std::ofstream out((path + "/particle_list.dat").c_str());
        std::string line;
        for (int i = 0; i < 17 && std::getline(header_in, line); i++) {
            out << line << '\n';
        }
        const double m_phi = 1.019, m_K = 0.493677, pT = 0.8;
        // Kaons back to back along y in the phi rest frame, then boosted
        // along x: each gets px = pT/2.
        const double p_star = std::sqrt(0.25 * m_phi * m_phi - m_K * m_K);
        const double E_K =
            std::sqrt(m_K * m_K + 0.25 * pT * pT + p_star * p_star);
        out << "           3        8000\n";
        out << "       3       0       0       0       0       0       0"
               "       0\n";
        out << urqmd_track(m_K, 106, 1, 1, 5.0, E_K, 0.5 * pT, p_star, 0.0)
            << '\n';
        out << urqmd_track(m_K, -106, -1, -1, 5.0, E_K, 0.5 * pT, -p_star, 0.0)
            << '\n';
        out << urqmd_track(
            m_phi, 109, 0, 0, 6.0, std::sqrt(m_phi * m_phi + pT * pT), 0.0, pT,
            0.0)
            << '\n';
    }
    ParameterReader params;
    params.readFromFile("parameters.dat");
    params.setVal("read_in_mode", 1);
    params.setVal("read_in_real_mixed_events", 0);
    params.setVal("particle_monval", 333);
    params.setVal("compute_correlation", 0);
    auto random = std::make_shared<RandomUtil::Random>(3);
    auto samples = std::make_shared<particleSamples>(params, path, random);
    samples->read_in_particle_samples_and_filter();
    ASSERT_EQ(samples->get_number_of_events(), 1);
    ASSERT_EQ(samples->get_number_of_particles(0), 2);

    singleParticleSpectra spectra(params, path, random);
    spectra.calculate_Qn_vector_shell(samples);
    spectra.output_spectra_and_Qn_results();

    const double weight = 1.0 / 0.489;
    const double yield = 1.0 + weight;
    const double error = std::sqrt(1.0 + weight * weight);
    const auto integrated =
        read_rows(path + "/particle_333_vndata_y_-0.5_0.5.dat");
    ASSERT_FALSE(integrated.empty());
    EXPECT_NEAR(integrated.front()[1], yield, 1e-7);  // dN/dy, |y| < 0.5
    EXPECT_NEAR(integrated.front()[2], error, 1e-7);
    EXPECT_NEAR(integrated.back()[1], yield, 1e-6);  // row 99: total_N

    double spectrum_total = 0.0;
    for (const auto &row :
         read_rows(path + "/particle_333_vndata_diff_y_-0.5_0.5.dat")) {
        spectrum_total += row.back();
    }
    EXPECT_NEAR(spectrum_total, yield, 1e-6);

    double rapidity_total = 0.0;
    for (const auto &row :
         read_rows(path + "/particle_333_dNdy_pT_0.2_3.dat")) {
        rapidity_total += row[1] * 0.1;  // 41 points on [-2, 2]
    }
    EXPECT_NEAR(rapidity_total, yield, 1e-6);
}

TEST(Analysis, MultiparticleErrorsUseTheRatioDeltaMethod) {
    // Three read_in_mode = 9 events of 4, 5 and 7 pions at eta = 0 with
    // fixed azimuths, so X_e and Y_e can be enumerated directly.
    const std::string path = "test_correlator_errors";
    mkdir(path.c_str(), 0755);
    const std::vector<std::vector<double>> azimuths = {
        {0.1, 0.9, 2.0, 4.4},
        {0.3, 0.35, 1.7, 3.1, 5.0},
        {0.0, 0.2, 0.4, 1.1, 2.5, 3.3, 6.0}};
    std::vector<std::vector<double>> stored_phi;
    {
        std::ofstream out(
            (path + "/particle_list.bin").c_str(), std::ios::binary);
        const float mass = 0.13957f, pT = 0.5f;
        for (const auto &event : azimuths) {
            const int count = event.size();
            out.write(reinterpret_cast<const char *>(&count), sizeof(int));
            std::vector<double> phis;
            for (std::size_t i = 0; i < event.size(); i++) {
                const int pdg = i % 2 == 0 ? 211 : -211;
                const float px = pT * std::cos(event[i]);
                const float py = pT * std::sin(event[i]);
                const float E = std::sqrt(mass * mass + pT * pT);
                const float record[9] = {mass, 0.f, 0.f, 0.f, 0.f,
                                         E,    px,  py,  0.f};
                out.write(reinterpret_cast<const char *>(&pdg), sizeof(int));
                out.write(
                    reinterpret_cast<const char *>(record), sizeof(record));
                phis.push_back(
                    std::atan2(
                        static_cast<double>(py), static_cast<double>(px)));
            }
            stored_phi.push_back(phis);
        }
    }
    ParameterReader params;
    params.readFromFile("parameters.dat");
    params.setVal("read_in_mode", 9);
    params.setVal("read_in_real_mixed_events", 0);
    params.setVal("particle_monval", 9999);
    params.setVal("compute_correlation", 1);
    params.setVal("rapidity_distribution", 0);
    params.setVal("rap_type", 0);
    params.setVal("rap_min", -1);
    params.setVal("rap_max", 1);
    auto random = std::make_shared<RandomUtil::Random>(3);
    auto samples = std::make_shared<particleSamples>(params, path, random);
    samples->read_in_particle_samples_and_filter();
    ASSERT_EQ(samples->get_number_of_events(), 3);

    singleParticleSpectra spectra(params, path, random);
    spectra.calculate_Qn_vector_shell(samples);
    spectra.output_spectra_and_Qn_results();

    // Per event: X = sum of the correlator over distinct ordered tuples,
    // Y = number of tuples. R = sum X / sum Y, error from the delta method.
    struct Moments {
        double X = 0.0, Y = 0.0;
    };
    const auto expect_ratio = [&](const std::vector<Moments> &events,
                                  const std::vector<double> &row) {
        double sum_X = 0.0, sum_Y = 0.0;
        for (const auto &event : events) {
            sum_X += event.X;
            sum_Y += event.Y;
        }
        const double R = sum_X / sum_Y;
        double residuals = 0.0;
        for (const auto &event : events) {
            residuals += (event.X - R * event.Y) * (event.X - R * event.Y);
        }
        EXPECT_NEAR(row[1], R, 1e-7 * std::abs(R) + 1e-12);
        EXPECT_NEAR(
            row[2], std::sqrt(residuals) / sum_Y,
            1e-7 * std::sqrt(residuals) / sum_Y);
    };

    std::vector<Moments> cn4(3), c112(3);
    for (std::size_t e = 0; e < stored_phi.size(); e++) {
        const auto &phi = stored_phi[e];
        const std::size_t M = phi.size();
        for (std::size_t a = 0; a < M; a++) {
            for (std::size_t b = 0; b < M; b++) {
                if (b == a) continue;
                for (std::size_t c = 0; c < M; c++) {
                    if (c == a || c == b) continue;
                    c112[e].X += std::cos(phi[a] + phi[b] - 2.0 * phi[c]);
                    c112[e].Y += 1.0;
                    for (std::size_t d = 0; d < M; d++) {
                        if (d == a || d == b || d == c) continue;
                        cn4[e].X +=
                            std::cos(2.0 * (phi[a] - phi[b] + phi[c] - phi[d]));
                        cn4[e].Y += 1.0;
                    }
                }
            }
        }
    }
    const auto cn4_rows = read_rows(path + "/particle_9999_Cn4_eta_-1_1.dat");
    const auto cmnk_rows = read_rows(path + "/particle_9999_Cmnk_eta_-1_1.dat");
    ASSERT_GT(cn4_rows.size(), 2u);
    ASSERT_GT(cmnk_rows.size(), 1u);
    expect_ratio(cn4, cn4_rows[2]);    // C_2{4} numerator <4>
    expect_ratio(c112, cmnk_rows[1]);  // C_112

    // Row 0 keeps the mean tuple count and its standard error.
    const double counts[] = {24.0, 120.0, 840.0};
    const double mean = (counts[0] + counts[1] + counts[2]) / 3.0;
    double variance = 0.0;
    for (double count : counts) variance += (count - mean) * (count - mean);
    variance /= 3.0;
    EXPECT_NEAR(cn4_rows[0][1], mean, 1e-6);
    EXPECT_NEAR(cn4_rows[0][2], std::sqrt(variance / 3.0), 1e-6);
}

}  // namespace
