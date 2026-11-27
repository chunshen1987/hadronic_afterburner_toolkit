// Copyright Chun Shen @ 2016
#include "single_particleSpectra.h"

#include <stdlib.h>

#include <algorithm>
#include <cmath>
#include <complex>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <vector>

#include "arsenal.h"
#include "parameters.h"

using AfterburnerUtil::resize2DVector;
using AfterburnerUtil::resize3DVector;
using std::ofstream;
using std::ostringstream;
using std::scientific;
using std::setprecision;
using std::setw;
using std::vector;

namespace {

// The spectra accumulator packs summed pT into Im(Q[0]). Correlation
// kernels use the mathematical Q[0] = multiplicity + 0i.
double correlation_imag(const vector<double> &imag, int harmonic) {
    return harmonic == 0 ? 0.0 : imag[harmonic];
}

// Transverse-energy convention used by ALICE (1603.04775) and
// PHENIX (nucl-ex/0409015).
double transverse_energy(const particle_info &track) {
    const double momentum = std::hypot(track.pT, track.pz);
    if (momentum == 0.0) return 0.0;
    const int baryon_sign = (track.baryon == 1) - (track.baryon == -1);
    return (track.E - baryon_sign * track.mass) * track.pT / momentum;
}

// The ratio uses event sums. Pair count and numerator fluctuate together
// because they come from the same event. Apply the delta method to sqrt(X/Y).
struct PairMoments {
    double sum;
    double squares;
    double cross;
};

void write_pair_moments(
    std::ostream &out, const PairMoments &moments, const PairMoments &pairs,
    int events, bool harmonic_zero) {
    const double mean = moments.sum / events;
    const double second = moments.squares / events;
    double value = 0.0;
    double error = 0.0;
    if (harmonic_zero) {
        value = mean;
        error = sqrt(std::max(0.0, second - mean * mean) / events);
    } else if (moments.sum > 0.0 && pairs.sum > 0.0) {
        const double ratio = moments.sum / pairs.sum;
        value = sqrt(ratio);
        // The residual X - (X/Y)Y has zero mean by construction.
        const double residual_squares = moments.squares
            - 2.0 * ratio * moments.cross + ratio * ratio * pairs.squares;
        error = sqrt(std::max(0.0, residual_squares))
                / (2.0 * value * pairs.sum);
    }
    out << "  " << value << "  " << error << "  " << mean << "  " << second;
}

// Multiparticle correlators are R = sum_e X_e / sum_e Y_e, where X_e is the
// event's tuple numerator and Y_e = local[0] its ordered tuple count. Store
// sum X, sum X^2 and sum X*Y per entry so the error can include the
// covariance between X and Y.
void add_event_moments(
    const vector<double> &local, vector<double> &sum, vector<double> &squares,
    vector<double> *cross) {
    for (std::size_t i = 0; i < local.size(); i++) {
        sum[i] += local[i];
        squares[i] += local[i] * local[i];
        if (cross) (*cross)[i] += local[i] * local[0];
    }
}

// The same for rapidity-binned correlators, local[entry][bin].
void add_event_moments(
    const vector<vector<double>> &local, vector<vector<double>> &sum,
    vector<vector<double>> &squares, vector<vector<double>> *cross) {
    for (std::size_t i = 0; i < local.size(); i++) {
        for (std::size_t bin = 0; bin < local[i].size(); bin++) {
            const double value = local[i][bin];
            sum[i][bin] += value;
            squares[i][bin] += value * value;
            if (cross) (*cross)[i][bin] += value * local[0][bin];
        }
    }
}

// Mean tuple count per event and its standard error.
void tuple_count(
    double sum, double squares, int events, double &mean, double &error) {
    mean = sum / events;
    const double variance = squares / events - mean * mean;
    error = variance > 0.0 ? sqrt(variance / events) : 0.0;
}

// R = sum X / sum Y and its delta-method error under independent events,
//   sigma_R = sqrt(sum (X_e - R Y_e)^2) / sum Y_e.
// Bins without tuples report zero for both, as before.
void tuple_ratio(
    double sum, double squares, double cross, double count_sum,
    double count_squares, double &value, double &error) {
    value = 0.0;
    error = 0.0;
    if (!(count_sum > 0.0)) return;
    value = sum / count_sum;
    const double residual_squares =
        squares - 2.0 * value * cross + value * value * count_squares;
    error = sqrt(std::max(0.0, residual_squares)) / count_sum;
}

}  // namespace

singleParticleSpectra::singleParticleSpectra(
    const ParameterReader &paraRdr, std::string path,
    std::shared_ptr<RandomUtil::Random> ran_gen)
    : path_(path) {
    ran_gen_ptr = ran_gen;

    particle_monval = paraRdr.getVal("particle_monval");
    bWeakFeedDown_ = false;
    if (paraRdr.getVal("resonance_weak_feed_down_flag") == 1) {
        bWeakFeedDown_ = true;
    }

    if (particle_monval == 333) {
        // phi(1020) reconstructed from (K^+, K^-) pairs is rescaled by the
        // branching ratio; phi stored explicitly in the input is not.
        reconst_branching_ratio = 0.489;
    } else {
        reconst_branching_ratio = 1.0;
    }
    reconstructed_weight_ = 1.0 / reconst_branching_ratio;

    if (paraRdr.getVal("ecoOutput", 0) == 1) {
        ecoOutput_ = true;
    } else {
        ecoOutput_ = false;
    }

    order_max = paraRdr.getVal("order_max");
    flag_correlation = paraRdr.getVal("compute_correlation");
    if (particle_monval != 9999 && flag_correlation != 1) order_max = 6;
    // The four-particle kernels require Q_8 for 2*4 and 5+3.
    if (order_max < 1 || (flag_correlation == 1 && order_max < 9)) {
        throw std::invalid_argument(
            "order_max must be positive and at least 9 for multiparticle "
            "correlations");
    }

    Qn_vector_real = vector<double>(order_max, 0.);
    Qn_vector_imag = vector<double>(order_max, 0.);
    Qn_vector_real_err = vector<double>(order_max, 0.);
    Qn_vector_imag_err = vector<double>(order_max, 0.);

    npT = paraRdr.getVal("npT");
    pT_min = paraRdr.getVal("pT_min");
    pT_max = paraRdr.getVal("pT_max");
    if (npT < 1 || !std::isfinite(pT_min) || !std::isfinite(pT_max)
        || pT_min < 0.0 || pT_max <= pT_min) {
        throw std::invalid_argument("Invalid pT spectrum grid");
    }
    dpT = (pT_max - pT_min) / (npT + 1e-15);
    pT_array = vector<double>(npT, 0.);
    pT_mean_array = vector<double>(npT, 0.);
    pT_mean_array_err = vector<double>(npT, 0.);
    reconstructed_count_pT_ = vector<double>(npT, 0.);
    for (int i = 0; i < npT; i++) {
        pT_array[i] = pT_min + dpT * (i + 0.5);
    }

    resize2DVector(Qn_diff_vector_real, order_max, npT, 0);
    resize2DVector(Qn_diff_vector_imag, order_max, npT, 0);
    resize2DVector(Qn_diff_vector_real_err, order_max, npT, 0);
    resize2DVector(Qn_diff_vector_imag_err, order_max, npT, 0);

    total_number_of_events = 0;

    rap_type = paraRdr.getVal("rap_type");
    rap_min = paraRdr.getVal("rap_min");
    rap_max = paraRdr.getVal("rap_max");

    if (particle_monval == 9999)  // use pseudo-rapidity for all charged hadrons
        rap_type = 0;

    if (!std::isfinite(rap_min) || !std::isfinite(rap_max)
        || rap_max <= rap_min) {
        throw std::invalid_argument("Invalid integrated rapidity acceptance");
    }
    rapidity_distribution_flag = paraRdr.getVal("rapidity_distribution");
    if (rapidity_distribution_flag == 1 || flag_correlation == 1) {
        vn_rapidity_dis_pT_min = paraRdr.getVal("vn_rapidity_dis_pT_min");
        vn_rapidity_dis_pT_max = paraRdr.getVal("vn_rapidity_dis_pT_max");
        if (!std::isfinite(vn_rapidity_dis_pT_min)
            || !std::isfinite(vn_rapidity_dis_pT_max)
            || vn_rapidity_dis_pT_min < 0.0
            || vn_rapidity_dis_pT_max <= vn_rapidity_dis_pT_min) {
            throw std::invalid_argument(
                "Invalid correlation reference pT acceptance");
        }
    }
    if (rapidity_distribution_flag == 1) {
        N_rap = paraRdr.getVal("n_rap");
        rapidity_dis_min = paraRdr.getVal("rapidity_dis_min");
        rapidity_dis_max = paraRdr.getVal("rapidity_dis_max");
        if (N_rap < 2 || !std::isfinite(rapidity_dis_min)
            || !std::isfinite(rapidity_dis_max)
            || rapidity_dis_max <= rapidity_dis_min) {
            throw std::invalid_argument("Invalid rapidity distribution grid");
        }
        drap = (rapidity_dis_max - rapidity_dis_min) / (N_rap - 1.);
        rapidity_array = vector<double>(N_rap, 0.);
        dNdy_array = vector<double>(N_rap, 0.);
        dETdy_array = vector<double>(N_rap, 0.);
        dNdy_covMatrix = vector<double>(N_rap * N_rap, 0.);  // N_rap x N_rap
        reconstructed_count_rap_ = vector<double>(N_rap, 0.);
        for (int i = 0; i < N_rap; i++)
            rapidity_array[i] = rapidity_dis_min + i * drap;

        resize2DVector(vn_real_rapidity_dis_array, N_rap, order_max, 0);
        resize2DVector(vn_imag_rapidity_dis_array, N_rap, order_max, 0);
        rapidityPTDistributionFlag_ =
            paraRdr.getVal("rapidityPTDistributionFlag", 0);
        if (rapidityPTDistributionFlag_ != 0) {
            resize3DVector(
                Qn_pTrapdiff_vector, N_rap, npT, 2 * order_max + 2, 0);
        }
    }

    // check dN/dtau distribution
    check_spatial_flag = paraRdr.getVal("check_spatial_dis");
    if (check_spatial_flag == 1) {
        // dN/dtau
        intrinsic_dtau = paraRdr.getVal("intrinsic_dtau");
        N_tau = 50;
        tau_min = 0.6;
        tau_max = 15.0;
        dtau = (tau_max - tau_min) / (N_tau - 1);
        tau_array = new double[N_tau];
        dNdtau_array = new double[N_tau];
        for (int i = 0; i < N_tau; i++) {
            tau_array[i] = tau_min + i * dtau;
            dNdtau_array[i] = 0.0;
        }

        // dN/dx
        intrinsic_dx = paraRdr.getVal("intrinsic_dx");
        N_xpt = 50;
        spatial_x_min = -10.0;
        spatial_x_max = 10.0;
        spatial_r_min = 0.0;
        spatial_r_max = 10.0;
        dspatial_x = (spatial_x_max - spatial_x_min) / (N_xpt - 1);
        dspatial_r = (spatial_r_max - spatial_r_min) / (N_xpt - 1);
        xpt_array = new double[N_xpt];
        ypt_array = new double[N_xpt];
        rpt_array = new double[N_xpt];
        dNdx1_array = new double[N_xpt];
        dNdx2_array = new double[N_xpt];
        dNdr_array = new double[N_xpt];
        for (int i = 0; i < N_xpt; i++) {
            xpt_array[i] = spatial_x_min + i * dspatial_x;
            ypt_array[i] = spatial_x_min + i * dspatial_x;
            rpt_array[i] = spatial_r_min + i * dspatial_r;
            dNdx1_array[i] = 0.0;
            dNdx2_array[i] = 0.0;
            dNdr_array[i] = 0.0;
        }

        // dN/(dtau dx)
        dNdtaudx1_array = new double *[N_tau];
        dNdtaudx2_array = new double *[N_tau];
        for (int i = 0; i < N_tau; i++) {
            dNdtaudx1_array[i] = new double[N_xpt];
            dNdtaudx2_array[i] = new double[N_xpt];
            for (int j = 0; j < N_xpt; j++) {
                dNdtaudx1_array[i][j] = 0.0;
                dNdtaudx2_array[i][j] = 0.0;
            }
        }

        // dN/deta_s
        intrinsic_detas = paraRdr.getVal("intrinsic_detas");
        N_eta_s = 40;
        eta_s_min = -3.0;
        eta_s_max = 3.0;
        deta_s = (eta_s_max - eta_s_min) / (N_eta_s - 1);
        eta_s_array = new double[N_eta_s];
        dNdetas_array = new double[N_eta_s];
        for (int i = 0; i < N_eta_s; i++) {
            eta_s_array[i] = eta_s_min + i * deta_s;
            dNdetas_array[i] = 0.0;
        }
    }

    if (flag_correlation == 1) {
        Qn2_vector = vector<double>(order_max, 0.);
        Qn2_vector_err = vector<double>(order_max, 0.);
        Qn2_vector_cross = vector<double>(order_max, 0.);
        if (N_rap > 0) {
            resize2DVector(QnSP_eta12, order_max, N_rap, 0.);
            resize2DVector(QnSP_eta12_err, order_max, N_rap, 0.);
            resize2DVector(QnSP_eta12_cross, order_max, N_rap, 0.);
        }

        num_corr = 9;
        C_nmk = vector<double>(num_corr, 0.);
        C_nmk_err = vector<double>(num_corr, 0.);
        C_nmk_cross = vector<double>(num_corr, 0.);
        if (N_rap > 0) {
            resize2DVector(C_nmk_eta12, num_corr, N_rap, 0.);
            resize2DVector(C_nmk_eta13, num_corr, N_rap, 0.);
            resize2DVector(C_nmk_eta12_err, num_corr, N_rap, 0.);
            resize2DVector(C_nmk_eta13_err, num_corr, N_rap, 0.);
            resize2DVector(C_nmk_eta12_cross, num_corr, N_rap, 0.);
            resize2DVector(C_nmk_eta13_cross, num_corr, N_rap, 0.);
        }

        flag_charge_dependence = paraRdr.getVal("flag_charge_dependence");
        if (particle_monval != 9999) flag_charge_dependence = 0;
        if (flag_charge_dependence == 1) {
            Cn2_ss = vector<double>(order_max, 0.);
            Cn2_ss_err = vector<double>(order_max, 0.);
            Cn2_ss_cross = vector<double>(order_max, 0.);
            Cn2_os = vector<double>(order_max, 0.);
            Cn2_os_err = vector<double>(order_max, 0.);
            Cn2_os_cross = vector<double>(order_max, 0.);

            if (N_rap > 0) {
                resize2DVector(Cn2_ss_eta12, order_max, N_rap, 0.);
                resize2DVector(Cn2_ss_eta12_err, order_max, N_rap, 0.);
                resize2DVector(Cn2_ss_eta12_cross, order_max, N_rap, 0.);
                resize2DVector(Cn2_os_eta12, order_max, N_rap, 0.);
                resize2DVector(Cn2_os_eta12_err, order_max, N_rap, 0.);
                resize2DVector(Cn2_os_eta12_cross, order_max, N_rap, 0.);
            }

            C_nmk_ss = vector<double>(num_corr, 0.);
            C_nmk_ss_err = vector<double>(num_corr, 0.);
            C_nmk_os = vector<double>(num_corr, 0.);
            C_nmk_os_err = vector<double>(num_corr, 0.);
            C_nmk_ss_13 = vector<double>(num_corr, 0.);
            C_nmk_ss_13_err = vector<double>(num_corr, 0.);
            C_nmk_os_13 = vector<double>(num_corr, 0.);
            C_nmk_os_13_err = vector<double>(num_corr, 0.);
            C_nmk_ss_cross = vector<double>(num_corr, 0.);
            C_nmk_os_cross = vector<double>(num_corr, 0.);
            C_nmk_ss_13_cross = vector<double>(num_corr, 0.);
            C_nmk_os_13_cross = vector<double>(num_corr, 0.);
            if (N_rap > 0) {
                resize2DVector(C_nmk_eta12_ss, num_corr, N_rap, 0.);
                resize2DVector(C_nmk_eta12_os, num_corr, N_rap, 0.);
                resize2DVector(C_nmk_eta13_ss, num_corr, N_rap, 0.);
                resize2DVector(C_nmk_eta13_os, num_corr, N_rap, 0.);
                resize2DVector(C_nmk_eta12_ss_err, num_corr, N_rap, 0.);
                resize2DVector(C_nmk_eta12_os_err, num_corr, N_rap, 0.);
                resize2DVector(C_nmk_eta13_ss_err, num_corr, N_rap, 0.);
                resize2DVector(C_nmk_eta13_os_err, num_corr, N_rap, 0.);
                resize2DVector(C_nmk_eta12_ss_cross, num_corr, N_rap, 0.);
                resize2DVector(C_nmk_eta12_os_cross, num_corr, N_rap, 0.);
                resize2DVector(C_nmk_eta13_ss_cross, num_corr, N_rap, 0.);
                resize2DVector(C_nmk_eta13_os_cross, num_corr, N_rap, 0.);
            }
        }
        SC_num_corr = 6;
        SC_mn = vector<double>(SC_num_corr, 0.);
        SC_mn_err = vector<double>(SC_num_corr, 0.);
        SC_mn_cross = vector<double>(SC_num_corr, 0.);

        num_Cn4 = 5;
        Cn4 = vector<double>(num_Cn4, 0.);
        Cn4_err = vector<double>(num_Cn4, 0.);
        Cn4_cross = vector<double>(num_Cn4, 0.);
    }
}

singleParticleSpectra::~singleParticleSpectra() {
    if (check_spatial_flag == 1) {
        delete[] tau_array;
        delete[] dNdtau_array;
        delete[] xpt_array;
        delete[] ypt_array;
        delete[] rpt_array;
        delete[] dNdx1_array;
        delete[] dNdx2_array;
        delete[] dNdr_array;
        delete[] eta_s_array;
        delete[] dNdetas_array;
        for (int i = 0; i < N_tau; i++) {
            delete[] dNdtaudx1_array[i];
            delete[] dNdtaudx2_array[i];
        }
        delete[] dNdtaudx1_array;
        delete[] dNdtaudx2_array;
    }
}

//! This is a driver function to compute the Qn flow vector
void singleParticleSpectra::calculate_Qn_vector_shell(
    std::shared_ptr<particleSamples> particle_list_in) {
    particle_list = particle_list_in;
    // initialize some temp arrays
    vector<double> event_pT_mean(npT, 0.);
    vector<double> event_pT_mean_err(npT, 0.);
    vector<double> event_Qn_real(order_max, 0.);
    vector<double> event_Qn_real_err(order_max, 0.);
    vector<double> event_Qn_imag(order_max, 0.);
    vector<double> event_Qn_imag_err(order_max, 0.);

    vector<vector<double>> event_Qn_diff_real, event_Qn_diff_real_err;
    vector<vector<double>> event_Qn_diff_imag, event_Qn_diff_imag_err;
    resize2DVector(event_Qn_diff_real, order_max, npT, 0.);
    resize2DVector(event_Qn_diff_real_err, order_max, npT, 0.);
    resize2DVector(event_Qn_diff_imag, order_max, npT, 0.);
    resize2DVector(event_Qn_diff_imag_err, order_max, npT, 0.);

    vector<double> event_Qn_p_real(order_max, 0.);
    vector<double> event_Qn_p_real_err(order_max, 0.);
    vector<double> event_Qn_p_imag(order_max, 0.);
    vector<double> event_Qn_p_imag_err(order_max, 0.);
    vector<double> event_Qn_m_real(order_max, 0.);
    vector<double> event_Qn_m_real_err(order_max, 0.);
    vector<double> event_Qn_m_imag(order_max, 0.);
    vector<double> event_Qn_m_imag_err(order_max, 0.);

    vector<vector<double>> event_Qn_rap_real;
    vector<vector<double>> event_Qn_rap_imag;
    if (rapidity_distribution_flag == 1) {
        resize2DVector(event_Qn_rap_real, N_rap, order_max, 0.);
        resize2DVector(event_Qn_rap_imag, N_rap, order_max, 0.);
    }

    vector<vector<double>> event_Qn_p_rap_real;
    vector<vector<double>> event_Qn_p_rap_imag;
    vector<vector<double>> event_Qn_m_rap_real;
    vector<vector<double>> event_Qn_m_rap_imag;
    if (flag_charge_dependence == 1 && rapidity_distribution_flag == 1) {
        resize2DVector(event_Qn_p_rap_real, N_rap, order_max, 0.);
        resize2DVector(event_Qn_p_rap_imag, N_rap, order_max, 0.);
        resize2DVector(event_Qn_m_rap_real, N_rap, order_max, 0.);
        resize2DVector(event_Qn_m_rap_imag, N_rap, order_max, 0.);
    }

    const bool rapidity_correlations =
        flag_correlation == 1 && rapidity_correlation_output_.enabled;
    const bool rapidity_pt =
        rapidityPTDistributionFlag_ == 1 && rapidity_pt_output_.enabled;
    const bool use_rapidity = rapidity_distribution_flag == 1
        && (rapidity_output_.enabled || covariance_output_.enabled
            || rapidity_correlations || rapidity_pt);
    // start the loop
    int nev = particle_list->get_number_of_events();
    for (int iev = 0; iev < nev; iev++) {
        if (spectra_output_.enabled) {
            calculate_Qn_vector(
                iev, pT_min, pT_max, event_pT_mean, event_pT_mean_err,
                event_Qn_real, event_Qn_real_err, event_Qn_imag,
                event_Qn_imag_err, event_Qn_diff_real, event_Qn_diff_real_err,
                event_Qn_diff_imag, event_Qn_diff_imag_err);
            for (int i = 0; i < npT; i++) {
                pT_mean_array[i] += event_pT_mean[i];
                pT_mean_array_err[i] += event_pT_mean_err[i];
            }
            for (int i = 0; i < order_max; i++) {
                Qn_vector_real[i] += event_Qn_real[i];
                Qn_vector_real_err[i] += event_Qn_real_err[i];
                Qn_vector_imag[i] += event_Qn_imag[i];
                Qn_vector_imag_err[i] += event_Qn_imag_err[i];
                for (int j = 0; j < npT; j++) {
                    Qn_diff_vector_real[i][j] += event_Qn_diff_real[i][j];
                    Qn_diff_vector_real_err[i][j] +=
                        event_Qn_diff_real_err[i][j];
                    Qn_diff_vector_imag[i][j] += event_Qn_diff_imag[i][j];
                    Qn_diff_vector_imag_err[i][j] +=
                        event_Qn_diff_imag_err[i][j];
                }
            }
        }

        if (flag_correlation == 1) {
            accumulate_event_qn(
                iev, 0, vn_rapidity_dis_pT_min, vn_rapidity_dis_pT_max,
                event_Qn_real, event_Qn_real_err, event_Qn_imag,
                event_Qn_imag_err);
            calculate_two_particle_correlation(event_Qn_real, event_Qn_imag);
            // The kernels add one event's values into zeroed buffers.
            vector<double> local_Cnmk(num_corr, 0.), local_SC(SC_num_corr, 0.);
            vector<double> local_Cn4(num_Cn4, 0.), unused(num_corr, 0.);
            calculate_three_particle_correlation(
                event_Qn_real, event_Qn_imag, event_Qn_real, event_Qn_imag,
                event_Qn_real, event_Qn_imag, 0, local_Cnmk, unused);
            add_event_moments(local_Cnmk, C_nmk, C_nmk_err, &C_nmk_cross);
            calculate_four_particle_correlation_SC(
                event_Qn_real, event_Qn_imag, event_Qn_real, event_Qn_imag,
                event_Qn_real, event_Qn_imag, event_Qn_real, event_Qn_imag, 0,
                local_SC, unused);
            add_event_moments(local_SC, SC_mn, SC_mn_err, &SC_mn_cross);
            calculate_four_particle_correlation_Cn4(
                event_Qn_real, event_Qn_imag, local_Cn4, unused);
            add_event_moments(local_Cn4, Cn4, Cn4_err, &Cn4_cross);
            if (flag_charge_dependence == 1) {
                accumulate_event_qn(
                    iev, 1, vn_rapidity_dis_pT_min, vn_rapidity_dis_pT_max,
                    event_Qn_p_real, event_Qn_p_real_err, event_Qn_p_imag,
                    event_Qn_p_imag_err);
                accumulate_event_qn(
                    iev, 2, vn_rapidity_dis_pT_min, vn_rapidity_dis_pT_max,
                    event_Qn_m_real, event_Qn_m_real_err, event_Qn_m_imag,
                    event_Qn_m_imag_err);

                calculate_two_particle_correlation_charge_dep(
                    event_Qn_p_real, event_Qn_p_imag, event_Qn_m_real,
                    event_Qn_m_imag, Cn2_ss, Cn2_ss_err, Cn2_os, Cn2_os_err,
                    &Cn2_ss_cross, &Cn2_os_cross);

                calculate_three_particle_correlation_charge_dep(
                    event_Qn_p_real, event_Qn_p_imag, event_Qn_m_real,
                    event_Qn_m_imag, event_Qn_real, event_Qn_imag, C_nmk_ss,
                    C_nmk_ss_err, C_nmk_os, C_nmk_os_err, C_nmk_ss_13,
                    C_nmk_ss_13_err, C_nmk_os_13, C_nmk_os_13_err,
                    &C_nmk_ss_cross, &C_nmk_os_cross, &C_nmk_ss_13_cross,
                    &C_nmk_os_13_cross);
            }
        }

        if (use_rapidity) {
            calculate_rapidity_distribution(
                iev, event_Qn_rap_real, event_Qn_rap_imag, 0);
            for (int i = 0; i < N_rap; i++) {
                for (int j = 0; j < order_max; j++) {
                    vn_real_rapidity_dis_array[i][j] +=
                        (event_Qn_rap_real[i][j]);
                    vn_imag_rapidity_dis_array[i][j] +=
                        (event_Qn_rap_imag[i][j]);
                }
            }
            if (rapidity_pt) {
                calculateRapidityPTDistribution(iev);
            }
            if (rapidity_correlations) {
                calculate_two_particle_correlation_deltaeta(
                    event_Qn_rap_real, event_Qn_rap_imag);
                vector<vector<double>> local_eta(
                    num_corr, vector<double>(N_rap, 0.));
                vector<vector<double>> unused_eta(
                    num_corr, vector<double>(N_rap, 0.));
                calculate_three_particle_correlation_deltaeta(
                    event_Qn_rap_real, event_Qn_rap_imag, event_Qn_rap_real,
                    event_Qn_rap_imag, event_Qn_rap_real, event_Qn_rap_imag, 1,
                    1, local_eta, unused_eta);
                add_event_moments(
                    local_eta, C_nmk_eta12, C_nmk_eta12_err,
                    &C_nmk_eta12_cross);
                for (auto &row : local_eta) {
                    std::fill(row.begin(), row.end(), 0.0);
                }
                calculate_three_particle_correlation_deltaeta(
                    event_Qn_rap_real, event_Qn_rap_imag, event_Qn_rap_real,
                    event_Qn_rap_imag, event_Qn_rap_real, event_Qn_rap_imag, 2,
                    1, local_eta, unused_eta);
                add_event_moments(
                    local_eta, C_nmk_eta13, C_nmk_eta13_err,
                    &C_nmk_eta13_cross);
                if (flag_charge_dependence == 1) {
                    calculate_rapidity_distribution(
                        iev, event_Qn_p_rap_real, event_Qn_p_rap_imag, 1);
                    calculate_rapidity_distribution(
                        iev, event_Qn_m_rap_real, event_Qn_m_rap_imag, 2);

                    calculate_two_particle_correlation_deltaeta_chdep(
                        event_Qn_p_rap_real, event_Qn_p_rap_imag,
                        event_Qn_m_rap_real, event_Qn_m_rap_imag, Cn2_ss_eta12,
                        Cn2_ss_eta12_err, Cn2_os_eta12, Cn2_os_eta12_err,
                        &Cn2_ss_eta12_cross, &Cn2_os_eta12_cross);

                    calculate_three_particle_correlation_deltaeta_chdep(
                        event_Qn_p_rap_real, event_Qn_p_rap_imag,
                        event_Qn_m_rap_real, event_Qn_m_rap_imag,
                        event_Qn_rap_real, event_Qn_rap_imag, C_nmk_eta12_ss,
                        C_nmk_eta12_ss_err, C_nmk_eta12_os, C_nmk_eta12_os_err,
                        C_nmk_eta13_ss, C_nmk_eta13_ss_err, C_nmk_eta13_os,
                        C_nmk_eta13_os_err, &C_nmk_eta12_ss_cross,
                        &C_nmk_eta12_os_cross, &C_nmk_eta13_ss_cross,
                        &C_nmk_eta13_os_cross);
                }
            }
        }

        if (check_spatial_flag == 1) check_dNdSV(iev);
    }

    total_number_of_events += nev;
}

std::string singleParticleSpectra::output_suffix(
    const OutputSelection &selection) const {
    return selection.suffix + (bWeakFeedDown_ ? "_weakFD.dat" : ".dat");
}

void singleParticleSpectra::plan_outputs(
    std::map<std::string, std::set<std::string>> &claimed) {
    ostringstream prefix, window, pt_grid, rapidity_grid, reference_pt;
    prefix << particle_monval << (bWeakFeedDown_ ? "_weakFD" : "");
    window << (rap_type == 0 ? "_eta_" : "_y_") << rap_min << "_" << rap_max;
    pt_grid << "_pT_" << pT_min << "_" << pT_max << "_" << npT;
    const auto claim = [&](OutputSelection &output, const std::string &key,
                           const std::string &signature,
                           const std::string &suffix) {
        auto &definitions = claimed[prefix.str() + key];
        const bool first = definitions.empty();
        output.enabled = definitions.insert(signature).second;
        if (output.enabled && !first) output.suffix = suffix;
    };

    claim(spectra_output_, "_vndata" + window.str(), pt_grid.str(),
          pt_grid.str());
    if (check_spatial_flag == 1) {
        claim(spatial_output_, "_spatial", window.str(), window.str());
    }
    if (rapidity_distribution_flag != 1) return;

    rapidity_grid << "_grid_" << N_rap << "_" << std::setprecision(17)
                  << rapidity_dis_min << "_" << rapidity_dis_max;
    reference_pt << "_pT_" << vn_rapidity_dis_pT_min << "_"
                 << vn_rapidity_dis_pT_max;
    const std::string coordinate = rap_type == 0 ? "_eta" : "_y";
    claim(rapidity_output_, "_density" + coordinate + reference_pt.str(),
          rapidity_grid.str(), rapidity_grid.str());
    // The count covariance has no reference pT or integrated rapidity cut.
    claim(covariance_output_, "_covariance" + coordinate,
          rapidity_grid.str(), rapidity_grid.str());
    if (rapidityPTDistributionFlag_ == 1) {
        const auto grid = rapidity_grid.str() + pt_grid.str();
        claim(rapidity_pt_output_, "_rapidity_pt" + coordinate, grid, grid);
    }
    if (flag_correlation == 1) {
        // Rapidity-difference correlations use the full rapidity grid.
        const auto grid = coordinate + rapidity_grid.str();
        claim(rapidity_correlation_output_,
              "_rapidity_corr" + reference_pt.str(), grid, grid);
    }
}

void singleParticleSpectra::output_spectra_and_Qn_results() {
    if (total_number_of_events == 0) return;
    if (spectra_output_.enabled) output_Qn_vectors();
    if (rapidity_distribution_flag == 1) {
        output_rapidity_distribution();
        if (rapidityPTDistributionFlag_ == 1 && rapidity_pt_output_.enabled) {
            outputRapidityPTDistribution();
        }
        if (flag_correlation == 1 && rapidity_correlation_output_.enabled) {
            output_two_particle_correlation_rap();
            output_three_particle_correlation_rap();
        }
    }
    if (check_spatial_flag == 1 && spatial_output_.enabled) output_dNdSV();
    if (flag_correlation == 1) {
        output_two_particle_correlation();
        output_three_particle_correlation();
        output_four_particle_Cn4_correlation();
        output_four_particle_SC_correlation();
    }
}

//! this function computes the pT-integrated and pT-differential Qn vector
//! within a given rapidity region in one event
void singleParticleSpectra::calculate_Qn_vector(
    int event_id, double pT_min_selected, double pT_max_selected,
    vector<double> &event_pT_mean, vector<double> &event_pT_mean_err,
    vector<double> &event_Qn_real, vector<double> &event_Qn_real_err,
    vector<double> &event_Qn_imag, vector<double> &event_Qn_imag_err,
    vector<vector<double>> &event_Qn_diff_real,
    vector<vector<double>> &event_Qn_diff_real_err,
    vector<vector<double>> &event_Qn_diff_imag,
    vector<vector<double>> &event_Qn_diff_imag_err) {
    DifferentialMoments differential(
        event_Qn_diff_real, event_Qn_diff_real_err,
        event_Qn_diff_imag, event_Qn_diff_imag_err);
    differential.pt = &event_pT_mean;
    differential.pt_squares = &event_pT_mean_err;
    accumulate_event_qn(
        event_id, 0, pT_min_selected, pT_max_selected,
        event_Qn_real, event_Qn_real_err, event_Qn_imag, event_Qn_imag_err,
        &differential);
}

void singleParticleSpectra::calculate_Qn_vector_positive_charge(
    int event_id,
    vector<double> &event_Qn_real, vector<double> &event_Qn_real_err,
    vector<double> &event_Qn_imag, vector<double> &event_Qn_imag_err,
    vector<vector<double>> &event_Qn_diff_real,
    vector<vector<double>> &event_Qn_diff_real_err,
    vector<vector<double>> &event_Qn_diff_imag,
    vector<vector<double>> &event_Qn_diff_imag_err) {
    DifferentialMoments differential(
        event_Qn_diff_real, event_Qn_diff_real_err,
        event_Qn_diff_imag, event_Qn_diff_imag_err);
    accumulate_event_qn(
        event_id, 1, vn_rapidity_dis_pT_min, vn_rapidity_dis_pT_max,
        event_Qn_real, event_Qn_real_err, event_Qn_imag, event_Qn_imag_err,
        &differential);
}

void singleParticleSpectra::calculate_Qn_vector_negative_charge(
    int event_id,
    vector<double> &event_Qn_real, vector<double> &event_Qn_real_err,
    vector<double> &event_Qn_imag, vector<double> &event_Qn_imag_err,
    vector<vector<double>> &event_Qn_diff_real,
    vector<vector<double>> &event_Qn_diff_real_err,
    vector<vector<double>> &event_Qn_diff_imag,
    vector<vector<double>> &event_Qn_diff_imag_err) {
    DifferentialMoments differential(
        event_Qn_diff_real, event_Qn_diff_real_err,
        event_Qn_diff_imag, event_Qn_diff_imag_err);
    accumulate_event_qn(
        event_id, 2, vn_rapidity_dis_pT_min, vn_rapidity_dis_pT_max,
        event_Qn_real, event_Qn_real_err, event_Qn_imag, event_Qn_imag_err,
        &differential);
}

int singleParticleSpectra::selected_track_count(int event_id, int charge) {
    if (charge == 1) {
        return particle_list->get_number_of_positive_particles(event_id);
    }
    if (charge == 2) {
        return particle_list->get_number_of_negative_particles(event_id);
    }
    return particle_list->get_number_of_particles(event_id);
}

particle_info singleParticleSpectra::selected_track(
    int event_id, int index, int charge) {
    if (charge == 1) {
        return particle_list->get_positive_particle(event_id, index);
    }
    if (charge == 2) {
        return particle_list->get_negative_particle(event_id, index);
    }
    return particle_list->get_particle(event_id, index);
}

void singleParticleSpectra::accumulate_event_qn(
    int event_id, int charge, double pT_min_selected, double pT_max_selected,
    vector<double> &event_Qn_real, vector<double> &event_Qn_real_err,
    vector<double> &event_Qn_imag, vector<double> &event_Qn_imag_err,
    DifferentialMoments *differential) {
    std::fill(event_Qn_real.begin(), event_Qn_real.end(), 0.0);
    std::fill(event_Qn_real_err.begin(), event_Qn_real_err.end(), 0.0);
    std::fill(event_Qn_imag.begin(), event_Qn_imag.end(), 0.0);
    std::fill(event_Qn_imag_err.begin(), event_Qn_imag_err.end(), 0.0);
    if (differential) {
        const auto matrices = {&differential->real, &differential->real_squares,
                               &differential->imag, &differential->imag_squares};
        for (auto *moments : matrices) {
            for (auto &row : *moments) std::fill(row.begin(), row.end(), 0.0);
        }
        if (differential->pt) {
            std::fill(differential->pt->begin(), differential->pt->end(), 0.0);
            std::fill(differential->pt_squares->begin(),
                      differential->pt_squares->end(), 0.0);
        }
    }

    const int number_of_particles = selected_track_count(event_id, charge);
    for (int i = 0; i < number_of_particles; i++) {
        const auto track = selected_track(event_id, i, charge);
        const double rap = rap_type == 0 ? track.rap_eta : track.rap_y;
        const double pt = track.pT;
        if (!(rap > rap_min && rap < rap_max
              && pt > pT_min_selected && pt < pT_max_selected)) continue;

        int bin = -1;
        if (differential && pt >= pT_min && pt < pT_max) {
            const int index = static_cast<int>((pt - pT_min) / dpT);
            if (index < npT) bin = index;
        }
        if (track.reconstructed && differential && differential->pt) {
            reconstructed_count_ += 1.0;
            if (bin >= 0) reconstructed_count_pT_[bin] += 1.0;
        }
        event_Qn_real[0] += 1.0;
        event_Qn_real_err[0] += 1.0;
        event_Qn_imag[0] += pt;  // spectra use this slot for mean pT
        event_Qn_imag_err[0] += pt * pt;
        if (bin >= 0) {
            differential->real[0][bin] += 1.0;
            differential->real_squares[0][bin] += 1.0;
            if (differential->pt) {
                (*differential->pt)[bin] += pt;
                (*differential->pt_squares)[bin] += pt * pt;
            }
        }
        for (int n = 1; n < order_max; n++) {
            const double real = cos(n * track.phi_p);
            const double imag = sin(n * track.phi_p);
            event_Qn_real[n] += real;
            event_Qn_real_err[n] += real * real;
            event_Qn_imag[n] += imag;
            event_Qn_imag_err[n] += imag * imag;
            if (bin >= 0) {
                differential->real[n][bin] += real;
                differential->real_squares[n][bin] += real * real;
                differential->imag[n][bin] += imag;
                differential->imag_squares[n][bin] += imag * imag;
            }
        }
    }
}

//! This function outputs the event averaged particle pT-spectra
//! and flow coefficients
void singleParticleSpectra::output_Qn_vectors() {
    double drapidity = rap_max - rap_min;
    // pT-integrated flow
    ostringstream filename;
    filename << path_ << "/particle_" << particle_monval << "_vndata";
    if (rap_type == 0) {
        filename << "_eta_";
    } else {
        filename << "_y_";
    }
    filename << rap_min << "_" << rap_max;
    filename << output_suffix(spectra_output_);
    ofstream output(filename.str().c_str());
    output << "# n  Qn_real  Qn_real_err  Qn_imag  Qn_imag_err" << '\n';

    // Each reconstructed track carries weight w = 1/branching ratio; stored
    // tracks carry one. The Poisson error uses the sum of squared weights.
    const double w = reconstructed_weight_;
    const double total_count = Qn_vector_real[0];
    const double total_N = total_count + (w - 1.0) * reconstructed_count_;
    const double total_N_sq =
        total_count + (w * w - 1.0) * reconstructed_count_;
    const double yield_norm = total_number_of_events * drapidity;
    const double dN_ev_avg = total_N / yield_norm;
    const double dN_ev_avg_err = sqrt(total_N_sq) / yield_norm;
    const double pT_avg =
        total_count > 0.0 ? Qn_vector_imag[0] / total_count : 0.0;
    const double pT_avg_err =
        total_count > 0.0 ? sqrt(
                                std::max(
                                    0.0, Qn_vector_imag_err[0] / total_count
                                             - pT_avg * pT_avg)
                                / total_count)
                          : 0.0;
    output << scientific << setw(18) << setprecision(8) << 0 << "   "
           << dN_ev_avg << "   " << dN_ev_avg_err << "   " << pT_avg << "   "
           << pT_avg_err << '\n';
    for (int iorder = 1; iorder < order_max; iorder++) {
        double vn_evavg_real = 0.0;
        double vn_evavg_imag = 0.0;
        double vn_real_err = 0.0;
        double vn_imag_err = 0.0;
        if (Qn_vector_real[0] > 0) {
            vn_evavg_real = Qn_vector_real[iorder] / Qn_vector_real[0];
            vn_evavg_imag = Qn_vector_imag[iorder] / Qn_vector_real[0];
            vn_real_err =
                (sqrt(
                     Qn_vector_real_err[iorder] / Qn_vector_real[0]
                     - vn_evavg_real * vn_evavg_real)
                 / sqrt(Qn_vector_real[0]));
            vn_imag_err =
                (sqrt(
                     Qn_vector_imag_err[iorder] / Qn_vector_real[0]
                     - vn_evavg_imag * vn_evavg_imag)
                 / sqrt(Qn_vector_real[0]));
            if (std::isnan(vn_real_err)) {
                vn_real_err = 0.0;
            }
            if (std::isnan(vn_imag_err)) {
                vn_imag_err = 0.0;
            }
        }
        output << scientific << setw(18) << setprecision(8) << iorder << "   "
               << vn_evavg_real << "   " << vn_real_err << "   "
               << vn_evavg_imag << "   " << vn_imag_err << '\n';
    }
    // output total number of particles in the last row
    // this quantities is useful when one wants to reconst the Qn vectors
    output << scientific << setw(18) << setprecision(8) << 99 << "   "
           << total_N << "   " << total_number_of_events << "   " << 0.0
           << "   " << 0.0 << '\n';
    output.close();

    // pT-differential flow
    ostringstream filename_diff;
    filename_diff << path_ << "/particle_" << particle_monval << "_vndata_diff";
    if (rap_type == 0) {
        filename_diff << "_eta_";
    } else {
        filename_diff << "_y_";
    }
    filename_diff << rap_min << "_" << rap_max;
    filename_diff << output_suffix(spectra_output_);
    ofstream output_diff(filename_diff.str().c_str());
    if (ecoOutput_) {
        output_diff << "# pT  EdN/d^3p  Qn_real(pT)  Qn_imag(pT) (n=1-"
                    << order_max - 1 << ")" << '\n';
    } else {
        output_diff << "# pT  pT_err EdN/d^3p  dN_err  "
                    << "Qn_real(pT)  Qn_real_err  Qn_imag(pT)  Qn_imag_err"
                    << "(n=1-" << order_max - 1 << ")" << '\n';
    }

    for (int ipT = 0; ipT < npT; ipT++) {
        const double bin_count = Qn_diff_vector_real[0][ipT];
        const double bin_reconstructed = reconstructed_count_pT_[ipT];
        const double total_NpT = bin_count + (w - 1.0) * bin_reconstructed;
        const double dNpT_ev_avg = total_NpT / total_number_of_events;
        const double dNpT_ev_avg_err =
            sqrt(bin_count + (w * w - 1.0) * bin_reconstructed)
            / total_number_of_events;
        double mean_pT, mean_pT_err;
        if (dNpT_ev_avg > 0.) {
            mean_pT = pT_mean_array[ipT] / Qn_diff_vector_real[0][ipT];
            const double variance =
                pT_mean_array_err[ipT] / Qn_diff_vector_real[0][ipT]
                - mean_pT * mean_pT;
            mean_pT_err =
                sqrt(std::max(0.0, variance) / Qn_diff_vector_real[0][ipT]);
        } else {
            mean_pT = pT_array[ipT];
            mean_pT_err = 0.0;
        }
        if (ecoOutput_) {
            output_diff << scientific << setw(18) << setprecision(8) << mean_pT
                        << "   "
                        << dNpT_ev_avg / mean_pT / dpT / (2 * M_PI) / drapidity
                        << "   ";
        } else {
            output_diff << scientific << setw(18) << setprecision(8) << mean_pT
                        << "   " << mean_pT_err << "   "
                        << dNpT_ev_avg / mean_pT / dpT / (2 * M_PI) / drapidity
                        << "   "
                        << dNpT_ev_avg_err / mean_pT / dpT / (2 * M_PI)
                               / drapidity;
        }
        for (int iorder = 1; iorder < order_max; iorder++) {
            if (dNpT_ev_avg > 0.) {
                double vn_evavg_real =
                    (Qn_diff_vector_real[iorder][ipT]
                     / Qn_diff_vector_real[0][ipT]);
                double vn_evavg_imag =
                    (Qn_diff_vector_imag[iorder][ipT]
                     / Qn_diff_vector_real[0][ipT]);
                double vn_evavg_real_err =
                    (sqrt(
                         Qn_diff_vector_real_err[iorder][ipT]
                             / Qn_diff_vector_real[0][ipT]
                         - vn_evavg_real * vn_evavg_real)
                     / sqrt(Qn_diff_vector_real[0][ipT]));
                double vn_evavg_imag_err =
                    (sqrt(
                         Qn_diff_vector_imag_err[iorder][ipT]
                             / Qn_diff_vector_real[0][ipT]
                         - vn_evavg_imag * vn_evavg_imag)
                     / sqrt(Qn_diff_vector_real[0][ipT]));
                if (std::isnan(vn_evavg_real_err)) {
                    vn_evavg_real_err = 0.0;
                }
                if (std::isnan(vn_evavg_imag_err)) {
                    vn_evavg_imag_err = 0.0;
                }
                if (ecoOutput_) {
                    output_diff << scientific << setw(18) << setprecision(8)
                                << vn_evavg_real << "   " << vn_evavg_imag
                                << "   ";
                } else {
                    output_diff << scientific << setw(18) << setprecision(8)
                                << vn_evavg_real << "   " << vn_evavg_real_err
                                << "   " << vn_evavg_imag << "   "
                                << vn_evavg_imag_err << "   ";
                }
            } else {
                if (ecoOutput_) {
                    output_diff << scientific << setw(18) << setprecision(8)
                                << 0.0e0 << "   " << 0.0e0 << "   ";
                } else {
                    output_diff << scientific << setw(18) << setprecision(8)
                                << 0.0e0 << "   " << 0.0e0 << "   " << 0.0e0
                                << "   " << 0.0e0 << "   ";
                }
            }
        }
        // output total number of particles in the last column
        // this quantities is useful when one wants to reconst the Qn vectors
        output_diff << scientific << setw(18) << setprecision(8) << total_NpT;
        output_diff << '\n';
    }
    output_diff.close();
}

//! This function computes the 2-particle correlation for Qn vectors
//! within one event with same sign pairs
void singleParticleSpectra::calculate_two_particle_correlation_charge_dep(
    vector<double> &event_Qn_p_real, vector<double> &event_Qn_p_imag,
    vector<double> &event_Qn_m_real, vector<double> &event_Qn_m_imag,
    vector<double> &corr_ss, vector<double> &corr_ss_err,
    vector<double> &corr_os, vector<double> &corr_os_err,
    vector<double> *cross_ss, vector<double> *cross_os) {
    vector<double> temp_Cn2_ss_1(order_max, 0.);
    vector<double> temp_Cn2_ss_2(order_max, 0.);
    vector<double> temp_Cn2_ss_1_err(order_max, 0.);
    vector<double> temp_Cn2_ss_2_err(order_max, 0.);
    calculate_two_particle_correlation_charge_base(
        event_Qn_p_real, event_Qn_p_imag, event_Qn_p_real, event_Qn_p_imag, 0,
        temp_Cn2_ss_1, temp_Cn2_ss_1_err);
    calculate_two_particle_correlation_charge_base(
        event_Qn_m_real, event_Qn_m_imag, event_Qn_m_real, event_Qn_m_imag, 0,
        temp_Cn2_ss_2, temp_Cn2_ss_2_err);
    for (int i = 0; i < order_max; i++) {
        double local_Cn2_ss = temp_Cn2_ss_1[i] + temp_Cn2_ss_2[i];
        corr_ss[i] += local_Cn2_ss;
        corr_ss_err[i] += local_Cn2_ss * local_Cn2_ss;
        if (cross_ss) {
            (*cross_ss)[i] += local_Cn2_ss
                * (temp_Cn2_ss_1[0] + temp_Cn2_ss_2[0]);
        }
    }

    vector<double> temp_Cn2_os_1(order_max, 0.);
    vector<double> temp_Cn2_os_2(order_max, 0.);
    vector<double> temp_Cn2_os_1_err(order_max, 0.);
    vector<double> temp_Cn2_os_2_err(order_max, 0.);
    calculate_two_particle_correlation_charge_base(
        event_Qn_p_real, event_Qn_p_imag, event_Qn_m_real, event_Qn_m_imag, 1,
        temp_Cn2_os_1, temp_Cn2_os_1_err);
    calculate_two_particle_correlation_charge_base(
        event_Qn_m_real, event_Qn_m_imag, event_Qn_p_real, event_Qn_p_imag, 1,
        temp_Cn2_os_2, temp_Cn2_os_2_err);
    for (int i = 0; i < order_max; i++) {
        double local_Cn2_os = temp_Cn2_os_1[i] + temp_Cn2_os_2[i];
        corr_os[i] += local_Cn2_os;
        corr_os_err[i] += local_Cn2_os * local_Cn2_os;
        if (cross_os) {
            (*cross_os)[i] += local_Cn2_os
                * (temp_Cn2_os_1[0] + temp_Cn2_os_2[0]);
        }
    }
}

//! This function computes the 2-particle correlation for Qn vectors
//! within one event with particle dependence
//!     Real(Qn*conj(Qn)) for n = 0, 1, ... , order_max
//! self correlation is subtracted when flag == 0 (full overlap)
void singleParticleSpectra::calculate_two_particle_correlation_charge_base(
    vector<double> &event_Q1_real, vector<double> &event_Q1_imag,
    vector<double> &event_Q2_real, vector<double> &event_Q2_imag, int flag,
    vector<double> &corr, vector<double> &corr_err) {
    for (int i = 0; i < order_max; i++) {
        double Q2_local =
            event_Q1_real[i] * event_Q2_real[i]
            + correlation_imag(event_Q1_imag, i)
                  * correlation_imag(event_Q2_imag, i);
        if (flag == 0) {
            // subtract the self correlation
            Q2_local -= event_Q1_real[0];
        }
        corr[i] += Q2_local;
        corr_err[i] += Q2_local * Q2_local;
    }
}

//! This function computes the 2-particle correlation for Qn vectors
//! within one event
//!     Real(Qn*conj(Qn)) for n = 0, 1, ... , order_max
//! self correlation is subtracted assuming full overlap
void singleParticleSpectra::calculate_two_particle_correlation(
    vector<double> &event_Qn_real, vector<double> &event_Qn_imag,
    vector<vector<double>> &, vector<vector<double>> &) {
    calculate_two_particle_correlation(event_Qn_real, event_Qn_imag);
}

void singleParticleSpectra::calculate_two_particle_correlation(
    const vector<double> &event_Qn_real,
    const vector<double> &event_Qn_imag) {
    const double multiplicity = event_Qn_real[0];
    const double pairs = multiplicity * (multiplicity - 1.0);
    for (int i = 0; i < order_max; i++) {
        const double imag = correlation_imag(event_Qn_imag, i);
        const double Q2_local =
            event_Qn_real[i] * event_Qn_real[i] + imag * imag - multiplicity;
        Qn2_vector[i] += Q2_local;
        Qn2_vector_err[i] += Q2_local * Q2_local;
        Qn2_vector_cross[i] += Q2_local * pairs;
    }
}

//! This function computes the 2-particle correlation for Qn vectors
//! as a function of \delta \eta within one event
//!     Real(Qn(eta_1)*conj(Qn(eta_2))) for n = 0, 1, ... , order_max
//! self correlation is subtracted when eta_1 = eta_2
void singleParticleSpectra::calculate_two_particle_correlation_deltaeta(
    vector<vector<double>> &event_Qn_diff_real,
    vector<vector<double>> &event_Qn_diff_imag) {
    vector<double> temp_corr(N_rap);
    vector<double> pair_count(N_rap);
    for (int i = 0; i < order_max; i++) {
        for (int j = 0; j < N_rap; j++) {
            temp_corr[j] = 0.0;
        }
        for (int ii = 0; ii < N_rap; ii++) {
            for (int jj = 0; jj < N_rap; jj++) {
                int rap_idx =
                    ((ii - jj) - static_cast<int>(rapidity_dis_min / drap));
                if (rap_idx >= 0 && rap_idx < N_rap) {
                    double QnSP_rap_local =
                        (event_Qn_diff_real[ii][i] * event_Qn_diff_real[jj][i]
                         + correlation_imag(event_Qn_diff_imag[ii], i)
                               * correlation_imag(event_Qn_diff_imag[jj], i));
                    if (ii == jj) {
                        QnSP_rap_local -= event_Qn_diff_real[ii][0];
                    }
                    temp_corr[rap_idx] += QnSP_rap_local;
                }
            }
        }
        for (int j = 0; j < N_rap; j++) {
            QnSP_eta12[i][j] += temp_corr[j];
            QnSP_eta12_err[i][j] += temp_corr[j] * temp_corr[j];
            if (i == 0) pair_count[j] = temp_corr[j];
            QnSP_eta12_cross[i][j] += temp_corr[j] * pair_count[j];
        }
    }
}

//! This function computes the 2-particle correlation for Qn vectors
//! as a function of \delta \eta within one event with charge depenedence
//!     Real(Qn(eta_1)*conj(Qn(eta_2))) for n = 0, 1, ... , order_max
//! self correlation is subtracted when eta_1 = eta_2 and flag == 0
void singleParticleSpectra::calculate_two_particle_correlation_deltaeta_chdep(
    vector<vector<double>> &event_Qn_p_rap_real,
    vector<vector<double>> &event_Qn_p_rap_imag,
    vector<vector<double>> &event_Qn_m_rap_real,
    vector<vector<double>> &event_Qn_m_rap_imag,
    vector<vector<double>> &C_n2_ss_eta12,
    vector<vector<double>> &C_n2_ss_eta12_err,
    vector<vector<double>> &C_n2_os_eta12,
    vector<vector<double>> &C_n2_os_eta12_err,
    vector<vector<double>> *cross_ss, vector<vector<double>> *cross_os) {
    vector<vector<double>> temp_Cn2_ss_1(order_max, vector<double>(N_rap, 0.));
    vector<vector<double>> temp_Cn2_ss_2(order_max, vector<double>(N_rap, 0.));
    vector<vector<double>> temp_Cn2_os_1(order_max, vector<double>(N_rap, 0.));
    vector<vector<double>> temp_Cn2_os_2(order_max, vector<double>(N_rap, 0.));
    vector<vector<double>> temp_Cn2_err(order_max, vector<double>(N_rap, 0.));

    calculate_two_particle_correlation_deltaeta_chdep_base(
        event_Qn_p_rap_real, event_Qn_p_rap_imag, event_Qn_p_rap_real,
        event_Qn_p_rap_imag, 0, temp_Cn2_ss_1, temp_Cn2_err);
    calculate_two_particle_correlation_deltaeta_chdep_base(
        event_Qn_m_rap_real, event_Qn_m_rap_imag, event_Qn_m_rap_real,
        event_Qn_m_rap_imag, 0, temp_Cn2_ss_2, temp_Cn2_err);
    calculate_two_particle_correlation_deltaeta_chdep_base(
        event_Qn_p_rap_real, event_Qn_p_rap_imag, event_Qn_m_rap_real,
        event_Qn_m_rap_imag, 1, temp_Cn2_os_1, temp_Cn2_err);
    calculate_two_particle_correlation_deltaeta_chdep_base(
        event_Qn_m_rap_real, event_Qn_m_rap_imag, event_Qn_p_rap_real,
        event_Qn_p_rap_imag, 1, temp_Cn2_os_2, temp_Cn2_err);

    for (int i = 0; i < order_max; i++) {
        for (int j = 0; j < N_rap; j++) {
            double temp_ss = temp_Cn2_ss_1[i][j] + temp_Cn2_ss_2[i][j];
            C_n2_ss_eta12[i][j] += temp_ss;
            C_n2_ss_eta12_err[i][j] += temp_ss * temp_ss;
            if (cross_ss) {
                (*cross_ss)[i][j] += temp_ss
                    * (temp_Cn2_ss_1[0][j] + temp_Cn2_ss_2[0][j]);
            }
            double temp_os = temp_Cn2_os_1[i][j] + temp_Cn2_os_2[i][j];
            C_n2_os_eta12[i][j] += temp_os;
            C_n2_os_eta12_err[i][j] += temp_os * temp_os;
            if (cross_os) {
                (*cross_os)[i][j] += temp_os
                    * (temp_Cn2_os_1[0][j] + temp_Cn2_os_2[0][j]);
            }
        }
    }
}

//! This function computes the 2-particle correlation for Qn vectors
//! as a function of \delta \eta within one event with charge depenedence
//!     Real(Qn(eta_1)*conj(Qn(eta_2))) for n = 0, 1, ... , order_max
//! self correlation is subtracted when eta_1 = eta_2 and flag == 0
void
singleParticleSpectra::calculate_two_particle_correlation_deltaeta_chdep_base(
    vector<vector<double>> &event_Q1_diff_real,
    vector<vector<double>> &event_Q1_diff_imag,
    vector<vector<double>> &event_Q2_diff_real,
    vector<vector<double>> &event_Q2_diff_imag, int flag,
    vector<vector<double>> &corr, vector<vector<double>> &corr_err) {
    vector<double> temp_corr(N_rap, 0.);
    for (int i = 0; i < order_max; i++) {
        for (int j = 0; j < N_rap; j++) {
            temp_corr[j] = 0.0;
        }
        for (int ii = 0; ii < N_rap; ii++) {
            for (int jj = 0; jj < N_rap; jj++) {
                int rap_idx =
                    ((ii - jj) - static_cast<int>(rapidity_dis_min / drap));
                if (rap_idx >= 0 && rap_idx < N_rap) {
                    double QnSP_rap_local =
                        (event_Q1_diff_real[ii][i] * event_Q2_diff_real[jj][i]
                         + correlation_imag(event_Q1_diff_imag[ii], i)
                               * correlation_imag(event_Q2_diff_imag[jj], i));
                    if (ii == jj && flag == 0) {
                        QnSP_rap_local -= event_Q1_diff_real[ii][0];
                    }
                    temp_corr[rap_idx] += QnSP_rap_local;
                }
            }
        }
        for (int j = 0; j < N_rap; j++) {
            corr[i][j] += temp_corr[j];
            corr_err[i][j] += temp_corr[j] * temp_corr[j];
        }
    }
}

//! This function outputs the event averaged two-particle flow correlation
void singleParticleSpectra::output_two_particle_correlation() {
    const auto write = [&](const char *observable, const char *columns,
                           const vector<double> &sum,
                           const vector<double> &squares,
                           const vector<double> &cross) {
        ostringstream filename;
        filename << path_ << "/particle_" << particle_monval << observable
                 << (rap_type == 0 ? "_eta_" : "_y_") << rap_min << "_"
                 << rap_max
                 << (bWeakFeedDown_ ? "_weakFD.dat" : ".dat");
        ofstream out(filename.str());
        out << "# n  vn{2}  vn{2}_err  " << columns << '\n';
        const PairMoments pairs{sum[0], squares[0], cross[0]};
        for (int n = 0; n < order_max; n++) {
            out << scientific << setw(18) << setprecision(8) << n;
            write_pair_moments(out, {sum[n], squares[n], cross[n]}, pairs,
                               total_number_of_events, n == 0);
            out << '\n';
        }
    };
    write("_vn2", "<Qn*conj(Qn)>  <(Qn*conj(Qn))^2>", Qn2_vector,
          Qn2_vector_err, Qn2_vector_cross);
    if (flag_charge_dependence == 1) {
        write("_Cn2_ss", "Cn2_ss  Cn2_ss_err", Cn2_ss, Cn2_ss_err,
              Cn2_ss_cross);
        write("_Cn2_os", "Cn2_os  Cn2_os_err", Cn2_os, Cn2_os_err,
              Cn2_os_cross);
    }
}

//! Output two-particle correlations as a function of the rapidity difference.
void singleParticleSpectra::output_two_particle_correlation_rap() {
    const auto write = [&](const char *observable,
                           const vector<vector<double>> &sum,
                           const vector<vector<double>> &squares,
                           const vector<vector<double>> &cross) {
        ostringstream filename;
        filename << path_ << "/particle_" << particle_monval << observable
                 << "_pT_" << vn_rapidity_dis_pT_min << "_"
                 << vn_rapidity_dis_pT_max
                 << output_suffix(rapidity_correlation_output_);
        ofstream out(filename.str());
        out << (rap_type == 0 ? "# eta" : "# y")
            << "  vn{2}  vn{2}_err  <Qn*conj(Qn)>  <(Qn*conj(Qn))^2>\n";
        for (int bin = 0; bin < N_rap; bin++) {
            out << scientific << setw(18) << setprecision(8)
                << rapidity_dis_min + bin * drap;
            const PairMoments pairs {
                sum[0][bin], squares[0][bin], cross[0][bin]};
            for (int n = 0; n < order_max; n++) {
                write_pair_moments(
                    out, {sum[n][bin], squares[n][bin], cross[n][bin]}, pairs,
                    total_number_of_events, n == 0);
            }
            out << '\n';
        }
    };
    write("_vn2_eta12", QnSP_eta12, QnSP_eta12_err, QnSP_eta12_cross);
    if (flag_charge_dependence == 1) {
        write("_vn2_eta12_ss", Cn2_ss_eta12, Cn2_ss_eta12_err,
              Cn2_ss_eta12_cross);
        write("_vn2_eta12_os", Cn2_os_eta12, Cn2_os_eta12_err,
              Cn2_os_eta12_cross);
    }
}

//! This function computes the 3-particle correlation for Qn vectors
//! within one event with same-sign and opposite-sign pairs
void singleParticleSpectra::calculate_three_particle_correlation_charge_dep(
    vector<double> &event_Qn_p_real, vector<double> &event_Qn_p_imag,
    vector<double> &event_Qn_m_real, vector<double> &event_Qn_m_imag,
    vector<double> &event_Qn_real, vector<double> &event_Qn_imag,
    vector<double> &corr_ss, vector<double> &corr_ss_err,
    vector<double> &corr_os, vector<double> &corr_os_err,
    vector<double> &corr_ss_13, vector<double> &corr_ss_13_err,
    vector<double> &corr_os_13, vector<double> &corr_os_13_err,
    vector<double> *cross_ss, vector<double> *cross_os,
    vector<double> *cross_ss_13, vector<double> *cross_os_13) {
    vector<double> temp_Cmnk_ss_1(num_corr, 0.);
    vector<double> temp_Cmnk_ss_2(num_corr, 0.);
    vector<double> temp_Cmnk_ss_3(num_corr, 0.);
    vector<double> temp_Cmnk_ss_4(num_corr, 0.);
    vector<double> temp_Cmnk_ss_1_err(num_corr, 0.);
    vector<double> temp_Cmnk_ss_2_err(num_corr, 0.);
    vector<double> temp_Cmnk_ss_3_err(num_corr, 0.);
    vector<double> temp_Cmnk_ss_4_err(num_corr, 0.);
    calculate_three_particle_correlation(
        event_Qn_p_real, event_Qn_p_imag, event_Qn_p_real, event_Qn_p_imag,
        event_Qn_real, event_Qn_imag, 0, temp_Cmnk_ss_1, temp_Cmnk_ss_1_err);
    calculate_three_particle_correlation(
        event_Qn_m_real, event_Qn_m_imag, event_Qn_m_real, event_Qn_m_imag,
        event_Qn_real, event_Qn_imag, 0, temp_Cmnk_ss_2, temp_Cmnk_ss_2_err);
    calculate_three_particle_correlation(
        event_Qn_p_real, event_Qn_p_imag, event_Qn_real, event_Qn_imag,
        event_Qn_p_real, event_Qn_p_imag, 3, temp_Cmnk_ss_3,
        temp_Cmnk_ss_3_err);
    calculate_three_particle_correlation(
        event_Qn_m_real, event_Qn_m_imag, event_Qn_real, event_Qn_imag,
        event_Qn_m_real, event_Qn_m_imag, 3, temp_Cmnk_ss_4,
        temp_Cmnk_ss_4_err);
    vector<double> local_ss(num_corr), local_ss_13(num_corr);
    for (int i = 0; i < num_corr; i++) {
        local_ss[i] = temp_Cmnk_ss_1[i] + temp_Cmnk_ss_2[i];
        local_ss_13[i] = temp_Cmnk_ss_3[i] + temp_Cmnk_ss_4[i];
    }
    add_event_moments(local_ss, corr_ss, corr_ss_err, cross_ss);
    add_event_moments(local_ss_13, corr_ss_13, corr_ss_13_err, cross_ss_13);

    vector<double> temp_Cmnk_os_1(num_corr, 0.);
    vector<double> temp_Cmnk_os_2(num_corr, 0.);
    vector<double> temp_Cmnk_os_3(num_corr, 0.);
    vector<double> temp_Cmnk_os_4(num_corr, 0.);
    vector<double> temp_Cmnk_os_1_err(num_corr, 0.);
    vector<double> temp_Cmnk_os_2_err(num_corr, 0.);
    vector<double> temp_Cmnk_os_3_err(num_corr, 0.);
    vector<double> temp_Cmnk_os_4_err(num_corr, 0.);
    calculate_three_particle_correlation(
        event_Qn_p_real, event_Qn_p_imag, event_Qn_m_real, event_Qn_m_imag,
        event_Qn_real, event_Qn_imag, 1, temp_Cmnk_os_1, temp_Cmnk_os_1_err);
    calculate_three_particle_correlation(
        event_Qn_m_real, event_Qn_m_imag, event_Qn_p_real, event_Qn_p_imag,
        event_Qn_real, event_Qn_imag, 1, temp_Cmnk_os_2, temp_Cmnk_os_2_err);
    calculate_three_particle_correlation(
        event_Qn_p_real, event_Qn_p_imag, event_Qn_real, event_Qn_imag,
        event_Qn_m_real, event_Qn_m_imag, 4, temp_Cmnk_os_3,
        temp_Cmnk_os_3_err);
    calculate_three_particle_correlation(
        event_Qn_m_real, event_Qn_m_imag, event_Qn_real, event_Qn_imag,
        event_Qn_p_real, event_Qn_p_imag, 4, temp_Cmnk_os_4,
        temp_Cmnk_os_4_err);
    vector<double> local_os(num_corr), local_os_13(num_corr);
    for (int i = 0; i < num_corr; i++) {
        local_os[i] = temp_Cmnk_os_1[i] + temp_Cmnk_os_2[i];
        local_os_13[i] = temp_Cmnk_os_3[i] + temp_Cmnk_os_4[i];
    }
    add_event_moments(local_os, corr_os, corr_os_err, cross_os);
    add_event_moments(local_os_13, corr_os_13, corr_os_13_err, cross_os_13);
}

//! This function computes the 3-particle correlation for Qn vectors
//! within one event
//!     C_nmk = Real(Q_n*Q_m*conj(Q_k)) for (112), (123), (224), (235)
//! self correlation is subtracted assuming Qk's sample >= Qn's and Qm's
//! flag = 0: Qn = Qm <= Qk, flag = 1: Qn != Qm, Qn \in Qk, Qm \in Qk
//! flag = 2: no overlap
//! flag = 3: Qn = Qk <= Qm, flag = 4: Qn != Qk, Qn \in Qm, Qk \in Qn
void singleParticleSpectra::calculate_three_particle_correlation(
    vector<double> &event_Q1_real, vector<double> &event_Q1_imag,
    vector<double> &event_Q2_real, vector<double> &event_Q2_imag,
    vector<double> &event_Q3_real, vector<double> &event_Q3_imag, int flag,
    vector<double> &corr, vector<double> &corr_err) {
    // C_nmk[0] = C_000 = N(N-1)(N-2) is the number of pairs
    // C_nmk[1] = C_112, C_nmk[2] = C_123, C_nmk[3] = C_224, C_nmk[4] = C_235
    // C_nmk[5] = C_134, C_nmk[6] = C_246, C_nmk[7] = C_336, C_nmk[8] = C_347
    int n[9] = {0, 1, 1, 2, 2, 1, 2, 3, 3};
    int m[9] = {0, 1, 2, 2, 3, 3, 4, 3, 4};
    for (int i = 0; i < num_corr; i++) {
        int k = n[i] + m[i];
        if (k >= order_max) {
            messager << "the vn needed for three-particle correlation "
                     << "is not computed, order_max = " << order_max
                     << "n+m = " << k << ". Make sure order_max > n+m!";
            messager.flush("error");
            exit(1);
        }
        const double q1_imag_n = correlation_imag(event_Q1_imag, n[i]);
        const double q2_imag_m = correlation_imag(event_Q2_imag, m[i]);
        const double q3_imag_k = correlation_imag(event_Q3_imag, k);
        const double q2_imag_n = correlation_imag(event_Q2_imag, n[i]);
        const double q1_imag_m = correlation_imag(event_Q1_imag, m[i]);
        const double q1_imag_k = correlation_imag(event_Q1_imag, k);
        const double q3_imag_n = correlation_imag(event_Q3_imag, n[i]);
        double Qn_Qm_Qkstar =
            ((event_Q1_real[n[i]] * event_Q2_real[m[i]]
              - q1_imag_n * q2_imag_m)
                 * event_Q3_real[k]
             + (event_Q1_real[n[i]] * q2_imag_m
                + q1_imag_n * event_Q2_real[m[i]])
                   * q3_imag_k);

        // full overlap between Q2 and Q3 (Q2 <= Q3) := Q1[n]*conj(Q2[n])
        double Qn_Qnstar =
            (event_Q1_real[n[i]] * event_Q2_real[n[i]]
             + q1_imag_n * q2_imag_n);
        // full overlap between Q1 and Q3 (Q1 <= Q3) := Q2[m]*conj(Q1[m])
        double Qm_Qmstar =
            (event_Q2_real[m[i]] * event_Q1_real[m[i]]
             + q2_imag_m * q1_imag_m);
        // full overlap between Q1 and Q2 (Q1 <= Q2) := Q1[k]*conj(Q3[k])
        double Qk_Qkstar =
            (event_Q1_real[k] * event_Q3_real[k]
             + q1_imag_k * q3_imag_k);

        // full overlap between Q2 and Q3 (Q3 <= Q2) := Q1[n]*conj(Q3[n])
        double Qn_Qnstar_3 =
            (event_Q1_real[n[i]] * event_Q3_real[n[i]]
             + q1_imag_n * q3_imag_n);
        double corr_local = 0.;
        if (flag == 0) {
            corr_local =
                (Qn_Qm_Qkstar - Qn_Qnstar - Qm_Qmstar - Qk_Qkstar
                 + 2. * event_Q1_real[0]);
        } else if (flag == 1) {
            corr_local = (Qn_Qm_Qkstar - Qn_Qnstar - Qm_Qmstar);
        } else if (flag == 2) {
            corr_local = Qn_Qm_Qkstar;
        } else if (flag == 3) {
            corr_local =
                (Qn_Qm_Qkstar - Qn_Qnstar_3 - Qm_Qmstar - Qk_Qkstar
                 + 2. * event_Q1_real[0]);
        } else if (flag == 4) {
            corr_local = (Qn_Qm_Qkstar - Qn_Qnstar_3 - Qk_Qkstar);
        }

        corr[i] += corr_local;
        corr_err[i] += corr_local * corr_local;
    }
}

//! This function computes the charge dependent 3-particle correlations
void singleParticleSpectra::calculate_three_particle_correlation_deltaeta_chdep(
    vector<vector<double>> &event_Qn_p_rap_real,
    vector<vector<double>> &event_Qn_p_rap_imag,
    vector<vector<double>> &event_Qn_m_rap_real,
    vector<vector<double>> &event_Qn_m_rap_imag,
    vector<vector<double>> &event_Qn_rap_real,
    vector<vector<double>> &event_Qn_rap_imag,
    vector<vector<double>> &Cmnk_ss_eta12,
    vector<vector<double>> &Cmnk_ss_eta12_err,
    vector<vector<double>> &Cmnk_os_eta12,
    vector<vector<double>> &Cmnk_os_eta12_err,
    vector<vector<double>> &Cmnk_ss_eta13,
    vector<vector<double>> &Cmnk_ss_eta13_err,
    vector<vector<double>> &Cmnk_os_eta13,
    vector<vector<double>> &Cmnk_os_eta13_err,
    vector<vector<double>> *cross_ss_eta12,
    vector<vector<double>> *cross_os_eta12,
    vector<vector<double>> *cross_ss_eta13,
    vector<vector<double>> *cross_os_eta13) {
    vector<vector<double>> temp_Cmnk_ss_1(num_corr, vector<double>(N_rap, 0.));
    vector<vector<double>> temp_Cmnk_ss_2(num_corr, vector<double>(N_rap, 0.));
    vector<vector<double>> temp_Cmnk_ss_3(num_corr, vector<double>(N_rap, 0.));
    vector<vector<double>> temp_Cmnk_ss_4(num_corr, vector<double>(N_rap, 0.));
    vector<vector<double>> temp_Cmnk_os_1(num_corr, vector<double>(N_rap, 0.));
    vector<vector<double>> temp_Cmnk_os_2(num_corr, vector<double>(N_rap, 0.));
    vector<vector<double>> temp_Cmnk_os_3(num_corr, vector<double>(N_rap, 0.));
    vector<vector<double>> temp_Cmnk_os_4(num_corr, vector<double>(N_rap, 0.));
    vector<vector<double>> temp_Cmnk_err(num_corr, vector<double>(N_rap, 0.));

    calculate_three_particle_correlation_deltaeta(
        event_Qn_p_rap_real, event_Qn_p_rap_imag, event_Qn_m_rap_real,
        event_Qn_m_rap_imag, event_Qn_rap_real, event_Qn_rap_imag, 1, 2,
        temp_Cmnk_os_1, temp_Cmnk_err);
    calculate_three_particle_correlation_deltaeta(
        event_Qn_m_rap_real, event_Qn_m_rap_imag, event_Qn_p_rap_real,
        event_Qn_p_rap_imag, event_Qn_rap_real, event_Qn_rap_imag, 1, 2,
        temp_Cmnk_os_2, temp_Cmnk_err);
    calculate_three_particle_correlation_deltaeta(
        event_Qn_p_rap_real, event_Qn_p_rap_imag, event_Qn_m_rap_real,
        event_Qn_m_rap_imag, event_Qn_rap_real, event_Qn_rap_imag, 2, 2,
        temp_Cmnk_os_3, temp_Cmnk_err);
    calculate_three_particle_correlation_deltaeta(
        event_Qn_m_rap_real, event_Qn_m_rap_imag, event_Qn_p_rap_real,
        event_Qn_p_rap_imag, event_Qn_rap_real, event_Qn_rap_imag, 2, 2,
        temp_Cmnk_os_4, temp_Cmnk_err);
    calculate_three_particle_correlation_deltaeta(
        event_Qn_p_rap_real, event_Qn_p_rap_imag, event_Qn_p_rap_real,
        event_Qn_p_rap_imag, event_Qn_rap_real, event_Qn_rap_imag, 1, 1,
        temp_Cmnk_ss_1, temp_Cmnk_err);
    calculate_three_particle_correlation_deltaeta(
        event_Qn_m_rap_real, event_Qn_m_rap_imag, event_Qn_m_rap_real,
        event_Qn_m_rap_imag, event_Qn_rap_real, event_Qn_rap_imag, 1, 1,
        temp_Cmnk_ss_2, temp_Cmnk_err);
    calculate_three_particle_correlation_deltaeta(
        event_Qn_p_rap_real, event_Qn_p_rap_imag, event_Qn_p_rap_real,
        event_Qn_p_rap_imag, event_Qn_rap_real, event_Qn_rap_imag, 2, 1,
        temp_Cmnk_ss_3, temp_Cmnk_err);
    calculate_three_particle_correlation_deltaeta(
        event_Qn_m_rap_real, event_Qn_m_rap_imag, event_Qn_m_rap_real,
        event_Qn_m_rap_imag, event_Qn_rap_real, event_Qn_rap_imag, 2, 1,
        temp_Cmnk_ss_4, temp_Cmnk_err);

    const auto add_pair = [](vector<vector<double>> &first,
                             const vector<vector<double>> &second) {
        for (std::size_t i = 0; i < first.size(); i++) {
            for (std::size_t j = 0; j < first[i].size(); j++) {
                first[i][j] += second[i][j];
            }
        }
    };
    add_pair(temp_Cmnk_ss_1, temp_Cmnk_ss_2);
    add_pair(temp_Cmnk_ss_3, temp_Cmnk_ss_4);
    add_pair(temp_Cmnk_os_1, temp_Cmnk_os_2);
    add_pair(temp_Cmnk_os_3, temp_Cmnk_os_4);
    add_event_moments(
        temp_Cmnk_ss_1, Cmnk_ss_eta12, Cmnk_ss_eta12_err, cross_ss_eta12);
    add_event_moments(
        temp_Cmnk_ss_3, Cmnk_ss_eta13, Cmnk_ss_eta13_err, cross_ss_eta13);
    add_event_moments(
        temp_Cmnk_os_1, Cmnk_os_eta12, Cmnk_os_eta12_err, cross_os_eta12);
    add_event_moments(
        temp_Cmnk_os_3, Cmnk_os_eta13, Cmnk_os_eta13_err, cross_os_eta13);
}

//! This function computes the 3-particle correlation as a function of relative
//! rapidity between particle 1 and particle 2 within one event
//!     C_nmk(eta_12) = Real(Q_n(eta1)*Q_m(eta2)*conj(Q_k))
//! self correlation is subtracted assuming Qk's sample >= Qn's and Qm's
//! Qn != Qm, Qn \in Qk, Qm \in Qk
//! flag == 1 : eta_12, flag == 2: eta_13
//! flag_ch == 1: same charge, flag_ch == 2: opposite charge
void singleParticleSpectra::calculate_three_particle_correlation_deltaeta(
    vector<vector<double>> &event_Q1_real,
    vector<vector<double>> &event_Q1_imag,
    vector<vector<double>> &event_Q2_real,
    vector<vector<double>> &event_Q2_imag,
    vector<vector<double>> &event_Q3_real,
    vector<vector<double>> &event_Q3_imag, int flag, int flag_ch,
    vector<vector<double>> &corr_rap, vector<vector<double>> &corr_rap_err) {
    // C_nmk[0] = C_000 = N(N-1)(N-2) is the number of ordered triplets
    // C_nmk[1] = C_112, C_nmk[2] = C_123, C_nmk[3] = C_224, C_nmk[4] = C_235
    // C_nmk[5] = C_134, C_nmk[6] = C_246, C_nmk[7] = C_336, C_nmk[8] = C_347
    int n[9] = {0, 1, 1, 2, 2, 1, 2, 3, 3};
    int m[9] = {0, 1, 2, 2, 3, 3, 4, 3, 4};
    using Complex = std::complex<double>;
    const auto q = [](const vector<vector<double>> &real,
                      const vector<vector<double>> &imag, int bin,
                      int harmonic) {
        return Complex(
            real[bin][harmonic], correlation_imag(imag[bin], harmonic));
    };
    const bool eta12 = flag == 1;
    const bool same_charge = flag_ch == 1;
    vector<double> temp_corr(N_rap);
    for (int i = 0; i < num_corr; i++) {
        const int k = n[i] + m[i];
        if (k >= order_max) {
            messager << "Three-particle correlation requires order_max > "
                     << k << ", got " << order_max;
            messager.flush("error");
            exit(1);
        }

        Complex integrated_q;
        for (int bin = 0; bin < N_rap; bin++) {
            integrated_q += eta12
                                ? q(event_Q3_real, event_Q3_imag, bin, k)
                                : q(event_Q2_real, event_Q2_imag, bin, m[i]);
        }

        std::fill(temp_corr.begin(), temp_corr.end(), 0.0);
        for (int ii = 0; ii < N_rap; ii++) {
            const Complex q1n = q(event_Q1_real, event_Q1_imag, ii, n[i]);
            const Complex q1m = q(event_Q1_real, event_Q1_imag, ii, m[i]);
            const Complex q1k = q(event_Q1_real, event_Q1_imag, ii, k);
            for (int jj = 0; jj < N_rap; jj++) {
                const int rap_idx =
                    (ii - jj) - static_cast<int>(rapidity_dis_min / drap);
                if (rap_idx < 0 || rap_idx >= N_rap) continue;

                const Complex q2n = q(event_Q2_real, event_Q2_imag, jj, n[i]);
                const Complex q2m =
                    eta12 ? q(event_Q2_real, event_Q2_imag, jj, m[i])
                          : integrated_q;
                const Complex q3k =
                    eta12 ? integrated_q
                          : q(event_Q3_real, event_Q3_imag, jj, k);
                // Subtract the 2=3 overlap, then the 1=2 and 1=3 overlaps
                // wherever their charge selections and rapidity bins coincide.
                Complex corr =
                    q1n * q2m * std::conj(q3k) - q1n * std::conj(q2n);
                if (same_charge && (!eta12 || ii == jj)) {
                    corr -= q1k * std::conj(q3k);
                }
                if (eta12 || ii == jj) {
                    corr -= std::conj(q1m) * q2m;
                }
                if (same_charge && ii == jj) {
                    corr += 2.0 * event_Q1_real[ii][0];
                }
                temp_corr[rap_idx] += corr.real();
            }
        }
        for (int bin = 0; bin < N_rap; bin++) {
            corr_rap[i][bin] += temp_corr[bin];
            corr_rap_err[i][bin] += temp_corr[bin] * temp_corr[bin];
        }
    }
}

//! This function computes the 4-particle correlation for C_n{4}
//! using Qn vectors within one event
//!     C_n{4} = <Q1_n*conj(Q2_n)*Q3_n*conj(Q4_n)>
//! for n = 0, 1, 2, 3, 4
//! self correlation is subtracted
void singleParticleSpectra::calculate_four_particle_correlation_Cn4(
    vector<double> &event_Qn_real, vector<double> &event_Qn_imag,
    vector<double> &corr, vector<double> &corr_err) {
    // C_n[0] = N(N-1)(N-2)(N-3) counts ordered quadruplets
    double dN = event_Qn_real[0];
    for (int i = 0; i < num_Cn4; i++) {
        int ii = 2 * i;
        if (ii >= order_max) {
            messager << "Error: the vn needed for C_n{4} is not "
                     << "computed, order_max = " << order_max << "2n = " << ii
                     << ". Make sure order_max > 2n!";
            messager.flush("error");
            exit(1);
        }
        const double imag_n = correlation_imag(event_Qn_imag, i);
        const double imag_2n = correlation_imag(event_Qn_imag, ii);
        double abs_Qn_2 =
            event_Qn_real[i] * event_Qn_real[i] + imag_n * imag_n;
        double abs_Q2n_2 =
            event_Qn_real[ii] * event_Qn_real[ii] + imag_2n * imag_2n;

        double Re_Q2n_QnQn =
            (event_Qn_real[i]
                 * (event_Qn_real[ii] * event_Qn_real[i]
                    + imag_2n * imag_n)
             + imag_n
                   * (event_Qn_real[i] * imag_2n
                      - event_Qn_real[ii] * imag_n));

        double corr_local =
            (abs_Qn_2 * abs_Qn_2 - 2. * Re_Q2n_QnQn - 4. * (dN - 2.) * abs_Qn_2
             + abs_Q2n_2 + 2. * dN * (dN - 3.));
        corr[i] += corr_local;
        corr_err[i] += corr_local * corr_local;
    }
}

//! This function computes the 4-particle correlation for symmetric cumulants
//! using Qn vectors within one event
//!     SC_mn = <Q1_m*conj(Q2_m)*Q3_n*conj(Q4_n)>
//! for (32), (42), (52), (43), (53)
//! self correlation is subtracted assuming Qk's sample >= Qn's and Qm's
void singleParticleSpectra::calculate_four_particle_correlation_SC(
    vector<double> &event_Q1_real, vector<double> &event_Q1_imag,
    vector<double> &event_Q2_real, vector<double> &event_Q2_imag,
    vector<double> &event_Q3_real, vector<double> &event_Q3_imag,
    vector<double> &event_Q4_real, vector<double> &event_Q4_imag, int flag,
    vector<double> &corr, vector<double> &corr_err) {
    // SC_mn[0] = SC_00 = N(N-1)(N-2)(N-3) counts ordered quadruplets
    // SC_mn[1] = SC_32, SC_mn[2] = SC_42
    // SC_mn[3] = SC_52, SC_mn[4] = SC_43, SC_mn[5] = SC_53
    int m[6] = {0, 3, 4, 5, 4, 5};
    int n[6] = {0, 2, 2, 2, 3, 3};
    for (int i = 0; i < SC_num_corr; i++) {
        int k = m[i] + n[i];
        if (k >= order_max) {
            messager << "Error: the vn needed for symmetric cumulants is not "
                     << "computed, order_max = " << order_max << "m+n = " << k
                     << ". Make sure order_max > m+n!";
            messager.flush("error");
            exit(1);
        }
        int l = m[i] - n[i];
        const double q1_imag_m = correlation_imag(event_Q1_imag, m[i]);
        const double q2_imag_m = correlation_imag(event_Q2_imag, m[i]);
        const double q3_imag_n = correlation_imag(event_Q3_imag, n[i]);
        const double q4_imag_n = correlation_imag(event_Q4_imag, n[i]);
        const double q1_imag_k = correlation_imag(event_Q1_imag, k);
        const double q1_imag_l = correlation_imag(event_Q1_imag, l);
        const double q2_imag_k = correlation_imag(event_Q2_imag, k);
        const double q2_imag_l = correlation_imag(event_Q2_imag, l);
        const double q1_imag_n = correlation_imag(event_Q1_imag, n[i]);
        double Q1_conjQ2 =
            (event_Q1_real[m[i]] * event_Q2_real[m[i]]
             + q1_imag_m * q2_imag_m);
        double Q3_conjQ4 =
            (event_Q3_real[n[i]] * event_Q4_real[n[i]]
             + q3_imag_n * q4_imag_n);
        // full overlap of Q1 and Q3
        double Q1Q3_overlap =
            (event_Q1_real[k] * event_Q2_real[m[i]] * event_Q4_real[n[i]]
             - event_Q1_real[k] * q2_imag_m * q4_imag_n
             + q1_imag_k * event_Q2_real[m[i]] * q4_imag_n
             + q1_imag_k * q2_imag_m * event_Q4_real[n[i]]);
        // full overlap of Q1 and Q4
        double Q1Q4_overlap =
            (event_Q1_real[l] * event_Q2_real[m[i]] * event_Q3_real[n[i]]
             + event_Q1_real[l] * q2_imag_m * q3_imag_n
             - q1_imag_l * event_Q2_real[m[i]] * q3_imag_n
             + q1_imag_l * q2_imag_m * event_Q3_real[n[i]]);
        // full overlap of Q2 and Q4
        double Q2Q4_overlap =
            (event_Q2_real[k] * event_Q1_real[m[i]] * event_Q3_real[n[i]]
             - event_Q2_real[k] * q1_imag_m * q3_imag_n
             + q2_imag_k * event_Q1_real[m[i]] * q3_imag_n
             + q2_imag_k * q1_imag_m * event_Q3_real[n[i]]);
        // full overlap of Q2 and Q3
        double Q2Q3_overlap =
            (event_Q2_real[l] * event_Q1_real[m[i]] * event_Q4_real[n[i]]
             + event_Q2_real[l] * q1_imag_m * q4_imag_n
             - q2_imag_l * event_Q1_real[m[i]] * q4_imag_n
             + q2_imag_l * q1_imag_m * event_Q4_real[n[i]]);
        // full overlap of Q1 with Q3 and Q2 with Q4
        double Q1Q3_Q2Q4_overlap =
            (event_Q1_real[k] * event_Q2_real[k]
             + q1_imag_k * q2_imag_k);
        // full overlap of Q1 with Q4 and Q2 with Q3
        double Q1Q4_Q2Q3_overlap =
            (event_Q1_real[l] * event_Q2_real[l]
             + q1_imag_l * q2_imag_l);
        // full overlap of Q1 with Q2 and Q3
        double Q1Q2Q3_overlap =
            (event_Q1_real[n[i]] * event_Q4_real[n[i]]
             + q1_imag_n * q4_imag_n);
        // full overlap of Q1 with Q2 and Q4
        double Q1Q2Q4_overlap =
            (event_Q1_real[n[i]] * event_Q3_real[n[i]]
             + q1_imag_n * q3_imag_n);
        // full overlap of Q1 with Q3 and Q4
        double Q1Q3Q4_overlap =
            (event_Q1_real[m[i]] * event_Q2_real[m[i]]
             + q1_imag_m * q2_imag_m);
        // full overlap of Q2 with Q3 and Q4
        double Q2Q3Q4_overlap =
            (event_Q1_real[m[i]] * event_Q2_real[m[i]]
             + q1_imag_m * q2_imag_m);
        // full overlap of Q1 with Q2, Q3, and Q4
        double Q1Q2Q3Q4_overlap = event_Q1_real[0];

        double corr_local = 0.;
        if (flag == 0) {
            corr_local =
                (Q1_conjQ2 * Q3_conjQ4 - Q1Q3_overlap - Q1Q4_overlap
                 - Q2Q3_overlap - Q2Q4_overlap + Q1Q3_Q2Q4_overlap
                 + Q1Q4_Q2Q3_overlap
                 - (event_Q1_real[0] - 4.) / 2. * Q1Q2Q3_overlap
                 - (event_Q1_real[0] - 4.) / 2. * Q1Q2Q4_overlap
                 - (event_Q1_real[0] - 4.) / 2. * Q1Q3Q4_overlap
                 - (event_Q1_real[0] - 4.) / 2. * Q2Q3Q4_overlap
                 + Q1Q2Q3Q4_overlap * (Q1Q2Q3Q4_overlap - 6.));
        }
        corr[i] += corr_local;
        corr_err[i] += corr_local * corr_local;
    }
}

//! This function outputs the event averaged three-particle correlation
void singleParticleSpectra::output_three_particle_correlation() {
    const auto filename = [&](const char *observable) {
        ostringstream name;
        name << path_ << "/particle_" << particle_monval << observable
             << (rap_type == 0 ? "_eta_" : "_y_") << rap_min << "_" << rap_max
             << (bWeakFeedDown_ ? "_weakFD.dat" : ".dat");
        return name.str();
    };
    // One column pair: tuple count (row 0) or per-tuple correlator (rows 1+).
    const auto write_entry = [&](ofstream &out, const vector<double> &sum,
                                 const vector<double> &squares,
                                 const vector<double> &cross, int i) {
        double value, error;
        if (i == 0) {
            tuple_count(
                sum[0], squares[0], total_number_of_events, value, error);
        } else {
            tuple_ratio(
                sum[i], squares[i], cross[i], sum[0], squares[0], value, error);
        }
        out << value << "  " << error;
    };

    ofstream output(filename("_Cmnk").c_str());
    output << "# n  C_nmk  C_nmk_err" << '\n';
    for (int i = 0; i < num_corr; i++) {
        output << scientific << setw(18) << setprecision(8) << i << "  ";
        write_entry(output, C_nmk, C_nmk_err, C_nmk_cross, i);
        output << '\n';
    }
    if (flag_charge_dependence != 1) return;

    const auto write_charge =
        [&](const char *observable, const char *label,
            const vector<double> &sum_12, const vector<double> &squares_12,
            const vector<double> &cross_12, const vector<double> &sum_13,
            const vector<double> &squares_13, const vector<double> &cross_13) {
            ofstream out(filename(observable).c_str());
            out << "# n  C_nmk_" << label << "  C_nmk_" << label
                << "_err  C_nmk_" << label << "_13  C_nmk_" << label
                << "_13_err" << '\n';
            for (int i = 0; i < num_corr; i++) {
                out << scientific << setw(18) << setprecision(8) << i << "  ";
                write_entry(out, sum_12, squares_12, cross_12, i);
                out << "  ";
                write_entry(out, sum_13, squares_13, cross_13, i);
                out << '\n';
            }
        };
    write_charge(
        "_Cmnk_ss", "ss", C_nmk_ss, C_nmk_ss_err, C_nmk_ss_cross, C_nmk_ss_13,
        C_nmk_ss_13_err, C_nmk_ss_13_cross);
    write_charge(
        "_Cmnk_os", "os", C_nmk_os, C_nmk_os_err, C_nmk_os_cross, C_nmk_os_13,
        C_nmk_os_13_err, C_nmk_os_13_cross);
}

//! This function outputs the rapidity dependent three-particle correlation
void singleParticleSpectra::output_three_particle_correlation_rap() {
    using Matrix = vector<vector<double>>;
    const auto write = [&](const char *observable, const char *header,
                           const Matrix &sum, const Matrix &squares,
                           const Matrix &cross) {
        ostringstream filename;
        filename << path_ << "/particle_" << particle_monval << observable
                 << "_pT_" << vn_rapidity_dis_pT_min << "_"
                 << vn_rapidity_dis_pT_max
                 << output_suffix(rapidity_correlation_output_);
        ofstream out(filename.str().c_str());
        out << header << '\n';
        for (int j = 0; j < N_rap; j++) {
            double count, count_error;
            tuple_count(
                sum[0][j], squares[0][j], total_number_of_events, count,
                count_error);
            out << scientific << setw(18) << setprecision(8)
                << rapidity_dis_min + j * drap << "  " << count << "  "
                << count_error << "  ";
            for (int i = 1; i < num_corr; i++) {
                double value, error;
                tuple_ratio(
                    sum[i][j], squares[i][j], cross[i][j], sum[0][j],
                    squares[0][j], value, error);
                out << scientific << setw(18) << setprecision(8) << value
                    << "  " << error << "  ";
            }
            out << '\n';
        }
    };
    const std::string gap = rap_type == 0 ? "eta" : "y";
    const std::string list = "(000, 112, 123, 224, 235, 134, 246, 336, 347)";
    write(
        "_Cmnk_eta12", ("# " + gap + "12  C_nmk  C_nmk_err " + list).c_str(),
        C_nmk_eta12, C_nmk_eta12_err, C_nmk_eta12_cross);
    write(
        "_Cmnk_eta13", ("# " + gap + "13  C_nmk  C_nmk_err " + list).c_str(),
        C_nmk_eta13, C_nmk_eta13_err, C_nmk_eta13_cross);
    if (flag_charge_dependence != 1) return;
    write(
        "_Cmnk_eta12_ss", ("# " + gap + "12  C_nmk_ss  C_nmk_ss_err").c_str(),
        C_nmk_eta12_ss, C_nmk_eta12_ss_err, C_nmk_eta12_ss_cross);
    write(
        "_Cmnk_eta12_os", ("# " + gap + "12  C_nmk_os  C_nmk_os_err").c_str(),
        C_nmk_eta12_os, C_nmk_eta12_os_err, C_nmk_eta12_os_cross);
    write(
        "_Cmnk_eta13_ss", ("# " + gap + "13  C_nmk_ss  C_nmk_ss_err").c_str(),
        C_nmk_eta13_ss, C_nmk_eta13_ss_err, C_nmk_eta13_ss_cross);
    write(
        "_Cmnk_eta13_os", ("# " + gap + "13  C_nmk_os  C_nmk_os_err").c_str(),
        C_nmk_eta13_os, C_nmk_eta13_os_err, C_nmk_eta13_os_cross);
}

//! This function outputs the event averaged four-particle correlation Cn{4}
void singleParticleSpectra::output_four_particle_Cn4_correlation() {
    write_four_particle_correlation(
        "_Cn4", "# n  Cn4  Cn4_err", Cn4, Cn4_err, Cn4_cross);
}

//! This function outputs the event averaged four-particle correlation
void singleParticleSpectra::output_four_particle_SC_correlation() {
    write_four_particle_correlation(
        "_SCmn", "# n  SC_mn  SC_mn_err", SC_mn, SC_mn_err, SC_mn_cross);
}

void singleParticleSpectra::write_four_particle_correlation(
    const char *observable, const char *header, const vector<double> &sum,
    const vector<double> &squares, const vector<double> &cross) {
    ostringstream filename;
    filename << path_ << "/particle_" << particle_monval << observable
             << (rap_type == 0 ? "_eta_" : "_y_") << rap_min << "_" << rap_max
             << (bWeakFeedDown_ ? "_weakFD.dat" : ".dat");
    ofstream output(filename.str().c_str());
    output << header << '\n';
    for (std::size_t i = 0; i < sum.size(); i++) {
        double value, error;
        if (i == 0) {
            tuple_count(
                sum[0], squares[0], total_number_of_events, value, error);
        } else {
            tuple_ratio(
                sum[i], squares[i], cross[i], sum[0], squares[0], value, error);
        }
        output << scientific << setw(18) << setprecision(8) << i << "  "
               << value << "  " << error << '\n';
    }
}

//! this function computes the pT and y differential Qn vector in one event
//! The array saves <y>, <pT>, <N>, <ET>, Qn_real, Qn_imag
//!                 (for harmonic order n = 1 to order_max - 1)
int singleParticleSpectra::rapidity_bin(double rapidity) const {
    const double coordinate = (rapidity - rapidity_dis_min) / drap + 0.5;
    if (!(coordinate >= 0.0 && coordinate < N_rap)) return -1;
    return static_cast<int>(coordinate);
}

void singleParticleSpectra::calculateRapidityPTDistribution(
    const int event_id) {
    int number_of_particles = particle_list->get_number_of_particles(event_id);
    for (int i = 0; i < number_of_particles; i++) {
        const auto track = particle_list->get_particle(event_id, i);
        const double rap_local = rap_type == 0 ? track.rap_eta : track.rap_y;
        const int rap_idx = rapidity_bin(rap_local);
        const double p_perp = track.pT;
        if (rap_idx < 0 || !(p_perp >= pT_min && p_perp < pT_max)) continue;

        const int pT_idx = static_cast<int>((p_perp - pT_min) / dpT);
        if (pT_idx >= npT) continue;
        const double E_trans = transverse_energy(track);

        Qn_pTrapdiff_vector[rap_idx][pT_idx][0] += rap_local;  // <eta>
        Qn_pTrapdiff_vector[rap_idx][pT_idx][1] += p_perp;     // <p_T>
        Qn_pTrapdiff_vector[rap_idx][pT_idx][2] += 1;          // N
        Qn_pTrapdiff_vector[rap_idx][pT_idx][3] += E_trans;    // E_T

        // Calculate vn.
        double p_phi = track.phi_p;
        for (int iorder = 1; iorder < order_max; iorder++) {
            double cos_nphi = cos(iorder * p_phi);
            double sin_nphi = sin(iorder * p_phi);
            Qn_pTrapdiff_vector[rap_idx][pT_idx][2 * iorder + 2] += cos_nphi;
            Qn_pTrapdiff_vector[rap_idx][pT_idx][2 * iorder + 3] += sin_nphi;
        }
    }
}

//! this function computes the pT-integrated Qn vector as a function of
//! rapidity in one event
void singleParticleSpectra::calculate_rapidity_distribution(
    int event_id, vector<vector<double>> &event_Qn_real,
    vector<vector<double>> &event_Qn_imag, int flag) {
    // first clean the results arrays
    for (int i = 0; i < N_rap; i++) {
        for (int j = 0; j < order_max; j++) {
            event_Qn_real[i][j] = 0.0;
            event_Qn_imag[i][j] = 0.0;
        }
    }

    if (flag < 0 || flag > 2) {
        throw std::invalid_argument("Invalid charge for rapidity distribution");
    }
    const int number_of_particles = selected_track_count(event_id, flag);
    vector<double> dNdy_event(flag == 0 ? N_rap : 0, 0.0);
    for (int i = 0; i < number_of_particles; i++) {
        const auto track = selected_track(event_id, i, flag);
        const double rap = rap_type == 0 ? track.rap_eta : track.rap_y;
        const int bin = rapidity_bin(rap);
        if (bin < 0) continue;

        // The particle yield has no pT cut.
        if (flag == 0) {
            rapidity_array[bin] += rap;
            dNdy_event[bin]++;
            if (track.reconstructed) reconstructed_count_rap_[bin] += 1.0;
            dETdy_array[bin] += transverse_energy(track);
        }
        if (track.pT > vn_rapidity_dis_pT_min
            && track.pT < vn_rapidity_dis_pT_max) {
            event_Qn_real[bin][0] += 1.0;
            event_Qn_imag[bin][0] += track.pT;
            for (int n = 1; n < order_max; n++) {
                event_Qn_real[bin][n] += cos(n * track.phi_p);
                event_Qn_imag[bin][n] += sin(n * track.phi_p);
            }
        }
    }
    if (flag != 0) return;
    for (int i = 0; i < N_rap; i++) {
        dNdy_array[i] += dNdy_event[i];
        if (!covariance_output_.enabled || dNdy_event[i] == 0.0) continue;
        for (int j = 0; j < N_rap; j++) {
            int idx = i * N_rap + j;
            dNdy_covMatrix[idx] += dNdy_event[i] * dNdy_event[j];
        }
    }
}

void singleParticleSpectra::outputRapidityPTDistribution() {
    std::string rapName = "eta";
    if (rap_type == 1) {
        rapName = "y";
    }
    ostringstream filename;
    filename << path_ << "/particle_" << particle_monval << "_pT" << rapName
             << "_distribution";
    filename << output_suffix(rapidity_pt_output_);
    ofstream output(filename.str().c_str());
    output << "# y  pT  N  E_T"
           << "  Qn_real  Qn_imag (n = 1 -" << order_max - 1 << ")" << '\n';
    output << "# Nev = " << total_number_of_events << " d" << rapName << " = "
           << drap << " dpT = " << dpT << " GeV" << '\n';

    for (int i = 0; i < N_rap; i++) {
        double rap_bin_mid = rapidity_array[i] / (dNdy_array[i] + 1.);
        for (int j = 0; j < npT; j++) {
            double pT_bin_mid = pT_array[j];
            if (Qn_pTrapdiff_vector[i][j][2] > 0) {
                rap_bin_mid =
                    (Qn_pTrapdiff_vector[i][j][0]
                     / Qn_pTrapdiff_vector[i][j][2]);
                pT_bin_mid =
                    (Qn_pTrapdiff_vector[i][j][1]
                     / Qn_pTrapdiff_vector[i][j][2]);
            }
            double dN = Qn_pTrapdiff_vector[i][j][2] / total_number_of_events;
            double dET = Qn_pTrapdiff_vector[i][j][3] / total_number_of_events;
            output << scientific << setprecision(8) << rap_bin_mid << "  "
                   << pT_bin_mid << "  " << dN << "  " << dET;
            // vn
            for (int iorder = 1; iorder < order_max; iorder++) {
                double Qn_real =
                    (Qn_pTrapdiff_vector[i][j][2 * iorder + 2]
                     / total_number_of_events);
                double Qn_imag =
                    (Qn_pTrapdiff_vector[i][j][2 * iorder + 3]
                     / total_number_of_events);
                output << "  " << Qn_real << "  " << Qn_imag;
            }
            output << '\n';
        }
    }
    output.close();
}

void singleParticleSpectra::output_rapidity_distribution() {
    if (rapidity_output_.enabled) output_rapidity_density();
    if (covariance_output_.enabled) output_rapidity_covariance();
}

void singleParticleSpectra::output_rapidity_density() {
    ostringstream filename;
    filename << path_ << "/particle_" << particle_monval;
    if (rap_type == 0) {
        filename << "_dNdeta";
    } else {
        filename << "_dNdy";
    }
    filename << "_pT_" << vn_rapidity_dis_pT_min << "_"
             << vn_rapidity_dis_pT_max;
    filename << output_suffix(rapidity_output_);
    ofstream output(filename.str().c_str());
    output << "# y  dN/dy  meanpT  vn_real  vn_imag  dET/dy"
           << "total event = " << total_number_of_events << '\n';
    for (int i = 0; i < N_rap; i++) {
        const double mean_rapidity = rapidity_array[i] / (dNdy_array[i] + 1.);
        const double weighted_count =
            dNdy_array[i]
            + (reconstructed_weight_ - 1.0) * reconstructed_count_rap_[i];
        const double yield_density =
            weighted_count / total_number_of_events / drap;
        const double energy_density =
            dETdy_array[i] / total_number_of_events / drap;

        double total_Nrap = std::max(1e-30, vn_real_rapidity_dis_array[i][0]);
        double meanpT = vn_imag_rapidity_dis_array[i][0] / total_Nrap;
        output << scientific << setw(18) << setprecision(8) << mean_rapidity
               << "   " << yield_density << "   " << meanpT << "   ";
        for (int iorder = 1; iorder < order_max; iorder++) {
            double vn_evavg_real =
                (vn_real_rapidity_dis_array[i][iorder] / total_Nrap);
            double vn_evavg_imag =
                (vn_imag_rapidity_dis_array[i][iorder] / total_Nrap);
            output << scientific << setw(18) << setprecision(8) << vn_evavg_real
                   << "   " << vn_evavg_imag << "   ";
        }
        // output dET/dy and total number of particles in the last column
        // this quantities is useful when one wants to reconst the Qn vectors
        output << scientific << setw(18) << setprecision(8)
               << energy_density << "   " << total_Nrap;
        output << '\n';
    }
    output.close();
}

void singleParticleSpectra::output_rapidity_covariance() {
    vector<double> mean_count(N_rap);
    vector<double> mean_rapidity(N_rap);
    for (int i = 0; i < N_rap; i++) {
        mean_count[i] = dNdy_array[i] / total_number_of_events;
        mean_rapidity[i] = rapidity_array[i] / (dNdy_array[i] + 1.);
    }
    // Population covariance of reconstructed per-event bin counts.
    ostringstream filename_cov;
    filename_cov << path_ << "/particle_" << particle_monval;
    if (rap_type == 0) {
        filename_cov << "_dNdeta_covMatrix";
    } else {
        filename_cov << "_dNdy_covMatrix";
    }
    filename_cov << output_suffix(covariance_output_);
    ofstream output_cov(filename_cov.str().c_str());
    output_cov << "# y1  y2  cov(N_bin1, N_bin2)  total event = "
               << total_number_of_events << '\n';
    for (int i = 0; i < N_rap; i++) {
        for (int j = 0; j < N_rap; j++) {
            int idx = i * N_rap + j;
            const double covariance =
                dNdy_covMatrix[idx] / total_number_of_events
                - mean_count[i] * mean_count[j];
            output_cov << scientific << setw(18) << setprecision(8)
                       << mean_rapidity[i] << "  " << mean_rapidity[j] << "  "
                       << covariance << '\n';
        }
    }
    output_cov.close();
}

void singleParticleSpectra::check_dNdSV(int event_id) {
    // Bin of a smeared value on a grid of n points spaced by width. A smear
    // wider than one bin can push the index past the last bin; such entries
    // are dropped, as for the r histogram.
    const auto in_grid = [](int idx, int n) { return idx >= 0 && idx < n; };
    int number_of_particles = particle_list->get_number_of_particles(event_id);
    for (int i = 0; i < number_of_particles; i++) {
        // double pz_local = particle_list->get_particle(event_id, i).pz;
        // double E_local = particle_list->get_particle(event_id, i).E;

        double rap_local;
        if (rap_type == 0) {
            // double mass = particle_list->get_particle(event_id, i).mass;
            // double pmag = sqrt(E_local * E_local - mass * mass);
            // rap_local = 0.5 * log((pmag + pz_local) / (pmag - pz_local));
            rap_local = particle_list->get_particle(event_id, i).rap_eta;
        } else {
            // rap_local = 0.5 * log((E_local + pz_local) / (E_local -
            // pz_local));
            rap_local = particle_list->get_particle(event_id, i).rap_y;
        }

        if (rap_local > rap_min && rap_local < rap_max) {
            spatial_count_ += 1.0;
            if (particle_list->get_particle(event_id, i).reconstructed) {
                spatial_reconstructed_count_ += 1.0;
            }
            // first dN/dtau
            double t_local = particle_list->get_particle(event_id, i).t;
            double z_local = particle_list->get_particle(event_id, i).z;
            double tau_local = sqrt(t_local * t_local - z_local * z_local);
            if (tau_local > tau_min && tau_local < tau_max) {
                double random = ran_gen_ptr->rand_uniform() * intrinsic_dtau;
                int idx = (int)((tau_local + random - tau_min) / dtau);
                if (in_grid(idx, N_tau)) {
                    tau_array[idx] += tau_local;
                    dNdtau_array[idx]++;
                }
            }

            // second dN/dx
            double x_local = particle_list->get_particle(event_id, i).x;
            double y_local = particle_list->get_particle(event_id, i).y;
            double r_local = sqrt(x_local * x_local + y_local * y_local);
            if (fabs(y_local) < 0.5) {
                if (x_local > spatial_x_min && x_local < spatial_x_max) {
                    double random = ran_gen_ptr->rand_uniform() * intrinsic_dx;
                    int idx =
                        (int)((x_local + random - spatial_x_min) / dspatial_x);
                    if (in_grid(idx, N_xpt)) {
                        xpt_array[idx] += x_local;
                        dNdx1_array[idx]++;
                    }
                }
            }
            if (fabs(x_local) < 0.5) {
                if (y_local > spatial_x_min && y_local < spatial_x_max) {
                    double random = ran_gen_ptr->rand_uniform() * intrinsic_dx;
                    int idx =
                        (int)((y_local + random - spatial_x_min) / dspatial_x);
                    if (in_grid(idx, N_xpt)) {
                        ypt_array[idx] += y_local;
                        dNdx2_array[idx]++;
                    }
                }
            }
            if (r_local > spatial_r_min && r_local < spatial_r_max) {
                double random = ran_gen_ptr->rand_uniform() * intrinsic_dx;
                int idx = static_cast<int>(
                    (r_local + random - spatial_r_min) / dspatial_r);
                if (in_grid(idx, N_xpt)) {
                    rpt_array[idx] += r_local;
                    dNdr_array[idx]++;
                }
            }

            // dN/(dtau dx)
            if (fabs(y_local) < 0.5) {
                if (tau_local > tau_min && tau_local < tau_max) {
                    double random =
                        ran_gen_ptr->rand_uniform() * intrinsic_dtau;
                    int idx_tau = (int)((tau_local + random - tau_min) / dtau);
                    if (x_local > spatial_x_min && x_local < spatial_x_max) {
                        double random =
                            ran_gen_ptr->rand_uniform() * intrinsic_dx;
                        int idx_x = (int)((x_local + random - spatial_x_min)
                                          / dspatial_x);
                        if (in_grid(idx_tau, N_tau) && in_grid(idx_x, N_xpt)) {
                            dNdtaudx1_array[idx_tau][idx_x]++;
                        }
                    }
                }
            }
            if (fabs(x_local) < 0.5) {
                if (tau_local > tau_min && tau_local < tau_max) {
                    double random =
                        ran_gen_ptr->rand_uniform() * intrinsic_dtau;
                    int idx_tau = (int)((tau_local + random - tau_min) / dtau);
                    if (y_local > spatial_x_min && y_local < spatial_x_max) {
                        double random =
                            ran_gen_ptr->rand_uniform() * intrinsic_dx;
                        int idx_x = (int)((y_local + random - spatial_x_min)
                                          / dspatial_x);
                        if (in_grid(idx_tau, N_tau) && in_grid(idx_x, N_xpt)) {
                            dNdtaudx2_array[idx_tau][idx_x]++;
                        }
                    }
                }
            }

            // third dN/deta_s
            double etas_local =
                (0.5 * log((t_local + z_local) / (t_local - z_local)));
            double y_minus_etas = rap_local - etas_local;
            if (y_minus_etas > eta_s_min && y_minus_etas < eta_s_max) {
                double random = ran_gen_ptr->rand_uniform() * intrinsic_detas;
                int idx = (int)((y_minus_etas + random - eta_s_min) / deta_s);
                if (in_grid(idx, N_eta_s)) {
                    eta_s_array[idx] += y_minus_etas;
                    dNdetas_array[idx]++;
                }
            }
        }
    }
}

void singleParticleSpectra::output_dNdSV() {
    // Reconstructed tracks carry weight 1/branching ratio. One factor for all
    // spatial histograms: exact when the sample is all reconstructed or all
    // stored, approximate for a mixture.
    const double yield_scale = spatial_count_ > 0.0
                                   ? 1.0
                                         + (reconstructed_weight_ - 1.0)
                                               * spatial_reconstructed_count_
                                               / spatial_count_
                                   : 1.0;
    // first dN/dtau
    ostringstream filename;
    filename << path_ << "/check_" << particle_monval << "_dNdtau"
         << output_suffix(spatial_output_);
    ofstream output(filename.str().c_str());
    for (int i = 0; i < N_tau; i++) {
        tau_array[i] = tau_array[i] / (dNdtau_array[i] + 1.);
        dNdtau_array[i] = dNdtau_array[i] / total_number_of_events;
        double dNdtau_err = sqrt(dNdtau_array[i] / total_number_of_events);
        dNdtau_array[i] = dNdtau_array[i] * yield_scale;
        dNdtau_err = dNdtau_err * yield_scale;
        output << scientific << setw(18) << setprecision(8) << tau_array[i]
               << "   " << dNdtau_array[i] / dtau << "   " << dNdtau_err / dtau
               << '\n';
    }
    output.close();

    // second dN/dx
    ostringstream filename2;
    filename2 << path_ << "/check_" << particle_monval << "_dNdx"
         << output_suffix(spatial_output_);
    ofstream output2(filename2.str().c_str());
    output2 << "# x  dN/dx  y  dN/dy  r  dN/dr" << '\n';
    for (int i = 0; i < N_xpt; i++) {
        xpt_array[i] = xpt_array[i] / (dNdx1_array[i] + 1.);
        ypt_array[i] = ypt_array[i] / (dNdx2_array[i] + 1.);
        rpt_array[i] = rpt_array[i] / (dNdr_array[i] + 1.);
        dNdx1_array[i] = dNdx1_array[i] / total_number_of_events;
        dNdx2_array[i] = dNdx2_array[i] / total_number_of_events;
        dNdr_array[i] = dNdr_array[i] / total_number_of_events;
        double dNdx1_err = sqrt(dNdx1_array[i] / total_number_of_events);
        double dNdx2_err = sqrt(dNdx2_array[i] / total_number_of_events);
        double dNdr_err = sqrt(dNdr_array[i] / total_number_of_events);
        dNdx1_array[i] = dNdx1_array[i] * yield_scale;
        dNdx2_array[i] = dNdx2_array[i] * yield_scale;
        dNdr_array[i] = dNdr_array[i] * yield_scale;
        dNdx1_err = dNdx1_err * yield_scale;
        dNdx2_err = dNdx2_err * yield_scale;
        dNdr_err = dNdr_err * yield_scale;
        output2 << scientific << setw(18) << setprecision(8) << xpt_array[i]
                << "   " << dNdx1_array[i] / dspatial_x << "   "
                << dNdx1_err / dspatial_x << "   " << ypt_array[i] << "   "
                << dNdx2_array[i] / dspatial_x << "   "
                << dNdx2_err / dspatial_x << "   " << rpt_array[i] << "   "
                << dNdr_array[i] / dspatial_r << "   " << dNdr_err / dspatial_r
                << '\n';
    }
    output2.close();

    // dN/(dtau dx)
    ostringstream filename2_1, filename2_2;
    filename2_1 << path_ << "/check_" << particle_monval << "_dNdtaudx1"
         << output_suffix(spatial_output_);
    filename2_2 << path_ << "/check_" << particle_monval << "_dNdtaudx2"
         << output_suffix(spatial_output_);
    ofstream output2_1(filename2_1.str().c_str());
    ofstream output2_2(filename2_2.str().c_str());
    for (int i = 0; i < N_tau; i++) {
        for (int j = 0; j < N_xpt; j++) {
            dNdtaudx1_array[i][j] =
                (dNdtaudx1_array[i][j] / total_number_of_events);
            dNdtaudx2_array[i][j] =
                (dNdtaudx2_array[i][j] / total_number_of_events);
            dNdtaudx1_array[i][j] = (dNdtaudx1_array[i][j] * yield_scale);
            dNdtaudx2_array[i][j] = (dNdtaudx2_array[i][j] * yield_scale);
            output2_1 << scientific << setw(18) << setprecision(8)
                      << dNdtaudx1_array[i][j] / dspatial_x / dtau << "   ";
            output2_2 << scientific << setw(18) << setprecision(8)
                      << dNdtaudx2_array[i][j] / dspatial_x / dtau << "   ";
        }
        output2_1 << '\n';
        output2_2 << '\n';
    }
    output2_1.close();
    output2_2.close();

    // third dN/detas
    ostringstream filename3;
    filename3 << path_ << "/check_" << particle_monval << "_dNdetas"
         << output_suffix(spatial_output_);
    ofstream output3(filename3.str().c_str());
    for (int i = 0; i < N_eta_s; i++) {
        eta_s_array[i] = eta_s_array[i] / (dNdetas_array[i] + 1.);
        dNdetas_array[i] = dNdetas_array[i] / total_number_of_events;
        double dNdetas_err = sqrt(dNdetas_array[i] / total_number_of_events);
        dNdetas_array[i] = dNdetas_array[i] * yield_scale;
        dNdetas_err = dNdetas_err * yield_scale;
        output3 << scientific << setw(18) << setprecision(8) << eta_s_array[i]
                << "   " << dNdetas_array[i] / deta_s << "   "
                << dNdetas_err / deta_s << '\n';
    }
    output3.close();
}
