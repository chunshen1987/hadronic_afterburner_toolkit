#include "HBT_correlation.h"

#include <stdlib.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "arsenal.h"
#include "parameters.h"

namespace {
// Track preselection for the HBT pair loops that removes no binned pair.
//
// Every binned pair has q_inv^2 < q_inv_max_sq. The 3D bins need
// |q_out|, |q_side|, |q_long| < Q, and q_inv^2 <= q_out^2 + q_side^2 + q_long^2
// in the lab or the LCMS, so q_inv^2 < 3 Q^2. The q_inv bins need q_inv < Q.
// For two tracks with transverse masses m_T1, m_T2 and masses m1, m2,
//   q_inv^2 = 2 p1.p2 - m1^2 - m2^2
//           >= 2 m_T1 m_T2 (cosh(y1 - y2) - 1) - (m1 - m2)^2,
// using p_T1.p_T2 <= p_T1 p_T2 and m_T1 m_T2 - p_T1 p_T2 >= m1 m2. The pair
// velocity K_z/K_E is an energy-weighted mean of the track velocities, so
// the pair rapidity lies between y1 and y2. A track at y can therefore
// enter a binned pair only if y_min - dy_max < y < y_max + dy_max, where
// dy_max follows from the q_inv bound with the smallest partner m_T.
struct TrackKinematics {
    double rapidity, mT, mass;
    bool valid;
};

TrackKinematics track_kinematics(const particle_info &particle) {
    TrackKinematics result {0.0, 0.0, 0.0, false};
    const double mT_sq = particle.E * particle.E - particle.pz * particle.pz;
    if (!(particle.E > std::abs(particle.pz)) || !(mT_sq > 0.0)) return result;
    result.mT = std::sqrt(mT_sq);
    result.rapidity =
        0.5 * std::log((particle.E + particle.pz) / (particle.E - particle.pz));
    const double pT_sq = particle.px * particle.px + particle.py * particle.py;
    // The bound needs a timelike four-momentum; otherwise keep every track.
    if (mT_sq < pT_sq) return result;
    result.mass = std::sqrt(mT_sq - pT_sq);
    result.valid = true;
    return result;
}

struct PartnerBounds {
    double mT_min = std::numeric_limits<double>::infinity();
    double mass_min = std::numeric_limits<double>::infinity();
    double mass_max = 0.0;
    bool usable = true;
};

void add_partner_bounds(
    const std::vector<particle_info> &tracks, PartnerBounds &bounds) {
    for (const auto &particle : tracks) {
        const TrackKinematics kin = track_kinematics(particle);
        if (!kin.valid) {
            bounds.usable = false;
            return;
        }
        bounds.mT_min = std::min(bounds.mT_min, kin.mT);
        bounds.mass_min = std::min(bounds.mass_min, kin.mass);
        bounds.mass_max = std::max(bounds.mass_max, kin.mass);
    }
}

// Returns false when every track must be kept. Otherwise fills selected with
// the tracks of `tracks` that can form a binned pair with any partner
// described by `partners`, in their original order.
bool select_pair_candidates(
    const std::vector<particle_info> &tracks, const PartnerBounds &partners,
    double mass_spread, double q_inv_max_sq, double rapidity_min,
    double rapidity_max, std::vector<particle_info> &selected) {
    if (!partners.usable || !(partners.mT_min > 0.0)) return false;
    // A relative margin of 1e-3 absorbs rounding in the bound.
    const double numerator =
        1.001 * (q_inv_max_sq + mass_spread * mass_spread) + 1e-12;
    selected.clear();
    selected.reserve(tracks.size());
    for (const auto &particle : tracks) {
        const TrackKinematics kin = track_kinematics(particle);
        if (!kin.valid) return false;
        const double dy_max =
            std::acosh(1.0 + numerator / (2.0 * kin.mT * partners.mT_min));
        if (kin.rapidity > rapidity_min - dy_max
            && kin.rapidity < rapidity_max + dy_max) {
            selected.push_back(particle);
        }
    }
    return true;
}
}  // namespace

HBT_correlation::HBT_correlation(
    ParameterReader &paraRdr, std::string path,
    std::shared_ptr<RandomUtil::Random> ran_gen)
    : paraRdr_(paraRdr), path_(path) {
    ran_gen_ptr_ = ran_gen;
    const double group_size = paraRdr.exist("number_of_oversample_events")
                                  ? paraRdr.getVal("number_of_oversample_events")
                                  : 1.0;
    if (!(group_size >= 1 && group_size <= std::numeric_limits<int>::max())
        || std::floor(group_size) != group_size) {
        throw std::runtime_error(
            "number_of_oversample_events must be a positive integer");
    }
    oversamples_per_collision_ = static_cast<int>(group_size);

    long_comoving_boost = false;
    if (paraRdr_.getVal("long_comoving_boost") == 1) long_comoving_boost = true;

    qnpts = paraRdr_.getVal("qnpts");
    q_min = paraRdr_.getVal("q_min");
    q_max = paraRdr_.getVal("q_max");
    delta_q = (q_max - q_min) / (qnpts - 1);

    for (int i = 0; i < qnpts; i++) {
        q_out.push_back(q_min + i * delta_q);
        q_side.push_back(q_min + i * delta_q);
        q_long.push_back(q_min + i * delta_q);
    }

    needed_number_of_pairs = paraRdr_.getVal("needed_number_of_pairs");

    azimuthal_flag_ = paraRdr_.getVal("azimuthal_flag");
    invariant_radius_flag_ = paraRdr_.getVal("invariant_radius_flag");

    n_KT = paraRdr_.getVal("n_KT");
    n_Kphi = paraRdr_.getVal("n_Kphi");
    KT_min = paraRdr_.getVal("KT_min");
    KT_max = paraRdr_.getVal("KT_max");
    Krap_min_ = paraRdr_.getVal("HBTrap_min");
    Krap_max_ = paraRdr_.getVal("HBTrap_max");

    dKT = (KT_max - KT_min) / (n_KT - 1);
    dKphi = 2 * M_PI / n_Kphi;  // does not need 0 and 2*pi

    for (int i = 0; i < n_KT; i++) {
        KT_array_.push_back(KT_min + i * dKT);
        number_of_pairs_numerator_KTdiff.push_back(0);
        number_of_pairs_denormenator_KTdiff.push_back(0);
        number_of_pairs_numerator_KTdiff_qinv_.push_back(0);
        number_of_pairs_denormenator_KTdiff_qinv_.push_back(0);
    }

    for (int i = 0; i < n_Kphi; i++) {
        Kphi_array_.push_back(i * dKphi);
    }

    if (azimuthal_flag_ == 1) {
        create_a_2D_array(number_of_pairs_numerator_KTKphidiff, n_KT, n_Kphi);
        create_a_2D_array(
            number_of_pairs_denormenator_KTKphidiff, n_KT, n_Kphi);
        for (int i = 0; i < n_KT; i++) {
            for (int j = 0; j < n_Kphi; j++) {
                number_of_pairs_numerator_KTKphidiff[i][j] = 0;
                number_of_pairs_denormenator_KTKphidiff[i][j] = 0;
            }
        }
    }

    number_pairs_num = 0;
    number_pairs_denorm = 0;
    psi_ref = 0.;

    if (invariant_radius_flag_ == 1) {
        create_a_2D_array(q_inv_mean, n_KT, qnpts);
        create_a_2D_array(correl_1d_inv_num, n_KT, qnpts);
        create_a_2D_array(correl_1d_inv_num_count, n_KT, qnpts);
        create_a_2D_array(correl_1d_inv_denorm, n_KT, qnpts);
        for (int iK = 0; iK < n_KT; iK++) {
            for (int k = 0; k < qnpts; k++) {
                q_inv_mean[iK][k] = 0.0;
                correl_1d_inv_num[iK][k] = 0.0;
                correl_1d_inv_num_count[iK][k] = 0.0;
                correl_1d_inv_denorm[iK][k] = 0.0;
            }
        }
    }

    if (azimuthal_flag_ == 0) {
        create_a_4D_array(q_out_mean, n_KT, qnpts, qnpts, qnpts);
        create_a_4D_array(q_side_mean, n_KT, qnpts, qnpts, qnpts);
        create_a_4D_array(q_long_mean, n_KT, qnpts, qnpts, qnpts);
        create_a_4D_array(correl_3d_num, n_KT, qnpts, qnpts, qnpts);
        create_a_4D_array(correl_3d_num_count, n_KT, qnpts, qnpts, qnpts);
        create_a_4D_array(correl_3d_denorm, n_KT, qnpts, qnpts, qnpts);
        for (int iK = 0; iK < n_KT; iK++) {
            for (int i = 0; i < qnpts; i++) {
                for (int j = 0; j < qnpts; j++) {
                    for (int k = 0; k < qnpts; k++) {
                        q_out_mean[iK][i][j][k] = 0.0;
                        q_side_mean[iK][i][j][k] = 0.0;
                        q_long_mean[iK][i][j][k] = 0.0;
                        correl_3d_num[iK][i][j][k] = 0.0;
                        correl_3d_num_count[iK][i][j][k] = 0.0;
                        correl_3d_denorm[iK][i][j][k] = 0.0;
                    }
                }
            }
        }
    } else {
        create_a_5D_array(q_out_diff_mean, n_KT, n_Kphi, qnpts, qnpts, qnpts);
        create_a_5D_array(q_side_diff_mean, n_KT, n_Kphi, qnpts, qnpts, qnpts);
        create_a_5D_array(q_long_diff_mean, n_KT, n_Kphi, qnpts, qnpts, qnpts);
        create_a_5D_array(
            correl_3d_Kphi_diff_num, n_KT, n_Kphi, qnpts, qnpts, qnpts);
        create_a_5D_array(
            correl_3d_Kphi_diff_num_count, n_KT, n_Kphi, qnpts, qnpts, qnpts);
        create_a_5D_array(
            correl_3d_Kphi_diff_denorm, n_KT, n_Kphi, qnpts, qnpts, qnpts);
        for (int iK = 0; iK < n_KT; iK++) {
            for (int iphi = 0; iphi < n_Kphi; iphi++) {
                for (int i = 0; i < qnpts; i++) {
                    for (int j = 0; j < qnpts; j++) {
                        for (int k = 0; k < qnpts; k++) {
                            q_out_diff_mean[iK][iphi][i][j][k] = 0.0;
                            q_side_diff_mean[iK][iphi][i][j][k] = 0.0;
                            q_long_diff_mean[iK][iphi][i][j][k] = 0.0;
                            correl_3d_Kphi_diff_num[iK][iphi][i][j][k] = 0.0;
                            correl_3d_Kphi_diff_num_count[iK][iphi][i][j][k] =
                                0.0;
                            correl_3d_Kphi_diff_denorm[iK][iphi][i][j][k] = 0.0;
                        }
                    }
                }
            }
        }
    }
}

HBT_correlation::~HBT_correlation() {
    if (invariant_radius_flag_ == 1) {
        delete_a_2D_array(q_inv_mean, n_KT);
        delete_a_2D_array(correl_1d_inv_num, n_KT);
        delete_a_2D_array(correl_1d_inv_num_count, n_KT);
        delete_a_2D_array(correl_1d_inv_denorm, n_KT);
    }

    if (azimuthal_flag_ == 0) {
        delete_a_4D_array(q_out_mean, n_KT, qnpts, qnpts);
        delete_a_4D_array(q_side_mean, n_KT, qnpts, qnpts);
        delete_a_4D_array(q_long_mean, n_KT, qnpts, qnpts);
        delete_a_4D_array(correl_3d_num, n_KT, qnpts, qnpts);
        delete_a_4D_array(correl_3d_num_count, n_KT, qnpts, qnpts);
        delete_a_4D_array(correl_3d_denorm, n_KT, qnpts, qnpts);
    } else {
        delete_a_5D_array(q_out_diff_mean, n_KT, n_Kphi, qnpts, qnpts);
        delete_a_5D_array(q_side_diff_mean, n_KT, n_Kphi, qnpts, qnpts);
        delete_a_5D_array(q_long_diff_mean, n_KT, n_Kphi, qnpts, qnpts);
        delete_a_5D_array(correl_3d_Kphi_diff_num, n_KT, n_Kphi, qnpts, qnpts);
        delete_a_5D_array(
            correl_3d_Kphi_diff_num_count, n_KT, n_Kphi, qnpts, qnpts);
        delete_a_5D_array(
            correl_3d_Kphi_diff_denorm, n_KT, n_Kphi, qnpts, qnpts);
    }

    if (azimuthal_flag_ == 1) {
        delete_a_2D_array(number_of_pairs_numerator_KTKphidiff, n_KT);
        delete_a_2D_array(number_of_pairs_denormenator_KTKphidiff, n_KT);
    }
}

void HBT_correlation::calculate_HBT_correlation_function(
    std::shared_ptr<particleSamples> particle_list_in) {
    set_particle_list(particle_list_in);
    const int nev = particle_list->get_number_of_events();
    if (nev == 0) return;
    const bool use_external = particle_list->uses_external_mixed_events();
    const int mixed_nev =
        use_external ? particle_list->get_number_of_mixed_events() : 0;
    if (use_external && mixed_nev == 0) {
        throw std::runtime_error("HBT analysis requires a mixed-event batch");
    }

    for (int iev = 0; iev < nev; iev++) {
        if (use_external && pending_events_ == 0) {
            for (int i = 0; i < mixed_nev / 2 + 1; i++) {
                const int ref = ran_gen_ptr_->rand_int_uniform() % mixed_nev;
                append_external_reference(ref, external_reference_);
            }
        }
        const int count = particle_list->get_number_of_particles(iev);
        for (int i = 0; i < count; i++) {
            pending_collision_.push_back(particle_list->get_particle(iev, i));
        }
        if (++pending_events_ == oversamples_per_collision_) finish_collision();
    }
}

void HBT_correlation::finish_collision() {
    if (pending_events_ == 0) return;
    if (azimuthal_flag_ == 1) set_event_plane(pending_collision_, 2);
    bin_same_event_pairs(pending_collision_);

    if (particle_list->uses_external_mixed_events()) {
        bin_mixed_event_pairs(pending_collision_, external_reference_);
        external_reference_.clear();
    } else {
        if (!previous_collision_.empty()) {
            rotate_reference(previous_collision_, 0);
            bin_mixed_event_pairs(pending_collision_, previous_collision_);
        }
        previous_collision_.swap(pending_collision_);
    }
    pending_collision_.clear();
    pending_events_ = 0;
}

void HBT_correlation::output_HBTcorrelation() {
    finish_collision();
    if (invariant_radius_flag_ == 1) output_correlation_function_inv();

    if (azimuthal_flag_ == 0) {
        output_correlation_function();
    } else {
        output_correlation_function_Kphi_differential();
    }
}

//! This function computes the flow event plane angle Psi_n
//! The results is stored in the variable psi_ref
//! This is needed when computing the azimuthal dependent HBT radii
void HBT_correlation::calculate_flow_event_plane_angle(int n_order) {
    std::vector<particle_info> particles;
    for (int iev = 0; iev < particle_list->get_number_of_events(); iev++) {
        const int count = particle_list->get_number_of_particles(iev);
        for (int i = 0; i < count; i++) {
            particles.push_back(particle_list->get_particle(iev, i));
        }
    }
    set_event_plane(particles, n_order);
}

void HBT_correlation::set_event_plane(
    const std::vector<particle_info> &particles, int harmonic) {
    double real = 0.0, imag = 0.0;
    for (const auto &particle : particles) {
        const double angle = harmonic * std::atan2(particle.py, particle.px);
        real += std::cos(angle);
        imag += std::sin(angle);
    }
    psi_ref = std::atan2(imag, real) / harmonic;
}

void HBT_correlation::combine_and_bin_particle_pairs(
    std::vector<int> event_list) {
    std::vector<particle_info> particles;
    for (int event : event_list) {
        const int count = particle_list->get_number_of_particles(event);
        for (int i = 0; i < count; i++) {
            particles.push_back(particle_list->get_particle(event, i));
        }
    }
    bin_same_event_pairs(particles);
}

double HBT_correlation::binned_q_inv_max_sq() const {
    const double q_bound = std::max(
        std::abs(q_min - 0.5 * delta_q), std::abs(q_max + 0.5 * delta_q));
    return 3.0 * q_bound * q_bound;
}

void HBT_correlation::bin_same_event_pairs(
    const std::vector<particle_info> &all_tracks) {
    PartnerBounds bounds;
    add_partner_bounds(all_tracks, bounds);
    std::vector<particle_info> candidates;
    const bool reduced = select_pair_candidates(
        all_tracks, bounds, bounds.mass_max - bounds.mass_min,
        binned_q_inv_max_sq(), Krap_min_, Krap_max_, candidates);
    const std::vector<particle_info> &temp_particle_list =
        reduced ? candidates : all_tracks;
    double hbarC_inv = 1. / hbarC;
    const double min_pair_beta = std::tanh(Krap_min_);
    const double max_pair_beta = std::tanh(Krap_max_);
    const auto number_of_particles = temp_particle_list.size();
    double KT_min_sq = KT_min * KT_min;
    double KT_max_sq = KT_max * KT_max;
    for (std::size_t i = 0; i < number_of_particles; i++) {
        double particle_1_px = temp_particle_list[i].px;
        double particle_1_py = temp_particle_list[i].py;
        double particle_1_pz = temp_particle_list[i].pz;
        double particle_1_E = temp_particle_list[i].E;
        double particle_1_x = temp_particle_list[i].x;
        double particle_1_y = temp_particle_list[i].y;
        double particle_1_z = temp_particle_list[i].z;
        double particle_1_t = temp_particle_list[i].t;

        for (std::size_t j = i + 1; j < number_of_particles; j++) {
            double particle_2_px = temp_particle_list[j].px;
            double particle_2_py = temp_particle_list[j].py;
            double particle_2_pz = temp_particle_list[j].pz;
            double particle_2_E = temp_particle_list[j].E;
            double particle_2_x = temp_particle_list[j].x;
            double particle_2_y = temp_particle_list[j].y;
            double particle_2_z = temp_particle_list[j].z;
            double particle_2_t = temp_particle_list[j].t;

            double K_z = 0.5 * (particle_1_pz + particle_2_pz);
            double K_E = 0.5 * (particle_1_E + particle_2_E);
            if (!(K_E > std::abs(K_z))) continue;
            const double K_z_over_K_E = K_z / K_E;
            if (!(K_z_over_K_E > min_pair_beta && K_z_over_K_E < max_pair_beta))
                continue;

            // K_T cut
            double K_x = 0.5 * (particle_1_px + particle_2_px);
            double K_y = 0.5 * (particle_1_py + particle_2_py);
            double K_perp_sq = K_x * K_x + K_y * K_y;
            if (K_perp_sq < KT_min_sq || K_perp_sq >= KT_max_sq) {
                continue;
            }

            double K_perp = sqrt(K_perp_sq);
            int Kperp_idx = static_cast<int>((K_perp - KT_min) / dKT);

            // calculate q_inv
            double q_x = particle_1_px - particle_2_px;
            double q_y = particle_1_py - particle_2_py;
            double q_z = particle_1_pz - particle_2_pz;
            double q_E = particle_1_E - particle_2_E;
            double local_q_inv =
                sqrt(-(q_E * q_E - q_x * q_x - q_y * q_y - q_z * q_z));

            double t_diff = particle_1_t - particle_2_t;
            double x_diff = particle_1_x - particle_2_x;
            double y_diff = particle_1_y - particle_2_y;
            double z_diff = particle_1_z - particle_2_z;

            if (invariant_radius_flag_ == 1) {
                if (number_of_pairs_numerator_KTdiff_qinv_[Kperp_idx]
                    < 50 * needed_number_of_pairs) {
                    if (local_q_inv > (q_min - delta_q / 2. + 1e-8)
                        && local_q_inv < (q_max + delta_q / 2. - 1e-8)) {
                        int qinv_idx = static_cast<int>(
                            (local_q_inv - (q_min - delta_q / 2.)) / delta_q);
                        number_of_pairs_numerator_KTdiff_qinv_[Kperp_idx]++;
                        double cos_qx =
                            cos(hbarC_inv
                                * (q_E * t_diff - q_x * x_diff - q_y * y_diff
                                   - q_z * z_diff));
                        correl_1d_inv_num_count[Kperp_idx][qinv_idx]++;
                        q_inv_mean[Kperp_idx][qinv_idx] += local_q_inv;
                        correl_1d_inv_num[Kperp_idx][qinv_idx] += cos_qx;
                    }
                }
            }

            // calculate qout and qside in lcms
            if (K_perp == 0.0) continue;
            double cos_K_phi = K_x / K_perp;
            double sin_K_phi = K_y / K_perp;

            double local_q_out = q_x * cos_K_phi + q_y * sin_K_phi;
            if (local_q_out < (q_min - delta_q / 2. + 1e-8)
                || local_q_out > (q_max + delta_q / 2. - 1e-8)) {
                continue;
            }

            int qout_idx = static_cast<int>(
                (local_q_out - (q_min - delta_q / 2.)) / delta_q);
            if (qout_idx >= qnpts) continue;

            double local_q_side = q_y * cos_K_phi - q_x * sin_K_phi;
            if (local_q_side < (q_min - delta_q / 2. + 1e-8)
                || local_q_side > (q_max + delta_q / 2. - 1e-8)) {
                continue;
            }

            int qside_idx = static_cast<int>(
                (local_q_side - (q_min - delta_q / 2.)) / delta_q);
            if (qside_idx >= qnpts) continue;

            double local_q_long = q_z;
            if (long_comoving_boost) {
                // calcualte qlong in the lcms
                double Mt = sqrt(K_E * K_E - K_z * K_z);
                double boost_gamma = K_E / Mt;
                double boost_beta = K_z_over_K_E;
                // boost qz to lcms
                local_q_long = boost_gamma * (q_z - boost_beta * q_E);
            }

            if (local_q_long < (q_min - delta_q / 2. + 1e-8)
                || local_q_long > (q_max + delta_q / 2. - 1e-8)) {
                continue;
            }
            int qlong_idx = static_cast<int>(
                (local_q_long - (q_min - delta_q / 2.)) / delta_q);
            if (qlong_idx >= qnpts) continue;

            int Kphi_idx;
            if (azimuthal_flag_ == 0) {
                if (number_of_pairs_numerator_KTdiff[Kperp_idx]
                    > needed_number_of_pairs) {
                    continue;
                }
                number_of_pairs_numerator_KTdiff[Kperp_idx]++;
            } else {
                double local_K_phi = atan2(K_y, K_x);
                double delta_phi = local_K_phi - psi_ref;
                while (delta_phi < 0.) {
                    delta_phi += 2. * M_PI;
                }
                while (delta_phi > 2. * M_PI) {
                    delta_phi -= 2. * M_PI;
                }
                Kphi_idx = static_cast<int>(delta_phi / dKphi);
                if (Kphi_idx < 0 || Kphi_idx >= n_Kphi) {
                    messager << "delta_phi = " << delta_phi
                             << " is out of bound of [0, 2pi]! "
                             << "Kphi_idx = " << Kphi_idx;
                    messager.flush("warning");
                    continue;
                }
                if (number_of_pairs_numerator_KTKphidiff[Kperp_idx][Kphi_idx]
                    > needed_number_of_pairs) {
                    continue;
                }
                number_of_pairs_numerator_KTKphidiff[Kperp_idx][Kphi_idx]++;
            }

            double cos_qx = cos(
                hbarC_inv
                * (q_E * t_diff - q_x * x_diff - q_y * y_diff - q_z * z_diff));

            // bin results
            if (azimuthal_flag_ == 0) {
                correl_3d_num_count[Kperp_idx][qout_idx][qside_idx]
                                   [qlong_idx]++;
                q_out_mean[Kperp_idx][qout_idx][qside_idx][qlong_idx] +=
                    local_q_out;
                q_side_mean[Kperp_idx][qout_idx][qside_idx][qlong_idx] +=
                    local_q_side;
                q_long_mean[Kperp_idx][qout_idx][qside_idx][qlong_idx] +=
                    local_q_long;
                correl_3d_num[Kperp_idx][qout_idx][qside_idx][qlong_idx] +=
                    cos_qx;
            } else {
                correl_3d_Kphi_diff_num_count[Kperp_idx][Kphi_idx][qout_idx]
                                             [qside_idx][qlong_idx]++;
                q_out_diff_mean[Kperp_idx][Kphi_idx][qout_idx][qside_idx]
                               [qlong_idx] += local_q_out;
                q_side_diff_mean[Kperp_idx][Kphi_idx][qout_idx][qside_idx]
                                [qlong_idx] += local_q_side;
                q_long_diff_mean[Kperp_idx][Kphi_idx][qout_idx][qside_idx]
                                [qlong_idx] += local_q_long;
                correl_3d_Kphi_diff_num[Kperp_idx][Kphi_idx][qout_idx]
                                       [qside_idx][qlong_idx] += cos_qx;
            }
        }
    }
}

void HBT_correlation::rotate_reference(
    std::vector<particle_info> &reference, std::size_t first) {
    const double angle = ran_gen_ptr_->rand_uniform() * 2. * M_PI;
    const double cosine = std::cos(angle), sine = std::sin(angle);
    for (std::size_t i = first; i < reference.size(); i++) {
        auto &particle = reference[i];
        const double px = particle.px;
        particle.px = px * cosine - particle.py * sine;
        particle.py = px * sine + particle.py * cosine;
    }
}

void HBT_correlation::append_external_reference(
    int event_id, std::vector<particle_info> &reference) {
    const auto start = reference.size();
    const int count = particle_list->get_number_of_particles_mixed_event(event_id);
    for (int i = 0; i < count; i++) {
        reference.push_back(particle_list->get_particle_from_mixed_event(event_id, i));
    }
    rotate_reference(reference, start);
}

void HBT_correlation::combine_and_bin_particle_pairs_mixed_events(
    int event_id, std::vector<int> mixed_event_list) {
    std::vector<particle_info> particles, reference;
    const int count = particle_list->get_number_of_particles(event_id);
    for (int i = 0; i < count; i++) {
        particles.push_back(particle_list->get_particle(event_id, i));
    }
    for (int ref : mixed_event_list) append_external_reference(ref, reference);
    bin_mixed_event_pairs(particles, reference);
}

void HBT_correlation::bin_mixed_event_pairs(
    const std::vector<particle_info> &all_tracks_1,
    const std::vector<particle_info> &all_tracks_2) {
    PartnerBounds bounds_1, bounds_2;
    add_partner_bounds(all_tracks_1, bounds_1);
    add_partner_bounds(all_tracks_2, bounds_2);
    const double mass_spread = std::max(bounds_1.mass_max, bounds_2.mass_max)
                               - std::min(bounds_1.mass_min, bounds_2.mass_min);
    const double q_inv_max_sq = binned_q_inv_max_sq();
    std::vector<particle_info> candidates_1, candidates_2;
    const bool reduced = select_pair_candidates(
                             all_tracks_1, bounds_2, mass_spread, q_inv_max_sq,
                             Krap_min_, Krap_max_, candidates_1)
                         && select_pair_candidates(
                             all_tracks_2, bounds_1, mass_spread, q_inv_max_sq,
                             Krap_min_, Krap_max_, candidates_2);
    const std::vector<particle_info> &temp_particle_list_1 =
        reduced ? candidates_1 : all_tracks_1;
    const std::vector<particle_info> &temp_particle_list_2 =
        reduced ? candidates_2 : all_tracks_2;
    const auto number_of_particles_1 = temp_particle_list_1.size();
    const auto number_of_particles_2 = temp_particle_list_2.size();
    const double min_pair_beta = std::tanh(Krap_min_);
    const double max_pair_beta = std::tanh(Krap_max_);
    double KT_min_sq = KT_min * KT_min;
    double KT_max_sq = KT_max * KT_max;
    for (std::size_t i = 0; i < number_of_particles_1; i++) {
        double particle_1_pz = temp_particle_list_1[i].pz;
        double particle_1_px = temp_particle_list_1[i].px;
        double particle_1_py = temp_particle_list_1[i].py;
        double particle_1_E = temp_particle_list_1[i].E;
        for (std::size_t j = 0; j < number_of_particles_2; j++) {
            double particle_2_pz = temp_particle_list_2[j].pz;
            double particle_2_px = temp_particle_list_2[j].px;
            double particle_2_py = temp_particle_list_2[j].py;
            double particle_2_E = temp_particle_list_2[j].E;

            double K_z = 0.5 * (particle_1_pz + particle_2_pz);
            double K_E = 0.5 * (particle_1_E + particle_2_E);
            if (!(K_E > std::abs(K_z))) continue;
            const double K_z_over_K_E = K_z / K_E;
            if (!(K_z_over_K_E > min_pair_beta && K_z_over_K_E < max_pair_beta))
                continue;

            double K_x = 0.5 * (particle_1_px + particle_2_px);
            double K_y = 0.5 * (particle_1_py + particle_2_py);
            double K_perp_sq = K_x * K_x + K_y * K_y;
            if (K_perp_sq < KT_min_sq || K_perp_sq >= KT_max_sq) {
                continue;
            }

            double K_perp = sqrt(K_perp_sq);
            int Kperp_idx = static_cast<int>((K_perp - KT_min) / dKT);

            double q_z = particle_1_pz - particle_2_pz;
            double q_E = particle_1_E - particle_2_E;
            double q_x = particle_1_px - particle_2_px;
            double q_y = particle_1_py - particle_2_py;
            double local_q_inv =
                sqrt(-(q_E * q_E - q_x * q_x - q_y * q_y - q_z * q_z));

            if (invariant_radius_flag_ == 1) {
                if (local_q_inv > (q_min - delta_q / 2. + 1e-8)
                    && local_q_inv < (q_max + delta_q / 2. - 1e-8)) {
                    int qinv_idx = static_cast<int>(
                        (local_q_inv - (q_min - delta_q / 2.)) / delta_q);
                    if (qinv_idx < qnpts
                        && number_of_pairs_denormenator_KTdiff_qinv_[Kperp_idx]
                               < 50 * needed_number_of_pairs) {
                        number_of_pairs_denormenator_KTdiff_qinv_[Kperp_idx]++;
                        correl_1d_inv_denorm[Kperp_idx][qinv_idx] += 1.0;
                    }
                }
            }

            if (K_perp == 0.0) continue;
            double cos_K_phi = K_x / K_perp;
            double sin_K_phi = K_y / K_perp;
            double local_q_out = q_x * cos_K_phi + q_y * sin_K_phi;
            if (local_q_out < (q_min - delta_q / 2. + 1e-8)
                || local_q_out >= (q_max + delta_q / 2. - 1e-8)) {
                continue;
            }
            int qout_idx = static_cast<int>(
                (local_q_out - (q_min - delta_q / 2.)) / delta_q);
            if (qout_idx >= qnpts) continue;

            double local_q_side = q_y * cos_K_phi - q_x * sin_K_phi;
            if (local_q_side < (q_min - delta_q / 2. + 1e-8)
                || local_q_side >= (q_max + delta_q / 2. - 1e-8)) {
                continue;
            }

            int qside_idx = static_cast<int>(
                (local_q_side - (q_min - delta_q / 2.)) / delta_q);
            if (qside_idx >= qnpts) continue;

            double local_q_long = q_z;
            if (long_comoving_boost) {
                // calcualte qlong in the lcms
                double Mt = sqrt(K_E * K_E - K_z * K_z);
                double boost_gamma = K_E / Mt;
                double boost_beta = K_z_over_K_E;
                // boost qz to lcms
                local_q_long = boost_gamma * (q_z - boost_beta * q_E);
            }

            if (local_q_long < (q_min - delta_q / 2. + 1e-8)
                || local_q_long >= (q_max + delta_q / 2. - 1e-8)) {
                continue;
            }

            int qlong_idx = static_cast<int>(
                (local_q_long - (q_min - delta_q / 2.)) / delta_q);
            if (qlong_idx >= qnpts) continue;

            int Kphi_idx;
            if (azimuthal_flag_ == 0) {
                if (number_of_pairs_denormenator_KTdiff[Kperp_idx]
                    > needed_number_of_pairs) {
                    continue;
                }
                number_of_pairs_denormenator_KTdiff[Kperp_idx]++;
            } else {
                double local_K_phi = atan2(K_y, K_x);
                double delta_phi = local_K_phi - psi_ref;
                while (delta_phi < 0.) {
                    delta_phi += 2. * M_PI;
                }
                while (delta_phi > 2. * M_PI) {
                    delta_phi -= 2. * M_PI;
                }
                Kphi_idx = static_cast<int>(delta_phi / dKphi);
                if (Kphi_idx < 0 || Kphi_idx >= n_Kphi) {
                    messager << "delta_phi = " << delta_phi
                             << " is out of bound of [0, 2pi]! "
                             << "Kphi_idx = " << Kphi_idx;
                    messager.flush("warning");
                    continue;
                }
                if (number_of_pairs_denormenator_KTKphidiff[Kperp_idx][Kphi_idx]
                    > needed_number_of_pairs) {
                    continue;
                }
                number_of_pairs_denormenator_KTKphidiff[Kperp_idx][Kphi_idx]++;
            }

            // bin results
            if (azimuthal_flag_ == 0) {
                correl_3d_denorm[Kperp_idx][qout_idx][qside_idx][qlong_idx] +=
                    1.0;
            } else {
                correl_3d_Kphi_diff_denorm[Kperp_idx][Kphi_idx][qout_idx]
                                          [qside_idx][qlong_idx] += 1.0;
            }
        }
    }
}

void HBT_correlation::output_correlation_function_inv() {
    for (int iK = 0; iK < n_KT - 1; iK++) {
        double npair_ratio =
            (static_cast<double>(number_of_pairs_numerator_KTdiff_qinv_[iK])
             / static_cast<double>(
                 number_of_pairs_denormenator_KTdiff_qinv_[iK]));
        std::ostringstream filename;
        filename << path_ << "/HBT_correlation_function_inv_KT_"
                 << KT_array_[iK] << "_" << KT_array_[iK + 1] << ".dat";
        std::ofstream output(filename.str().c_str());
        for (int iqinv = 0; iqinv < qnpts; iqinv++) {
            double q_inv_local =
                (q_inv_mean[iK][iqinv] / correl_1d_inv_num_count[iK][iqinv]);
            double correl_fun_num = correl_1d_inv_num[iK][iqinv];
            double correl_fun_denorm =
                (correl_1d_inv_denorm[iK][iqinv] * npair_ratio);

            if (paraRdr_.getVal("ecoOutput", 0) == 1) {
                output << std::scientific << std::setw(18)
                       << std::setprecision(8) << correl_fun_num << "    "
                       << correl_fun_denorm << std::endl;
            } else {
                output << std::scientific << std::setw(18)
                       << std::setprecision(8) << q_inv_local << "    "
                       << correl_fun_num << "    " << correl_fun_denorm
                       << std::endl;
            }
        }
        output.close();
    }
}

void HBT_correlation::output_correlation_function() {
    for (int iK = 0; iK < n_KT - 1; iK++) {
        std::ostringstream suffix;
        suffix << KT_array_[iK] << "_" << KT_array_[iK + 1] << ".dat";
        output_3d_correlation(
            suffix.str(), number_of_pairs_numerator_KTdiff[iK],
            number_of_pairs_denormenator_KTdiff[iK], q_out_mean[iK],
            q_side_mean[iK], q_long_mean[iK], correl_3d_num_count[iK],
            correl_3d_num[iK], correl_3d_denorm[iK]);
    }
}

void HBT_correlation::output_correlation_function_Kphi_differential() {
    for (int iK = 0; iK < n_KT - 1; iK++) {
        for (int iKphi = 0; iKphi < n_Kphi; iKphi++) {
            std::ostringstream suffix;
            suffix << KT_array_[iK] << "_" << KT_array_[iK + 1] << "_Kphi_"
                   << Kphi_array_[iKphi] << ".dat";
            output_3d_correlation(
                suffix.str(), number_of_pairs_numerator_KTKphidiff[iK][iKphi],
                number_of_pairs_denormenator_KTKphidiff[iK][iKphi],
                q_out_diff_mean[iK][iKphi], q_side_diff_mean[iK][iKphi],
                q_long_diff_mean[iK][iKphi],
                correl_3d_Kphi_diff_num_count[iK][iKphi],
                correl_3d_Kphi_diff_num[iK][iKphi],
                correl_3d_Kphi_diff_denorm[iK][iKphi]);
        }
    }
}

void HBT_correlation::output_3d_correlation(
    const std::string &suffix, double total_same, double total_mixed,
    double ***out_sum, double ***side_sum, double ***long_sum,
    double ***same_count, double ***cos_sum, double ***mixed_count) {
    const bool compact = paraRdr_.getVal("ecoOutput", 0) == 1;
    const bool write_raw = paraRdr_.getVal("hbt_raw_output", 1) == 1;
    const double pair_ratio = total_mixed > 0.0 ? total_same / total_mixed : 0.0;
    std::ofstream output(path_ + "/HBT_correlation_function_KT_" + suffix);
    if (!output) throw std::runtime_error("Cannot write HBT output " + suffix);
    output << std::scientific << std::setprecision(8);

    std::ofstream raw;
    if (write_raw) {
        raw.open(path_ + "/HBT_raw_correlation_function_KT_" + suffix);
        if (!raw) {
            throw std::runtime_error("Cannot write raw HBT output " + suffix);
        }
        raw << "# q_out q_side q_long same_pairs sum_cos mixed_pairs\n"
            << std::scientific << std::setprecision(16);
    }

    for (int ilong = 0; ilong < qnpts; ilong++) {
        for (int iout = 0; iout < qnpts; iout++) {
            for (int iside = 0; iside < qnpts; iside++) {
                const double same = same_count[iout][iside][ilong];
                const double mixed = mixed_count[iout][iside][ilong];
                const double numerator = cos_sum[iout][iside][ilong];
                if (write_raw) {
                    raw << q_out[iout] << ' ' << q_side[iside] << ' '
                        << q_long[ilong] << ' ' << same << ' ' << numerator
                        << ' ' << mixed << '\n';
                }

                const bool populated = same >= 2.0 && mixed >= 2.0;
                if (!compact) {
                    const double out =
                        populated ? out_sum[iout][iside][ilong] / same
                                  : q_out[iout];
                    const double side =
                        populated ? side_sum[iout][iside][ilong] / same
                                  : q_side[iside];
                    const double longitudinal =
                        populated ? long_sum[iout][iside][ilong] / same
                                  : q_long[ilong];
                    output << std::setw(18) << out << "    " << side << "    "
                           << longitudinal << "    ";
                }
                output << (populated ? numerator : 0.0) << "    "
                       << (populated ? pair_ratio * mixed : mixed) << '\n';
            }
        }
    }
    output.close();
    if (write_raw) raw.close();
    if (!output || (write_raw && !raw)) {
        throw std::runtime_error("Failed to finish HBT output " + suffix);
    }
}

template <typename T>
void HBT_correlation::create_a_2D_array(T **&arr2D, int nx, int ny) {
    arr2D = new T *[nx];
    for (int i = 0; i < nx; i++) {
        arr2D[i] = new T[ny];
    }
}

template <typename T>
void HBT_correlation::delete_a_2D_array(T **&arr2D, int nx) {
    for (int i = 0; i < nx; i++) delete[] arr2D[i];
    delete[] arr2D;
}

template <typename T>
void HBT_correlation::create_a_3D_array(T ***&arr3D, int nx, int ny, int nz) {
    arr3D = new T **[nx];
    for (int i = 0; i < nx; i++) {
        create_a_2D_array(arr3D[i], ny, nz);
    }
}

template <typename T>
void HBT_correlation::delete_a_3D_array(T ***&arr3D, int nx, int ny) {
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) delete[] arr3D[i][j];
        delete[] arr3D[i];
    }
    delete[] arr3D;
}

template <typename T>
void HBT_correlation::create_a_4D_array(
    T ****&arr4D, int n1, int n2, int n3, int n4) {
    arr4D = new T ***[n1];
    for (int i = 0; i < n1; i++) {
        create_a_3D_array(arr4D[i], n2, n3, n4);
    }
}

template <typename T>
void HBT_correlation::delete_a_4D_array(T ****&arr4D, int n1, int n2, int n3) {
    for (int i = 0; i < n1; i++) {
        for (int j = 0; j < n2; j++) {
            for (int k = 0; k < n3; k++) delete[] arr4D[i][j][k];
            delete[] arr4D[i][j];
        }
        delete[] arr4D[i];
    }
    delete[] arr4D;
}

template <typename T>
void HBT_correlation::create_a_5D_array(
    T *****&arr5D, int n1, int n2, int n3, int n4, int n5) {
    arr5D = new T ****[n1];
    for (int i = 0; i < n1; i++) {
        create_a_4D_array(arr5D[i], n2, n3, n4, n5);
    }
}

template <typename T>
void HBT_correlation::delete_a_5D_array(
    T *****&arr5D, int n1, int n2, int n3, int n4) {
    for (int i = 0; i < n1; i++) {
        for (int j = 0; j < n2; j++) {
            for (int k = 0; k < n3; k++) {
                for (int l = 0; l < n4; l++) delete[] arr5D[i][j][k][l];
                delete[] arr5D[i][j][k];
            }
            delete[] arr5D[i][j];
        }
        delete[] arr5D[i];
    }
    delete[] arr5D;
}
