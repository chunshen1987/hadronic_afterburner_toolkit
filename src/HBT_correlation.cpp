#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <vector>
#include <iomanip>
#include <string>
#include "arsenal.h"
#include "parameters.h"
#include "HBT_correlation.h"

HBT_correlation::HBT_correlation(
            ParameterReader &paraRdr, std::string path,
            std::shared_ptr<RandomUtil::Random> ran_gen) :
        paraRdr_(paraRdr), path_(path) {

    ran_gen_ptr_ = ran_gen;

    qnpts   = paraRdr_.getVal("qnpts");
    q_min   = paraRdr_.getVal("q_min");
    q_max   = paraRdr_.getVal("q_max");
    delta_q = (q_max - q_min)/(qnpts - 1);

    for (int i = 0; i < qnpts; i++) {
        q_out.push_back(q_min + i*delta_q);
        q_side.push_back(q_min + i*delta_q);
        q_long.push_back(q_min + i*delta_q);
    }

    needed_number_of_pairs = paraRdr_.getVal("needed_number_of_pairs");

    azimuthal_flag_ = paraRdr_.getVal("azimuthal_flag");
    invariant_radius_flag_ = paraRdr_.getVal("invariant_radius_flag");

    n_KT            = paraRdr_.getVal("n_KT");
    n_Kphi          = paraRdr_.getVal("n_Kphi");
    KT_min          = paraRdr_.getVal("KT_min");
    KT_max          = paraRdr_.getVal("KT_max");
    Krap_min        = paraRdr_.getVal("Krap_min");
    Krap_max        = paraRdr_.getVal("Krap_max");
    buffer_rapidity = paraRdr_.getVal("buffer_rapidity");

    dKT   = (KT_max - KT_min)/(n_KT - 1);
    dKphi = 2*M_PI/n_Kphi;                  // does not need 0 and 2*pi

    for (int i = 0; i < n_KT; i++) {
        KT_array_.push_back(KT_min + i*dKT);
        number_of_pairs_numerator_KTdiff.push_back(0);
        number_of_pairs_denormenator_KTdiff.push_back(0);
    }

    for (int i = 0; i < n_Kphi; i++) {
        Kphi_array_.push_back(i*dKphi);
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
        create_a_5D_array(correl_3d_Kphi_diff_num, n_KT, n_Kphi,
                          qnpts, qnpts, qnpts);
        create_a_5D_array(correl_3d_Kphi_diff_num_count, n_KT, n_Kphi,
                          qnpts, qnpts, qnpts);
        create_a_5D_array(correl_3d_Kphi_diff_denorm, n_KT, n_Kphi,
                          qnpts, qnpts, qnpts);
        for (int iK = 0; iK < n_KT; iK++) {
            for (int iphi = 0; iphi < n_Kphi; iphi++) {
                for (int i = 0; i < qnpts; i++) {
                    for (int j = 0; j < qnpts; j++) {
                        for (int k = 0; k < qnpts; k++) {
                            q_out_diff_mean[iK][iphi][i][j][k] = 0.0;
                            q_side_diff_mean[iK][iphi][i][j][k] = 0.0;
                            q_long_diff_mean[iK][iphi][i][j][k] = 0.0;
                            correl_3d_Kphi_diff_num[iK][iphi][i][j][k] = 0.0;
                            correl_3d_Kphi_diff_num_count[iK][iphi][i][j][k] = 0.0;
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
        delete_a_5D_array(correl_3d_Kphi_diff_num_count,
                          n_KT, n_Kphi, qnpts, qnpts);
        delete_a_5D_array(correl_3d_Kphi_diff_denorm,
                          n_KT, n_Kphi, qnpts, qnpts);
    }

    if (azimuthal_flag_ == 1) {
        delete_a_2D_array(number_of_pairs_numerator_KTKphidiff, n_KT);
        delete_a_2D_array(number_of_pairs_denormenator_KTKphidiff, n_KT);
    }
}


void HBT_correlation::calculate_HBT_correlation_function(
                std::shared_ptr<particleSamples> particle_list_in) {
    set_particle_list(particle_list_in);
    int nev = particle_list->get_number_of_events();
    if (azimuthal_flag_ == 1) {
        calculate_flow_event_plane_angle(2);
    }

    // first pairs from the same event
    number_of_oversample_events_ = nev;
    messager.info("Compute pairs from the same event ...");
    messager.flush("info");

    std::vector<int> event_list(number_of_oversample_events_, 0);
    for (int isample = 0; isample < number_of_oversample_events_; isample++) {
        event_list[isample] = isample;
    }
    combine_and_bin_particle_pairs(event_list);

    // then pairs from mixed events
    int mixed_nev = particle_list->get_number_of_mixed_events();
    messager << "nev_mixed = " << mixed_nev;
    messager.flush("info");
    number_of_mixed_events_ = static_cast<int>(mixed_nev/2) + 1;
    messager.info("Compute pairs from the mixed event ...");
    for (int iev = 0; iev < nev; iev++) {
        messager << "progess: " << iev << "/" << nev;
        messager.flush("info");

        std::vector<int> mixed_event_list(number_of_mixed_events_, 0);
        int count = 0;
        while (count < number_of_mixed_events_) {
            int mixed_event_id = (
                    ran_gen_ptr_->rand_int_uniform() % mixed_nev);
            while (iev == mixed_event_id && mixed_nev != 1) {
                mixed_event_id = (
                    ran_gen_ptr_->rand_int_uniform() % mixed_nev);
            }
            mixed_event_list[count] = mixed_event_id;
            count++;
        }
        combine_and_bin_particle_pairs_mixed_events(iev, mixed_event_list);
    }
}

void HBT_correlation::output_HBTcorrelation() {
    if (invariant_radius_flag_ == 1)
        output_correlation_function_inv();

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
    int nev = particle_list->get_number_of_events();
    double vn_real = 0.0;
    double vn_imag = 0.0;
    for (int iev = 0; iev < nev; iev++) {
        int event_number_particle = (
                                particle_list->get_number_of_particles(iev));
        for (int i = 0; i < event_number_particle; i++) {
            double temp_px = particle_list->get_particle(iev, i).px;
            double temp_py = particle_list->get_particle(iev, i).py;
            double temp_phi = atan2(temp_py, temp_px);
            vn_real += cos(n_order*temp_phi);
            vn_imag += sin(n_order*temp_phi);
        }
    }
    psi_ref = atan2(vn_imag, vn_real)/n_order;
}


void HBT_correlation::combine_and_bin_particle_pairs(
                                                std::vector<int> event_list) {
    double hbarC_inv = 1./hbarC;
    std::vector<particle_info> temp_particle_list;
    double particle_list_rapidity_cut_max = tanh(Krap_max + buffer_rapidity);
    double particle_list_rapidity_cut_min = tanh(Krap_min - buffer_rapidity);
    for (int j = 0; j < number_of_oversample_events_; j++) {
        int event_id = event_list[j];
        int event_number_particle = (
                        particle_list->get_number_of_particles(event_list[j]));
        for (int i = 0; i < event_number_particle; i++) {
            double temp_E = particle_list->get_particle(event_id, i).E;
            double temp_pz = particle_list->get_particle(event_id, i).pz;
            double temp_ratio = temp_pz/temp_E;
            if (temp_ratio > particle_list_rapidity_cut_min 
                && temp_ratio < particle_list_rapidity_cut_max) {
                particle_info temp_particle;
                temp_particle.E = temp_E;
                temp_particle.pz = temp_pz;
                temp_particle.px = particle_list->get_particle(event_id, i).px;
                temp_particle.py = particle_list->get_particle(event_id, i).py;
                temp_particle.mass = (
                        particle_list->get_particle(event_id, i).mass);
                temp_particle.x = particle_list->get_particle(event_id, i).x;
                temp_particle.y = particle_list->get_particle(event_id, i).y;
                temp_particle.z = particle_list->get_particle(event_id, i).z;
                temp_particle.t = particle_list->get_particle(event_id, i).t;
                temp_particle_list.push_back(temp_particle);
            }
        }
    }
    long long int number_of_particles = temp_particle_list.size();
    unsigned long long int number_of_pairs = (
                    number_of_particles*(number_of_particles - 1)/2);

    // nested pair loop
    messager << "number of pairs: " << number_of_pairs;
    messager.flush("info");
    double rapidity_cut_max = tanh(Krap_max);
    double rapidity_cut_min = tanh(Krap_min);
    double KT_min_sq = KT_min*KT_min;
    double KT_max_sq = KT_max*KT_max;
    for (int i = 0; i < number_of_particles; i++) {
        double particle_1_px = temp_particle_list[i].px;
        double particle_1_py = temp_particle_list[i].py;
        double particle_1_pz = temp_particle_list[i].pz;
        double particle_1_E = temp_particle_list[i].E;
        double particle_1_x = temp_particle_list[i].x;
        double particle_1_y = temp_particle_list[i].y;
        double particle_1_z = temp_particle_list[i].z;
        double particle_1_t = temp_particle_list[i].t;

        for (int j = i+1; j < number_of_particles; j++) {
            double particle_2_px = temp_particle_list[j].px;
            double particle_2_py = temp_particle_list[j].py;
            double particle_2_pz = temp_particle_list[j].pz;
            double particle_2_E = temp_particle_list[j].E;
            double particle_2_x = temp_particle_list[j].x;
            double particle_2_y = temp_particle_list[j].y;
            double particle_2_z = temp_particle_list[j].z;
            double particle_2_t = temp_particle_list[j].t;

            double K_z = 0.5*(particle_1_pz + particle_2_pz);
            double K_E = 0.5*(particle_1_E + particle_2_E);
            double K_z_over_K_E = K_z/K_E;

            // rapidity cut
            if (K_z_over_K_E < rapidity_cut_min
                || K_z_over_K_E > rapidity_cut_max) {
                continue;
            }

            // K_T cut
            double K_x = 0.5*(particle_1_px + particle_2_px);
            double K_y = 0.5*(particle_1_py + particle_2_py);
            double K_perp_sq = K_x*K_x + K_y*K_y;
            if (K_perp_sq < KT_min_sq || K_perp_sq > KT_max_sq) {
                continue;
            }

            double K_perp = sqrt(K_perp_sq);
            int Kperp_idx = static_cast<int>((K_perp - KT_min)/dKT);

            // calculate q_inv
            double q_x = particle_1_px - particle_2_px;
            double q_y = particle_1_py - particle_2_py;
            double q_z = particle_1_pz - particle_2_pz;
            double q_E = particle_1_E - particle_2_E;
            double local_q_inv = sqrt(
                        - (q_E*q_E - q_x*q_x - q_y*q_y - q_z*q_z));

            double t_diff = particle_1_t - particle_2_t;
            double x_diff = particle_1_x - particle_2_x;
            double y_diff = particle_1_y - particle_2_y;
            double z_diff = particle_1_z - particle_2_z;

            if (invariant_radius_flag_ == 1) {
                if (local_q_inv > (q_min - delta_q/2. + 1e-8)
                    && local_q_inv > (q_max + delta_q/2. - 1e-8)) {
                    int qinv_idx = static_cast<int>(
                        (local_q_inv - (q_min - delta_q/2.))/delta_q);
                    if (qinv_idx < qnpts
                        && number_of_pairs_numerator_KTdiff[Kperp_idx]
                                    > needed_number_of_pairs) {
                        number_of_pairs_numerator_KTdiff[Kperp_idx]++;
                        double cos_qx = cos(hbarC_inv*(
                            q_E*t_diff - q_x*x_diff - q_y*y_diff - q_z*z_diff));
                        correl_1d_inv_num_count[Kperp_idx][qinv_idx]++;
                        q_inv_mean[Kperp_idx][qinv_idx] += local_q_inv;
                        correl_1d_inv_num[Kperp_idx][qinv_idx] += cos_qx;
                    }
                }
            }

            // calculate qout and qside in lcms
            double cos_K_phi = K_x/K_perp;
            double sin_K_phi = K_y/K_perp;

            double local_q_out = q_x*cos_K_phi + q_y*sin_K_phi;
            if (local_q_out < (q_min - delta_q/2. + 1e-8)
                || local_q_out > (q_max + delta_q/2. - 1e-8)) {
                continue;
            }

            int qout_idx = static_cast<int>(
                    (local_q_out - (q_min - delta_q/2.))/delta_q);
            if (qout_idx >= qnpts) continue;

            double local_q_side = q_y*cos_K_phi - q_x*sin_K_phi;
            if (local_q_side < (q_min - delta_q/2. + 1e-8)
                || local_q_side > (q_max + delta_q/2. - 1e-8)) {
                continue;
            }

            int qside_idx = static_cast<int>(
                    (local_q_side - (q_min - delta_q/2.))/delta_q);
            if (qside_idx >= qnpts) continue;

            // calcualte qlong in the lcms
            double Mt = sqrt(K_E*K_E - K_z*K_z);
            double boost_gamma = K_E/Mt;
            double boost_beta = K_z_over_K_E;
            // boost qz to lcms
            double local_q_long = boost_gamma*(q_z - boost_beta*q_E);

            if (local_q_long < (q_min - delta_q/2. + 1e-8)
                || local_q_long > (q_max + delta_q/2. - 1e-8)) {
                continue;
            }
            int qlong_idx = static_cast<int>(
                    (local_q_long - (q_min - delta_q/2.))/delta_q);
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
                    delta_phi += 2.*M_PI;
                }
                while (delta_phi > 2.*M_PI) {
                    delta_phi -= 2.*M_PI;
                }
                Kphi_idx = static_cast<int>(delta_phi/dKphi);
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

            double cos_qx = cos(hbarC_inv*(q_E*t_diff - q_x*x_diff
                                           - q_y*y_diff - q_z*z_diff));

            // bin results
            if (azimuthal_flag_ == 0) {
                correl_3d_num_count[Kperp_idx][qout_idx][qside_idx][qlong_idx]++;
                q_out_mean[Kperp_idx][qout_idx][qside_idx][qlong_idx] += local_q_out;
                q_side_mean[Kperp_idx][qout_idx][qside_idx][qlong_idx] += local_q_side;
                q_long_mean[Kperp_idx][qout_idx][qside_idx][qlong_idx] += local_q_long;
                correl_3d_num[Kperp_idx][qout_idx][qside_idx][qlong_idx] += cos_qx;
            } else {
                correl_3d_Kphi_diff_num_count[Kperp_idx][Kphi_idx][qout_idx][qside_idx][qlong_idx]++;
                q_out_diff_mean[Kperp_idx][Kphi_idx][qout_idx][qside_idx][qlong_idx] += local_q_out;
                q_side_diff_mean[Kperp_idx][Kphi_idx][qout_idx][qside_idx][qlong_idx] += local_q_side;
                q_long_diff_mean[Kperp_idx][Kphi_idx][qout_idx][qside_idx][qlong_idx] += local_q_long;
                correl_3d_Kphi_diff_num[Kperp_idx][Kphi_idx][qout_idx][qside_idx][qlong_idx] += cos_qx;
            }
        }
    }
    temp_particle_list.clear();
}

void HBT_correlation::combine_and_bin_particle_pairs_mixed_events(
                                        int event_id1,
                                        std::vector<int> mixed_event_list) {
    long long int number_of_particles_1 = (
            particle_list->get_number_of_particles(event_id1));
    double particle_list_rapidity_cut_max = tanh(Krap_max + buffer_rapidity);
    double particle_list_rapidity_cut_min = tanh(Krap_min - buffer_rapidity);
    std::vector<particle_info> temp_particle_list_1;
    for (int i = 0; i < number_of_particles_1; i++) {
        double temp_E = particle_list->get_particle(event_id1, i).E;
        double temp_pz = particle_list->get_particle(event_id1, i).pz;
        double temp_ratio = temp_pz/temp_E;
        if(temp_ratio > particle_list_rapidity_cut_min 
                && temp_ratio < particle_list_rapidity_cut_max) {
            particle_info temp_particle;
            temp_particle.E = temp_E;
            temp_particle.pz = temp_pz;
            temp_particle.px = particle_list->get_particle(event_id1, i).px;
            temp_particle.py = particle_list->get_particle(event_id1, i).py;
            temp_particle.mass = particle_list->get_particle(event_id1, i).mass;
            temp_particle.x = particle_list->get_particle(event_id1, i).x;
            temp_particle.y = particle_list->get_particle(event_id1, i).y;
            temp_particle.t = particle_list->get_particle(event_id1, i).t;
            temp_particle.z = particle_list->get_particle(event_id1, i).z;
            temp_particle_list_1.push_back(temp_particle);
        }
    }

    // prepare for the mixed events
    std::vector<particle_info> temp_particle_list_2;
    for (int iev = 0; iev < number_of_mixed_events_; iev++) {
        // introduce a random rotation for the mixed event
        const double random_rotation = ran_gen_ptr_->rand_uniform()*2*M_PI;
        const double cos_phi = cos(random_rotation);
        const double sin_phi = sin(random_rotation);
        int mixed_event_id = mixed_event_list[iev];
        int event_number_particle = (
            particle_list->get_number_of_particles_mixed_event(mixed_event_id));
        for (int i = 0; i < event_number_particle; i++) {
            double temp_E = (
                particle_list->get_particle_from_mixed_event(
                                                        mixed_event_id, i).E);
            double temp_pz = (
                particle_list->get_particle_from_mixed_event(
                                                        mixed_event_id, i).pz);
            double temp_ratio = temp_pz/temp_E;
            if(temp_ratio > particle_list_rapidity_cut_min 
               && temp_ratio < particle_list_rapidity_cut_max) {
                particle_info temp_particle;

                double px_temp = (
                    particle_list->get_particle_from_mixed_event(
                                                        mixed_event_id, i).px);
                double py_temp = (
                    particle_list->get_particle_from_mixed_event(
                                                        mixed_event_id, i).py);
                temp_particle.px = px_temp*cos_phi - py_temp*sin_phi;
                temp_particle.py = px_temp*sin_phi + py_temp*cos_phi;
                double x_temp = (
                        particle_list->get_particle_from_mixed_event(
                                                        mixed_event_id, i).x);
                double y_temp = (
                        particle_list->get_particle_from_mixed_event(
                                                        mixed_event_id, i).y);
                temp_particle.x = x_temp*cos_phi - y_temp*sin_phi;
                temp_particle.y = x_temp*sin_phi + y_temp*cos_phi;

                temp_particle.pz = temp_pz;
                temp_particle.E = temp_E;
                temp_particle.mass = (
                    particle_list->get_particle_from_mixed_event(
                                                    mixed_event_id, i).mass);
                temp_particle.z = (
                    particle_list->get_particle_from_mixed_event(
                                                    mixed_event_id, i).z);
                temp_particle.t = (
                    particle_list->get_particle_from_mixed_event(
                                                    mixed_event_id, i).t);
                temp_particle_list_2.push_back(temp_particle);
            }
        }
    }

    // nested pair loop
    number_of_particles_1 = temp_particle_list_1.size();
    long long int number_of_particles_2 = temp_particle_list_2.size();
    unsigned long long int number_of_pairs = (
            number_of_particles_1*number_of_particles_2);
    messager << "number of mixed pairs: " << number_of_pairs;
    messager.flush("info");
    double rapidity_cut_min = tanh(Krap_min);
    double rapidity_cut_max = tanh(Krap_max);
    double KT_min_sq = KT_min*KT_min;
    double KT_max_sq = KT_max*KT_max;
    for (int i = 0; i < number_of_particles_1; i++) {
        double particle_1_pz = temp_particle_list_1[i].pz;
        double particle_1_px = temp_particle_list_1[i].px;
        double particle_1_py = temp_particle_list_1[i].py;
        double particle_1_E = temp_particle_list_1[i].E;
        for (int j = 0; j < number_of_particles_2; j++) {
            double particle_2_pz = temp_particle_list_2[j].pz;
            double particle_2_px = temp_particle_list_2[j].px;
            double particle_2_py = temp_particle_list_2[j].py;
            double particle_2_E = temp_particle_list_2[j].E;

            double K_z = 0.5*(particle_1_pz + particle_2_pz);
            double K_E = 0.5*(particle_1_E + particle_2_E);
            double K_z_over_K_E = K_z/K_E;
            if (K_z_over_K_E < rapidity_cut_min 
                || K_z_over_K_E > rapidity_cut_max) {
                continue;
            }

            double K_x = 0.5*(particle_1_px + particle_2_px);
            double K_y = 0.5*(particle_1_py + particle_2_py);
            double K_perp_sq = K_x*K_x + K_y*K_y;
            if (K_perp_sq < KT_min_sq || K_perp_sq > KT_max_sq) {
                continue;
            }

            double K_perp = sqrt(K_perp_sq);
            int Kperp_idx = static_cast<int>((K_perp - KT_min)/dKT);

            double q_z = particle_1_pz - particle_2_pz;
            double q_E = particle_1_E - particle_2_E;
            double q_x = particle_1_px - particle_2_px;
            double q_y = particle_1_py - particle_2_py;
            double local_q_inv = sqrt(
                        - (q_E*q_E - q_x*q_x - q_y*q_y - q_z*q_z));

            if (invariant_radius_flag_ == 1) {
                if (local_q_inv > (q_min - delta_q/2. + 1e-8)
                    && local_q_inv < (q_max + delta_q/2. - 1e-8)) {
                    int qinv_idx = static_cast<int>(
                        (local_q_inv - (q_min - delta_q/2.))/delta_q);
                    if (qinv_idx < qnpts
                        && number_of_pairs_numerator_KTdiff[Kperp_idx]
                                > needed_number_of_pairs) {
                        number_of_pairs_denormenator_KTdiff[Kperp_idx]++;
                        correl_1d_inv_denorm[Kperp_idx][qinv_idx] += 1.0;
                    }
                }
            }

            double cos_K_phi = K_x/K_perp;
            double sin_K_phi = K_y/K_perp;
            double local_q_out = q_x*cos_K_phi + q_y*sin_K_phi;
            if (local_q_out < (q_min - delta_q/2. + 1e-8)
                || local_q_out >= (q_max + delta_q/2. - 1e-8)) {
                continue;
            }
            int qout_idx = static_cast<int>(
                    (local_q_out - (q_min - delta_q/2.))/delta_q);
            if (qout_idx >= qnpts) continue;

            double local_q_side = q_y*cos_K_phi - q_x*sin_K_phi;
            if (local_q_side < (q_min - delta_q/2. + 1e-8)
                || local_q_side >= (q_max + delta_q/2. - 1e-8)) {
                continue;
            }

            int qside_idx = static_cast<int>(
                    (local_q_side - (q_min - delta_q/2.))/delta_q);
            if (qside_idx >= qnpts) continue;

            // calcualte qlong in the lcms
            double Mt = sqrt(K_E*K_E - K_z*K_z);
            double boost_gamma = K_E/Mt;
            double boost_beta = K_z_over_K_E;
            // boost qz to lcms
            double local_q_long = boost_gamma*(q_z - boost_beta*q_E);

            if (local_q_long < (q_min - delta_q/2. + 1e-8)
                || local_q_long >= (q_max + delta_q/2. - 1e-8)) {
                continue;
            }

            int qlong_idx = static_cast<int>(
                    (local_q_long - (q_min - delta_q/2.))/delta_q);
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
                    delta_phi += 2.*M_PI;
                }
                while (delta_phi > 2.*M_PI) {
                    delta_phi -= 2.*M_PI;
                }
                Kphi_idx = static_cast<int>(delta_phi/dKphi);
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
                correl_3d_denorm[Kperp_idx][qout_idx][qside_idx][qlong_idx] += 1.0;
            } else {
                correl_3d_Kphi_diff_denorm[Kperp_idx][Kphi_idx][qout_idx][qside_idx][qlong_idx] += 1.0;
            }
        }
    }
    temp_particle_list_1.clear();
    temp_particle_list_2.clear();
}


void HBT_correlation::output_correlation_function_inv() {
    for (int iK = 0; iK < n_KT - 1; iK++) {
        double npair_ratio = (
                static_cast<double>(number_of_pairs_numerator_KTdiff[iK])
                /static_cast<double>(number_of_pairs_denormenator_KTdiff[iK]));
        std::ostringstream filename;
        filename << path_ << "/HBT_correlation_function_inv_KT_"
                 << KT_array_[iK] << "_" << KT_array_[iK+1] << ".dat";
        std::ofstream output(filename.str().c_str());
        for (int iqinv = 0; iqinv < qnpts; iqinv++) {
            int npart_num = correl_1d_inv_num_count[iK][iqinv];
            int npart_denorm = correl_1d_inv_denorm[iK][iqinv];
            double q_inv_local;
            double correl_fun_num, correl_fun_denorm, correl_fun_val;
            if (npart_num < 2 || npart_denorm < 2) {
                q_inv_local = q_out[iqinv];
                correl_fun_num = 0.0;
                correl_fun_denorm = npart_denorm;
                correl_fun_val = 0.0;
            } else {
                q_inv_local = q_inv_mean[iK][iqinv]/npart_num;
                correl_fun_num = correl_1d_inv_num[iK][iqinv];
                correl_fun_denorm = correl_1d_inv_denorm[iK][iqinv];
                correl_fun_val = (correl_fun_num
                                  /(correl_fun_denorm*npair_ratio));
            }

            if (q_out[iqinv] > 0.) {
                output << std::scientific << std::setw(18)
                       << std::setprecision(8) 
                       << q_inv_local << "    "
                       << npart_num << "    " << correl_fun_num << "    "
                       << correl_fun_denorm << "    "
                       << correl_fun_val << "    " << 0.0 << std::endl;
            }
        }
        output.close();
    }
}

void HBT_correlation::output_correlation_function() {
    for (int iK = 0; iK < n_KT - 1; iK++) {
        double npair_ratio = (
                static_cast<double>(number_of_pairs_numerator_KTdiff[iK])
                /static_cast<double>(number_of_pairs_denormenator_KTdiff[iK]));
        std::ostringstream filename;
        filename << path_ << "/HBT_correlation_function_KT_"
                 << KT_array_[iK] << "_" << KT_array_[iK+1] << ".dat";
        std::ofstream output(filename.str().c_str());
        for (int iqlong = 0; iqlong < qnpts; iqlong++) {
            for (int iqout = 0; iqout < qnpts; iqout++) {
                for (int iqside = 0; iqside < qnpts; iqside++) {
                    int npart_num = 
                        correl_3d_num_count[iK][iqout][iqside][iqlong];
                    int npart_denorm = 
                        correl_3d_denorm[iK][iqout][iqside][iqlong];
                    double q_out_local, q_side_local, q_long_local;
                    double correl_fun_num, correl_fun_denorm, correl_fun_val;
                    if (npart_num < 2 || npart_denorm < 2) {
                        q_out_local = q_out[iqout];
                        q_side_local = q_side[iqside];
                        q_long_local = q_long[iqlong];
                        correl_fun_num = 0.0;
                        correl_fun_denorm = npart_denorm;
                        correl_fun_val = 0.0;
                    } else {
                        q_out_local = (q_out_mean[iK][iqout][iqside][iqlong]
                                       /npart_num);
                        q_side_local = (q_side_mean[iK][iqout][iqside][iqlong]
                                        /npart_num);
                        q_long_local = (
                            q_long_mean[iK][iqout][iqside][iqlong]/npart_num);

                        correl_fun_num = (
                                correl_3d_num[iK][iqout][iqside][iqlong]);
                        correl_fun_denorm = (
                                correl_3d_denorm[iK][iqout][iqside][iqlong]);
                        correl_fun_val = (correl_fun_num
                                          /(correl_fun_denorm*npair_ratio));
                    }

                    output << std::scientific << std::setw(18)
                           << std::setprecision(8) 
                           << q_out_local << "    " << q_side_local << "    "
                           << q_long_local << "    "
                           << npart_num << "    " << correl_fun_num << "    "
                           << correl_fun_denorm << "    "
                           << correl_fun_val << "    " << 0.0 << std::endl;
                }
            }
        }
        output.close();
    }
}

void HBT_correlation::output_correlation_function_Kphi_differential() {
    for (int iK = 0; iK < n_KT - 1; iK++) {
        for (int iKphi = 0; iKphi < n_Kphi; iKphi++) {
            double npair_ratio = (
                static_cast<double>(number_of_pairs_numerator_KTKphidiff[iK][iKphi])
                /static_cast<double>(number_of_pairs_denormenator_KTKphidiff[iK][iKphi]));
            std::ostringstream filename;
            filename << path_ << "/HBT_correlation_function_KT_"
                     << KT_array_[iK] << "_" << KT_array_[iK+1] << "_Kphi_"
                     << Kphi_array_[iKphi] << ".dat";
            std::ofstream output(filename.str().c_str());
            for (int iqlong = 0; iqlong < qnpts; iqlong++) {
                for (int iqout = 0; iqout < qnpts; iqout++) {
                    for (int iqside = 0; iqside < qnpts; iqside++) {
                        int npart_num = 
                            correl_3d_Kphi_diff_num_count[iK][iKphi][iqout][iqside][iqlong];
                        int npart_denorm = 
                            correl_3d_Kphi_diff_denorm[iK][iKphi][iqout][iqside][iqlong];
                        double q_out_local, q_side_local, q_long_local;
                        double correl_fun_num, correl_fun_denorm;
                        double correl_fun_val;
                        if (npart_num < 2 || npart_denorm < 2) {
                            q_out_local = q_out[iqout];
                            q_side_local = q_side[iqside];
                            q_long_local = q_long[iqlong];
                            correl_fun_num = 0.0;
                            correl_fun_denorm = npart_denorm;
                            correl_fun_val = 0.0;
                        } else {
                            q_out_local = (
                                q_out_diff_mean[iK][iKphi][iqout][iqside][iqlong]
                                /npart_num);
                            q_side_local = (
                                q_side_diff_mean[iK][iKphi][iqout][iqside][iqlong]
                                /npart_num);
                            q_long_local = (
                                q_long_diff_mean[iK][iKphi][iqout][iqside][iqlong]
                                /npart_num);

                            correl_fun_num = correl_3d_Kphi_diff_num[iK][iKphi][iqout][iqside][iqlong];
                            correl_fun_denorm = correl_3d_Kphi_diff_denorm[iK][iKphi][iqout][iqside][iqlong];
                            correl_fun_val = (correl_fun_num
                                              /correl_fun_denorm/npair_ratio);
                        }

                        output << std::scientific << std::setw(18)
                               << std::setprecision(8)
                               << q_out_local << "    "
                               << q_side_local << "    "
                               << q_long_local << "    "
                               << npart_num << "    "
                               << correl_fun_num << "    "
                               << correl_fun_denorm << "    "
                               << correl_fun_val << "    " << 0.0 << std::endl;
                    }
                }
            }
            output.close();
        }
    }
}


template <typename T>
void HBT_correlation::create_a_2D_array(T **&arr2D, int nx, int ny) {
    arr2D = new T* [nx];
    for (int i = 0; i < nx; i++) {
        arr2D[i] = new T [ny];
    }
}

template <typename T>
void HBT_correlation::delete_a_2D_array(T **&arr2D, int nx) {
    for (int i = 0; i < nx; i++)
        delete [] arr2D[i];
    delete [] arr2D;
}

template <typename T>
void HBT_correlation::create_a_3D_array(T ***&arr3D, int nx, int ny, int nz) {
    arr3D = new T** [nx];
    for (int i = 0; i < nx; i++) {
        create_a_2D_array(arr3D[i], ny, nz);
    }
}

template <typename T>
void HBT_correlation::delete_a_3D_array(T ***&arr3D, int nx, int ny) {
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < nx; j++)
            delete[] arr3D[i][j];
        delete[] arr3D[i];
    }
    delete[] arr3D;
}

template <typename T>
void HBT_correlation::create_a_4D_array(
                T ****&arr4D, int n1, int n2, int n3, int n4) {
    arr4D = new T*** [n1];
    for (int i = 0; i < n1; i++) {
        create_a_3D_array(arr4D[i], n2, n3, n4);
    }
}

template <typename T>
void HBT_correlation::delete_a_4D_array(T ****&arr4D, int n1, int n2, int n3) {
    for (int i = 0; i < n1; i++) {
        for (int j = 0; j < n2; j++) {
            for (int k = 0; k < n3; k++)
                delete[] arr4D[i][j][k];
            delete[] arr4D[i][j];
        }
        delete[] arr4D[i];
    }
    delete[] arr4D;
}

template <typename T>
void HBT_correlation::create_a_5D_array(
                T *****&arr5D, int n1, int n2, int n3, int n4, int n5) {
    arr5D = new T**** [n1];
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
                for (int l = 0; l < n4; l++)
                    delete[] arr5D[i][j][k][l];
                delete[] arr5D[i][j][k];
            }
            delete[] arr5D[i][j];
        }
        delete[] arr5D[i];
    }
    delete[] arr5D;
}
