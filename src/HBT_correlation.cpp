#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <iomanip>
#include "parameters.h"
#include "HBT_correlation.h"
using namespace std;

HBT_correlation::HBT_correlation(ParameterReader* paraRdr_in, string path_in, 
                                 particleSamples *particle_list_in) {
    paraRdr = paraRdr_in;
    path = path_in;
    particle_list = particle_list_in;

    qnpts = paraRdr->getVal("qnpts");
    q_min = paraRdr->getVal("q_min");
    q_max = paraRdr->getVal("q_max");
    delta_q = (q_max - q_min)/(qnpts - 1);
    
    q_out = new double [qnpts];
    q_side = new double [qnpts];
    q_long = new double [qnpts];
    for (int i = 0; i < qnpts; i++) {
        q_out[i] = q_min + i*delta_q;
        q_side[i] = q_min + i*delta_q;
        q_long[i] = q_min + i*delta_q;
    }

    needed_number_of_pairs = paraRdr->getVal("needed_number_of_pairs");

    azimuthal_flag = paraRdr->getVal("azimuthal_flag");
    n_KT = paraRdr->getVal("n_KT");
    n_Kphi = paraRdr->getVal("n_Kphi");
    KT_min = paraRdr->getVal("KT_min");
    KT_max = paraRdr->getVal("KT_max");
    Krap_min = paraRdr->getVal("Krap_min");
    Krap_max = paraRdr->getVal("Krap_max");
    buffer_rapidity = paraRdr->getVal("buffer_rapidity");

    dKT = (KT_max - KT_min)/(n_KT - 1);
    dKphi = 2*M_PI/(n_Kphi - 1);
    
    KT_array = new double [n_KT];
    number_of_pairs_numerator_KTdiff = new unsigned long long int [n_KT];
    number_of_pairs_denormenator_KTdiff = new unsigned long long int [n_KT];
    for (int i = 0; i < n_KT; i++) {
        KT_array[i] = KT_min + i*dKT;
        number_of_pairs_numerator_KTdiff[i] = 0;
        number_of_pairs_denormenator_KTdiff[i] = 0;
    }
      
    Kphi_array = new double [n_Kphi];
    for (int i = 0; i < n_Kphi; i++)
        Kphi_array[i] = 0.0 + i*dKphi;

    number_of_pairs_numerator_KTKphidiff = new unsigned long long int* [n_KT];
    number_of_pairs_denormenator_KTKphidiff = (
                                          new unsigned long long int* [n_KT]);
    for (int i = 0; i < n_KT; i++) {
        number_of_pairs_numerator_KTKphidiff[i] = (
                                          new unsigned long long int [n_Kphi]);
        number_of_pairs_denormenator_KTKphidiff[i] = (
                                          new unsigned long long int [n_Kphi]);
        for (int j = 0; j < n_Kphi; j++) {
            number_of_pairs_numerator_KTKphidiff[i] = 0;
            number_of_pairs_denormenator_KTKphidiff[i] = 0;
        }
    }
    
    number_of_mixed_events = paraRdr->getVal("number_of_mixed_events");
    number_of_oversample_events = 
                             paraRdr->getVal("number_of_oversample_events");
    number_pairs_num = 0;
    number_pairs_denorm = 0;
    if (azimuthal_flag == 0) {
        q_out_mean = new double ***[n_KT];
        q_side_mean = new double ***[n_KT];
        q_long_mean = new double ***[n_KT];
        correl_3d_num = new double *** [n_KT];
        correl_3d_num_count = new double *** [n_KT];
        correl_3d_denorm = new double *** [n_KT];
        for (int iK = 0; iK < n_KT; iK++) {
            q_out_mean[iK] = new double ** [qnpts];
            q_side_mean[iK] = new double ** [qnpts];
            q_long_mean[iK] = new double ** [qnpts];
            correl_3d_num[iK] = new double ** [qnpts];
            correl_3d_num_count[iK] = new double ** [qnpts];
            correl_3d_denorm[iK] = new double ** [qnpts];
            for (int i = 0; i < qnpts; i++) {
                q_out_mean[iK][i] = new double * [qnpts];
                q_side_mean[iK][i] = new double * [qnpts];
                q_long_mean[iK][i] = new double * [qnpts];
                correl_3d_num[iK][i] = new double * [qnpts];
                correl_3d_num_count[iK][i] = new double * [qnpts];
                correl_3d_denorm[iK][i] = new double * [qnpts];
                for (int j = 0; j < qnpts; j++) {
                    q_out_mean[iK][i][j] = new double [qnpts];
                    q_side_mean[iK][i][j] = new double [qnpts];
                    q_long_mean[iK][i][j] = new double [qnpts];
                    correl_3d_num[iK][i][j] = new double [qnpts];
                    correl_3d_num_count[iK][i][j] = new double [qnpts];
                    correl_3d_denorm[iK][i][j] = new double [qnpts];
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
        q_out_diff_mean = new double **** [n_KT];
        q_side_diff_mean = new double **** [n_KT];
        q_long_diff_mean = new double **** [n_KT];
        correl_3d_Kphi_diff_num = new double **** [n_KT];
        correl_3d_Kphi_diff_num_count = new double **** [n_KT];
        correl_3d_Kphi_diff_denorm = new double **** [n_KT];
        for (int iK = 0; iK < n_KT; iK++) {
            q_out_diff_mean[iK] = new double *** [n_Kphi];
            q_side_diff_mean[iK] = new double *** [n_Kphi];
            q_long_diff_mean[iK] = new double *** [n_Kphi];
            correl_3d_Kphi_diff_num[iK] = new double *** [n_Kphi];
            correl_3d_Kphi_diff_num_count[iK] = new double *** [n_Kphi];
            correl_3d_Kphi_diff_denorm[iK] = new double *** [n_Kphi];
            for (int iphi = 0; iphi < n_Kphi; iphi++) {
                q_out_diff_mean[iK][iphi] = new double ** [qnpts];
                q_side_diff_mean[iK][iphi] = new double ** [qnpts];
                q_long_diff_mean[iK][iphi] = new double ** [qnpts];
                correl_3d_Kphi_diff_num[iK][iphi] = new double ** [qnpts];
                correl_3d_Kphi_diff_num_count[iK][iphi] = new double ** [qnpts];
                correl_3d_Kphi_diff_denorm[iK][iphi] = new double ** [qnpts];
                for (int i = 0; i < qnpts; i++) {
                    q_out_diff_mean[iK][iphi][i] = new double * [qnpts];
                    q_side_diff_mean[iK][iphi][i] = new double * [qnpts];
                    q_long_diff_mean[iK][iphi][i] = new double * [qnpts];
                    correl_3d_Kphi_diff_num[iK][iphi][i] = (
                                                        new double * [qnpts]);
                    correl_3d_Kphi_diff_num_count[iK][iphi][i] = (
                                                        new double * [qnpts]);
                    correl_3d_Kphi_diff_denorm[iK][iphi][i] = (
                                                        new double * [qnpts]);
                    for (int j = 0; j < qnpts; j++) {
                        q_out_diff_mean[iK][iphi][i][j] = new double [qnpts];
                        q_side_diff_mean[iK][iphi][i][j] = new double [qnpts];
                        q_long_diff_mean[iK][iphi][i][j] = new double [qnpts];
                        correl_3d_Kphi_diff_num[iK][iphi][i][j] = (
                                                        new double [qnpts]);
                        correl_3d_Kphi_diff_num_count[iK][iphi][i][j] = (
                                                        new double [qnpts]);
                        correl_3d_Kphi_diff_denorm[iK][iphi][i][j] = (
                                                        new double [qnpts]);
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
    delete[] q_out;
    delete[] q_side;
    delete[] q_long;
    if (azimuthal_flag == 0) {
        for (int iK = 0; iK < n_KT; iK++) {
            for (int i = 0; i < qnpts; i++) {
                for (int j = 0; j < qnpts; j++) {
                    delete[] q_out_mean[iK][i][j];
                    delete[] q_side_mean[iK][i][j];
                    delete[] q_long_mean[iK][i][j];
                    delete[] correl_3d_num[iK][i][j];
                    delete[] correl_3d_num_count[iK][i][j];
                    delete[] correl_3d_denorm[iK][i][j];
                }
                delete[] q_out_mean[iK][i];
                delete[] q_side_mean[iK][i];
                delete[] q_long_mean[iK][i];
                delete[] correl_3d_num[iK][i];
                delete[] correl_3d_num_count[iK][i];
                delete[] correl_3d_denorm[iK][i];
            }
            delete[] q_out_mean[iK];
            delete[] q_side_mean[iK];
            delete[] q_long_mean[iK];
            delete[] correl_3d_num[iK];
            delete[] correl_3d_num_count[iK];
            delete[] correl_3d_denorm[iK];
        }
        delete[] q_out_mean;
        delete[] q_side_mean;
        delete[] q_long_mean;
        delete[] correl_3d_num;
        delete[] correl_3d_num_count;
        delete[] correl_3d_denorm;
    } else {
        for (int iK = 0; iK < n_KT; iK++) {
            for (int iphi = 0; iphi < n_Kphi; iphi++) {
                for (int i = 0; i < qnpts; i++) {
                    for (int j = 0; j < qnpts; j++) {
                        delete[] q_out_diff_mean[iK][iphi][i][j];
                        delete[] q_side_diff_mean[iK][iphi][i][j];
                        delete[] q_long_diff_mean[iK][iphi][i][j];
                        delete[] correl_3d_Kphi_diff_num[iK][iphi][i][j];
                        delete[] correl_3d_Kphi_diff_num_count[iK][iphi][i][j];
                        delete[] correl_3d_Kphi_diff_denorm[iK][iphi][i][j];
                    }
                    delete[] q_out_diff_mean[iK][iphi][i];
                    delete[] q_side_diff_mean[iK][iphi][i];
                    delete[] q_long_diff_mean[iK][iphi][i];
                    delete[] correl_3d_Kphi_diff_num[iK][iphi][i];
                    delete[] correl_3d_Kphi_diff_num_count[iK][iphi][i];
                    delete[] correl_3d_Kphi_diff_denorm[iK][iphi][i];
                }
                delete[] q_out_diff_mean[iK][iphi];
                delete[] q_side_diff_mean[iK][iphi];
                delete[] q_long_diff_mean[iK][iphi];
                delete[] correl_3d_Kphi_diff_num[iK][iphi];
                delete[] correl_3d_Kphi_diff_num_count[iK][iphi];
                delete[] correl_3d_Kphi_diff_denorm[iK][iphi];
            }
            delete[] q_out_diff_mean[iK];
            delete[] q_side_diff_mean[iK];
            delete[] q_long_diff_mean[iK];
            delete[] correl_3d_Kphi_diff_num[iK];
            delete[] correl_3d_Kphi_diff_num_count[iK];
            delete[] correl_3d_Kphi_diff_denorm[iK];
        }
        delete[] q_out_diff_mean;
        delete[] q_side_diff_mean;
        delete[] q_long_diff_mean;
        delete[] correl_3d_Kphi_diff_num;
        delete[] correl_3d_Kphi_diff_num_count;
        delete[] correl_3d_Kphi_diff_denorm;
    }
    delete [] KT_array;
    for (int i = 0; i < n_KT; i++) {
        delete[] number_of_pairs_numerator_KTKphidiff[i];
        delete[] number_of_pairs_denormenator_KTKphidiff[i];
    }
    delete[] number_of_pairs_numerator_KTKphidiff;
    delete[] number_of_pairs_denormenator_KTKphidiff;
    delete[] number_of_pairs_numerator_KTdiff;
    delete[] number_of_pairs_denormenator_KTdiff;
}

void HBT_correlation::calculate_HBT_correlation_function() {
    int event_id = 0;
    int buffer_size = particle_list->get_event_buffer_size();
    while (!particle_list->end_of_file()) {
        cout << "Reading event: " << event_id+1 << "-" 
             << event_id + buffer_size << " ... " << flush;
        particle_list->read_in_particle_samples();
        particle_list->read_in_particle_samples_mixed_event();
        cout << " processing ..." << endl;
        int nev = particle_list->get_number_of_events();
        int n_skip_ev = (int)(nev/number_of_oversample_events);
        // first pairs from the same event
        cout << "Compute pairs from the same event ..." << endl;
        for (int iev = 0; iev < n_skip_ev; iev++) {
            cout << "progess: " << iev << "/" << n_skip_ev << endl;
            int *event_list = new int [number_of_oversample_events];
            for (int isample = 0; isample < number_of_oversample_events;
                    isample++) {
                event_list[isample] = iev + isample*n_skip_ev;
            }
            combine_and_bin_particle_pairs(event_list);
            delete[] event_list;
        }

        // then pairs from mixed events
        cout << "Compute pairs from the mixed event ..." << endl;
        for (int iev = 0; iev < nev; iev++) {
            cout << "progess: " << iev << "/" << nev << endl;
            int *mixed_event_list = new int [number_of_mixed_events];
            int count = 0;
            while (count < number_of_mixed_events) {
                int mixed_event_id = rand() % nev;
                mixed_event_list[count] = mixed_event_id;
                count++;
            }
            combine_and_bin_particle_pairs_mixed_events(iev, mixed_event_list);
            event_id++;
            delete [] mixed_event_list;
        }
    }
    if (azimuthal_flag == 0)
        output_correlation_function();
    else
        output_correlation_function_Kphi_differential();
}

void HBT_correlation::combine_and_bin_particle_pairs(int* event_list) {
    double hbarC_inv = 1./hbarC;
    vector<particle_info> temp_particle_list;
    double particle_list_rapidity_cut_max = tanh(Krap_max + buffer_rapidity);
    double particle_list_rapidity_cut_min = tanh(Krap_min - buffer_rapidity);
    for (int j = 0; j < number_of_oversample_events; j++) {
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
    cout << "number of pairs: " << number_of_pairs << endl;
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
            if (K_z_over_K_E > rapidity_cut_min 
                && K_z_over_K_E < rapidity_cut_max) {
                // rapidity cut
                double K_x = 0.5*(particle_1_px + particle_2_px);
                double K_y = 0.5*(particle_1_py + particle_2_py);
                double K_perp_sq = K_x*K_x + K_y*K_y;
                if (K_perp_sq > KT_min_sq && K_perp_sq < KT_max_sq) {
                    // K_T cut
                    double K_perp = sqrt(K_perp_sq);
                    int Kperp_idx = (int)((K_perp - KT_min)/dKT);

                    // calculate qout and qside in lcms
                    double cos_K_phi = K_x/K_perp;
                    double sin_K_phi = K_y/K_perp;
                    double q_x = particle_1_px - particle_2_px;
                    double q_y = particle_1_py - particle_2_py;

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
                    double q_z = particle_1_pz - particle_2_pz;
                    double q_E = particle_1_E - particle_2_E;
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

                    double local_K_phi;
                    int Kphi_idx;
                    if (azimuthal_flag == 0) {
                        if (number_of_pairs_numerator_KTdiff[Kperp_idx]
                            > needed_number_of_pairs) {
                            continue;
                        }
                        number_of_pairs_numerator_KTdiff[Kperp_idx]++;
                    } else {
                        local_K_phi = atan2(K_y, K_x);
                        Kphi_idx = (int)((local_K_phi)/dKphi);
                        if (number_of_pairs_numerator_KTKphidiff[Kperp_idx][Kphi_idx]
                                > needed_number_of_pairs) {
                            continue;
                        }
                        number_of_pairs_numerator_KTKphidiff[Kperp_idx][Kphi_idx]++;
                    }
                    
                    double t_diff = particle_1_t - particle_2_t;
                    double x_diff = particle_1_x - particle_2_x;
                    double y_diff = particle_1_y - particle_2_y;
                    double z_diff = particle_1_z - particle_2_z;

                    double cos_qx = cos(hbarC_inv
                        *(q_E*t_diff - q_x*x_diff - q_y*y_diff - q_z*z_diff));

                    // bin results
                    if (azimuthal_flag == 0) {
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
        }
    }
    temp_particle_list.clear();
}

void HBT_correlation::combine_and_bin_particle_pairs_mixed_events(
                                        int event_id1, int* mixed_event_list) {
    long long int number_of_particles_1 = (
            particle_list->get_number_of_particles(event_id1));
    double particle_list_rapidity_cut_max = tanh(Krap_max + buffer_rapidity);
    double particle_list_rapidity_cut_min = tanh(Krap_min - buffer_rapidity);
    vector<particle_info> temp_particle_list_1;
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
    vector<particle_info> temp_particle_list_2;
    for (int iev = 0; iev < number_of_mixed_events; iev++) {
        // introduce a random rotation for the mixed event
        double random_rotation = drand48()*2*M_PI;
        double cos_phi = cos(random_rotation);
        double sin_phi = sin(random_rotation);
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
    cout << "number of mixed pairs: " << number_of_pairs << endl;
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
            if (K_z_over_K_E > rapidity_cut_min 
                && K_z_over_K_E < rapidity_cut_max) {
                double K_x = 0.5*(particle_1_px + particle_2_px);
                double K_y = 0.5*(particle_1_py + particle_2_py);
                double K_perp_sq = K_x*K_x + K_y*K_y;
                if (K_perp_sq > KT_min_sq && K_perp_sq < KT_max_sq) {
                    // calcualte qout, qside first
                    double K_perp = sqrt(K_perp_sq);
                    int Kperp_idx = (int)((K_perp - KT_min)/dKT);
                    double cos_K_phi = K_x/K_perp;
                    double sin_K_phi = K_y/K_perp;
                    double q_x = particle_1_px - particle_2_px;
                    double q_y = particle_1_py - particle_2_py;

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
                    double q_z = particle_1_pz - particle_2_pz;
                    double q_E = particle_1_E - particle_2_E;
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
                    
                    double local_K_phi;
                    int Kphi_idx;
                    if (azimuthal_flag == 0) {
                        if (number_of_pairs_denormenator_KTdiff[Kperp_idx]
                                        > needed_number_of_pairs) {
                            continue;
                        }
                        number_of_pairs_denormenator_KTdiff[Kperp_idx]++;
                    } else {
                        local_K_phi = atan2(K_y, K_x);
                        Kphi_idx = (int)((local_K_phi)/dKphi);
                        if (number_of_pairs_denormenator_KTKphidiff[Kperp_idx][Kphi_idx]
                                        > needed_number_of_pairs) {
                            continue;
                        }
                        number_of_pairs_denormenator_KTKphidiff[Kperp_idx][Kphi_idx]++;
                    }
                    
                    // bin results
                    if (azimuthal_flag == 0) {
                        correl_3d_denorm[Kperp_idx][qout_idx][qside_idx][qlong_idx] += 1.0;
                    } else {
                        correl_3d_Kphi_diff_denorm[Kperp_idx][Kphi_idx][qout_idx][qside_idx][qlong_idx] += 1.0;
                    }
                }
            }
        }
    }
    temp_particle_list_1.clear();
    temp_particle_list_2.clear();
}

void HBT_correlation::output_correlation_function()
{
    for(int iK = 0; iK < n_KT - 1; iK++)
    {
        double npair_ratio = ((double)number_of_pairs_numerator_KTdiff[iK]
                            /(double)number_of_pairs_denormenator_KTdiff[iK]);
        ostringstream filename;
        filename << path << "/HBT_correlation_function_KT_" 
                 << KT_array[iK] << "_" << KT_array[iK+1] << ".dat";
        ofstream output(filename.str().c_str());
        for(int iqlong = 0; iqlong < qnpts; iqlong++)
        {
            for(int iqout = 0; iqout < qnpts; iqout++)
            {
                for(int iqside = 0; iqside < qnpts; iqside++)
                {
                    int npart_num = 
                        correl_3d_num_count[iK][iqout][iqside][iqlong];
                    int npart_denorm = 
                        correl_3d_denorm[iK][iqout][iqside][iqlong];
                    double q_out_local, q_side_local, q_long_local;
                    double correl_fun_num, correl_fun_denorm, correl_fun_val;
                    if(npart_num < 2 || npart_denorm < 2)
                    {
                        q_out_local = q_out[iqout];
                        q_side_local = q_side[iqside];
                        q_long_local = q_long[iqlong];
                        correl_fun_num = 0.0;
                        correl_fun_denorm = npart_denorm;
                        correl_fun_val = 0.0;
                    }
                    else
                    {
                        q_out_local = (q_out_mean[iK][iqout][iqside][iqlong]
                                       /npart_num);
                        q_side_local = (q_side_mean[iK][iqout][iqside][iqlong]
                                        /npart_num);
                        q_long_local = (
                            q_long_mean[iK][iqout][iqside][iqlong]/npart_num);

                        correl_fun_num = correl_3d_num[iK][iqout][iqside][iqlong];
                        correl_fun_denorm = correl_3d_denorm[iK][iqout][iqside][iqlong];
                        correl_fun_val = (correl_fun_num
                                          /(correl_fun_denorm*npair_ratio));
                    }

                    output << scientific << setw(18) << setprecision(8) 
                           << q_out_local << "    " << q_side_local << "    " 
                           << q_long_local << "    "
                           << npart_num << "    " << correl_fun_num << "    "  
                           << correl_fun_denorm << "    "
                           << correl_fun_val << "    " << 0.0 << endl;
                }
            }
        }
        output.close();
    }
}

void HBT_correlation::output_correlation_function_Kphi_differential()
{
    for(int iK = 0; iK < n_KT - 1; iK++)
    {
        for(int iKphi = 0; iKphi < n_Kphi; iKphi++)
        {
            double npair_ratio = (
                (double)number_of_pairs_numerator_KTKphidiff[iK][iKphi]
                /(double)number_of_pairs_denormenator_KTKphidiff[iK][iKphi]);
            ostringstream filename;
            filename << path << "/HBT_correlation_function_KT_" 
                     << KT_array[iK] << "_" << KT_array[iK+1] << "_Kphi_" 
                     << Kphi_array[iKphi] << ".dat";
            ofstream output(filename.str().c_str());
            for(int iqlong = 0; iqlong < qnpts; iqlong++)
            {
                for(int iqout = 0; iqout < qnpts; iqout++)
                {
                    for(int iqside = 0; iqside < qnpts; iqside++)
                    {
                        int npart_num = 
                            correl_3d_Kphi_diff_num_count[iK][iKphi][iqout][iqside][iqlong];
                        int npart_denorm = 
                            correl_3d_Kphi_diff_denorm[iK][iKphi][iqout][iqside][iqlong];
                        double q_out_local, q_side_local, q_long_local;
                        double correl_fun_num, correl_fun_denorm;
                        double correl_fun_val;
                        if(npart_num < 2 || npart_denorm < 2)
                        {
                            q_out_local = q_out[iqout];
                            q_side_local = q_side[iqside];
                            q_long_local = q_long[iqlong];
                            correl_fun_num = 0.0;
                            correl_fun_denorm = npart_denorm;
                            correl_fun_val = 0.0;
                        }
                        else
                        {
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

                        output << scientific << setw(18) << setprecision(8) 
                               << q_out_local << "    " 
                               << q_side_local << "    "
                               << q_long_local << "    "
                               << npart_num << "    "
                               << correl_fun_num << "    " 
                               << correl_fun_denorm << "    "
                               << correl_fun_val << "    " << 0.0 << endl;
                    }
                }
            }
            output.close();
        }
    }
}
