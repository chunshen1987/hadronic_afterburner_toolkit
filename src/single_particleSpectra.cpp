// Copyright Chun Shen @ 2016
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <iomanip>
#include "./parameters.h"
#include "./single_particleSpectra.h"

using namespace std;

singleParticleSpectra::singleParticleSpectra(
        ParameterReader *paraRdr_in, string path_in, 
        particleSamples *particle_list_in) {
    paraRdr = paraRdr_in;
    path = path_in;
    particle_list = particle_list_in;

    particle_monval = paraRdr->getVal("particle_monval");

    if (particle_monval == 333) {
        // phi(1020) is reconstructed from (K^+, K^-) pairs
        reconst_branching_ratio = 0.489;
    } else {
        reconst_branching_ratio = 1.0;
    }

    order_max = paraRdr->getVal("order_max");
    Qn_vector_real = new double[order_max];
    Qn_vector_imag = new double[order_max];
    Qn_vector_real_err = new double[order_max];
    Qn_vector_imag_err = new double[order_max];
    Qn_diff_vector_real = new double* [order_max];
    Qn_diff_vector_imag = new double* [order_max];
    Qn_diff_vector_real_err = new double* [order_max];
    Qn_diff_vector_imag_err = new double* [order_max];

    npT = paraRdr->getVal("npT");
    pT_min = paraRdr->getVal("pT_min");
    pT_max = paraRdr->getVal("pT_max");
    dpT = (pT_max - pT_min)/(npT - 1 + 1e-15);
    pT_array = new double [npT];
    pT_mean_array = new double [npT];
    pT_mean_array_err = new double [npT];
    for (int i = 0; i < npT; i++) {
        pT_array[i] = pT_min + dpT*i;
        pT_mean_array[i] = 0.0;
        pT_mean_array_err[i] = 0.0;
    }
    for (int i = 0; i < order_max; i++) {
        Qn_vector_real[i] = 0.0;
        Qn_vector_imag[i] = 0.0;
        Qn_vector_real_err[i] = 0.0;
        Qn_vector_imag_err[i] = 0.0;
        Qn_diff_vector_real[i] = new double[npT];
        Qn_diff_vector_imag[i] = new double[npT];
        Qn_diff_vector_real_err[i] = new double[npT];
        Qn_diff_vector_imag_err[i] = new double[npT];
        for (int j = 0; j < npT; j++) {
            Qn_diff_vector_real[i][j] = 0.0;
            Qn_diff_vector_imag[i][j] = 0.0;
            Qn_diff_vector_real_err[i][j] = 0.0;
            Qn_diff_vector_imag_err[i][j] = 0.0;
        }
    }
    total_number_of_events = 0;

    rap_type = paraRdr->getVal("rap_type");
    rap_min = paraRdr->getVal("rap_min");
    rap_max = paraRdr->getVal("rap_max");

    if (particle_monval == 9999)  // use pseudo-rapidity for all charged hadrons
        rap_type = 0;

    rapidity_distribution_flag = paraRdr->getVal("rapidity_distribution");
    if (rapidity_distribution_flag == 1) {
        N_rap = paraRdr->getVal("n_rap");
        rapidity_dis_min = paraRdr->getVal("rapidity_dis_min");
        rapidity_dis_max = paraRdr->getVal("rapidity_dis_max");
        drap = (rapidity_dis_max - rapidity_dis_min)/(N_rap - 1.);
        rapidity_array = new double [N_rap];
        dNdy_array = new double [N_rap];
        for (int i = 0; i < N_rap; i++) {
            rapidity_array[i] = rapidity_dis_min + i*drap;
            dNdy_array[i] = 0.0;
        }
        vn_rapidity_dis_pT_min = paraRdr->getVal("vn_rapidity_dis_pT_min");
        vn_rapidity_dis_pT_max = paraRdr->getVal("vn_rapidity_dis_pT_max");
        vn_real_rapidity_dis_array = new double* [N_rap];
        vn_imag_rapidity_dis_array = new double* [N_rap];
        vn_real_rapidity_dis_array_err = new double* [N_rap];
        vn_imag_rapidity_dis_array_err = new double* [N_rap];
        for (int i = 0; i < N_rap; i++) {
            vn_real_rapidity_dis_array[i] = new double [order_max];
            vn_imag_rapidity_dis_array[i] = new double [order_max];
            vn_real_rapidity_dis_array_err[i] = new double [order_max];
            vn_imag_rapidity_dis_array_err[i] = new double [order_max];
            for (int j = 0; j < order_max; j++) {
                vn_real_rapidity_dis_array[i][j] = 0.0;
                vn_imag_rapidity_dis_array[i][j] = 0.0;
                vn_real_rapidity_dis_array_err[i][j] = 0.0;
                vn_imag_rapidity_dis_array_err[i][j] = 0.0;
            }
        }
    }

    // check dN/dtau distribution
    check_spatial_flag = paraRdr->getVal("check_spatial_dis");
    if (check_spatial_flag == 1) {
        // dN/dtau
        intrinsic_dtau = paraRdr->getVal("intrinsic_dtau");
        N_tau = 50;
        tau_min = 0.6;
        tau_max = 15.0;
        dtau = (tau_max - tau_min)/(N_tau - 1);
        tau_array = new double [N_tau];
        dNdtau_array = new double [N_tau];
        for (int i = 0; i < N_tau; i++) {
            tau_array[i] = tau_min + i*dtau;
            dNdtau_array[i] = 0.0;
        }

        // dN/dx
        intrinsic_dx = paraRdr->getVal("intrinsic_dx");
        N_xpt = 50;
        spatial_x_min = -10.0;
        spatial_x_max = 10.0;
        dspatial_x = (spatial_x_max - spatial_x_min)/(N_xpt - 1);
        xpt_array = new double [N_xpt];
        ypt_array = new double [N_xpt];
        dNdx1_array = new double [N_xpt];
        dNdx2_array = new double [N_xpt];
        for (int i = 0; i < N_xpt; i++) {
            xpt_array[i] = spatial_x_min + i*dspatial_x;
            ypt_array[i] = spatial_x_min + i*dspatial_x;
            dNdx1_array[i] = 0.0;
            dNdx2_array[i] = 0.0;
        }

        // dN/(dtau dx)
        dNdtaudx1_array = new double * [N_tau];
        dNdtaudx2_array = new double * [N_tau];
        for(int i = 0; i < N_tau; i++) {
            dNdtaudx1_array[i] = new double [N_xpt];
            dNdtaudx2_array[i] = new double [N_xpt];
            for (int j = 0; j < N_xpt; j++) {
                dNdtaudx1_array[i][j] = 0.0;
                dNdtaudx2_array[i][j] = 0.0;
            }
        }

        // dN/deta_s
        intrinsic_detas = paraRdr->getVal("intrinsic_detas");
        N_eta_s = 40;
        eta_s_min = - 3.0;
        eta_s_max = 3.0;
        deta_s = (eta_s_max - eta_s_min)/(N_eta_s - 1);
        eta_s_array = new double [N_eta_s];
        dNdetas_array = new double [N_eta_s];
        for (int i = 0; i < N_eta_s; i++) {
            eta_s_array[i] = eta_s_min + i*deta_s;
            dNdetas_array[i] = 0.0;
        }
    }

    flag_correlation = paraRdr->getVal("compute_correlation");
    if (flag_correlation == 1) {
        Qn2_vector = new double[order_max];
        Qn2_vector_err = new double[order_max];
        QnSP_diff_vector = new double* [order_max];
        QnSP_diff_vector_err = new double* [order_max];
        QnSP_eta12 = new double* [order_max];
        QnSP_eta12_err = new double* [order_max];
        for (int i = 0; i < order_max; i++) {
            Qn2_vector[i] = 0.0;
            Qn2_vector_err[i] = 0.0;
            QnSP_diff_vector[i] = new double[npT];
            QnSP_diff_vector_err[i] = new double[npT];
            for (int j = 0; j < npT; j++) {
                QnSP_diff_vector[i][j] = 0.0;
                QnSP_diff_vector_err[i][j] = 0.0;
            }
            QnSP_eta12[i] = new double[N_rap];
            QnSP_eta12_err[i] = new double[N_rap];
            for (int j = 0; j < N_rap; j++) {
                QnSP_eta12[i][j] = 0.0;
                QnSP_eta12_err[i][j] = 0.0;
            }
        }

        num_corr = 9;
        C_nmk = new double[num_corr];
        C_nmk_err = new double[num_corr];
        C_nmk_eta12 = new double* [num_corr];
        C_nmk_eta12_err = new double* [num_corr];
        C_nmk_eta13 = new double* [num_corr];
        C_nmk_eta13_err = new double* [num_corr];
        for (int i = 0; i < num_corr; i++) {
            C_nmk[i] = 0.0;
            C_nmk_err[i] = 0.0;
            C_nmk_eta12[i] = new double[N_rap];
            C_nmk_eta12_err[i] = new double[N_rap];
            C_nmk_eta13[i] = new double[N_rap];
            C_nmk_eta13_err[i] = new double[N_rap];
            for (int j = 0; j < N_rap; j++) {
                C_nmk_eta12[i][j] = 0.0;
                C_nmk_eta12_err[i][j] = 0.0;
                C_nmk_eta13[i][j] = 0.0;
                C_nmk_eta13_err[i][j] = 0.0;
            }
        }
        flag_charge_dependence = paraRdr->getVal("flag_charge_dependence");
        if (flag_charge_dependence == 1) {
            C_nmk_ss = new double[num_corr];
            C_nmk_ss_err = new double[num_corr];
            C_nmk_os = new double[num_corr];
            C_nmk_os_err = new double[num_corr];
            C_nmk_eta12_ss = new double* [num_corr];
            C_nmk_eta12_ss_err = new double* [num_corr];
            C_nmk_eta12_os = new double* [num_corr];
            C_nmk_eta12_os_err = new double* [num_corr];
            C_nmk_eta13_ss = new double* [num_corr];
            C_nmk_eta13_ss_err = new double* [num_corr];
            C_nmk_eta13_os = new double* [num_corr];
            C_nmk_eta13_os_err = new double* [num_corr];
            for (int i = 0; i < num_corr; i++) {
                C_nmk_ss[i] = 0.0;
                C_nmk_ss_err[i] = 0.0;
                C_nmk_os[i] = 0.0;
                C_nmk_os_err[i] = 0.0;
                C_nmk_eta12_ss[i] = new double[N_rap];
                C_nmk_eta12_ss_err[i] = new double[N_rap];
                C_nmk_eta12_os[i] = new double[N_rap];
                C_nmk_eta12_os_err[i] = new double[N_rap];
                C_nmk_eta13_ss[i] = new double[N_rap];
                C_nmk_eta13_ss_err[i] = new double[N_rap];
                C_nmk_eta13_os[i] = new double[N_rap];
                C_nmk_eta13_os_err[i] = new double[N_rap];
                for (int j = 0; j < N_rap; j++) {
                    C_nmk_eta12_ss[i][j] = 0.0;
                    C_nmk_eta12_ss_err[i][j] = 0.0;
                    C_nmk_eta12_os[i][j] = 0.0;
                    C_nmk_eta12_os_err[i][j] = 0.0;
                    C_nmk_eta13_ss[i][j] = 0.0;
                    C_nmk_eta13_ss_err[i][j] = 0.0;
                    C_nmk_eta13_os[i][j] = 0.0;
                    C_nmk_eta13_os_err[i][j] = 0.0;
                }
            }
        }
        SC_num_corr = 6;
        SC_mn = new double[SC_num_corr];
        SC_mn_err = new double[SC_num_corr];
        for (int i = 0; i < SC_num_corr; i++) {
            SC_mn[i] = 0.0;
            SC_mn_err[i] = 0.0;
        }
    }
}

singleParticleSpectra::~singleParticleSpectra() {
    delete [] pT_array;
    delete [] pT_mean_array;
    delete [] pT_mean_array_err;
    delete [] Qn_vector_real;
    delete [] Qn_vector_imag;
    delete [] Qn_vector_real_err;
    delete [] Qn_vector_imag_err;
    for (int i = 0; i < order_max; i++) {
        delete [] Qn_diff_vector_real[i];
        delete [] Qn_diff_vector_imag[i];
        delete [] Qn_diff_vector_real_err[i];
        delete [] Qn_diff_vector_imag_err[i];
    }
    delete [] Qn_diff_vector_real;
    delete [] Qn_diff_vector_imag;
    delete [] Qn_diff_vector_real_err;
    delete [] Qn_diff_vector_imag_err;

    if (rapidity_distribution_flag == 1) {
        delete [] rapidity_array;
        delete [] dNdy_array;
        for (int i = 0; i < N_rap; i++) {
            delete [] vn_real_rapidity_dis_array[i];
            delete [] vn_imag_rapidity_dis_array[i];
            delete [] vn_real_rapidity_dis_array_err[i];
            delete [] vn_imag_rapidity_dis_array_err[i];
        }
        delete [] vn_real_rapidity_dis_array;
        delete [] vn_imag_rapidity_dis_array;
        delete [] vn_real_rapidity_dis_array_err;
        delete [] vn_imag_rapidity_dis_array_err;
    }

    if (check_spatial_flag == 1) {
        delete [] tau_array;
        delete [] dNdtau_array;
        delete [] xpt_array;
        delete [] ypt_array;
        delete [] dNdx1_array;
        delete [] dNdx2_array;
        delete [] eta_s_array;
        delete [] dNdetas_array;
        for (int i = 0; i < N_tau; i++) {
            delete [] dNdtaudx1_array[i];
            delete [] dNdtaudx2_array[i];
        }
        delete [] dNdtaudx1_array;
        delete [] dNdtaudx2_array;
    }

    if (flag_correlation == 1) {
        for (int i = 0; i < order_max; i++) {
            delete[] QnSP_diff_vector[i];
            delete[] QnSP_diff_vector_err[i];
            delete[] QnSP_eta12[i];
            delete[] QnSP_eta12_err[i];
        }
        delete[] QnSP_diff_vector;
        delete[] QnSP_diff_vector_err;
        delete[] QnSP_eta12;
        delete[] QnSP_eta12_err;
        delete[] Qn2_vector;
        delete[] Qn2_vector_err;

        delete[] C_nmk;
        delete[] C_nmk_err;
        for (int i = 0; i < num_corr; i++) {
            delete[] C_nmk_eta12[i];
            delete[] C_nmk_eta12_err[i];
            delete[] C_nmk_eta13[i];
            delete[] C_nmk_eta13_err[i];
        }
        delete[] C_nmk_eta12;
        delete[] C_nmk_eta12_err;
        delete[] C_nmk_eta13;
        delete[] C_nmk_eta13_err;
        if (flag_charge_dependence == 1) {
            delete[] C_nmk_ss;
            delete[] C_nmk_ss_err;
            delete[] C_nmk_os;
            delete[] C_nmk_os_err;
            for (int i = 0; i < num_corr; i++) {
                delete[] C_nmk_eta12_ss[i];
                delete[] C_nmk_eta12_ss_err[i];
                delete[] C_nmk_eta13_ss[i];
                delete[] C_nmk_eta13_ss_err[i];
                delete[] C_nmk_eta12_os[i];
                delete[] C_nmk_eta12_os_err[i];
                delete[] C_nmk_eta13_os[i];
                delete[] C_nmk_eta13_os_err[i];
            }
            delete[] C_nmk_eta12_ss;
            delete[] C_nmk_eta12_ss_err;
            delete[] C_nmk_eta13_ss;
            delete[] C_nmk_eta13_ss_err;
            delete[] C_nmk_eta12_os;
            delete[] C_nmk_eta12_os_err;
            delete[] C_nmk_eta13_os;
            delete[] C_nmk_eta13_os_err;
        }
        delete[] SC_mn;
        delete[] SC_mn_err;
    }
}

//! This is a driver function to compute the Qn flow vector
void singleParticleSpectra::calculate_Qn_vector_shell() {
    int event_id = 0;
    int buffer_size = particle_list->get_event_buffer_size();
    // initialize some temp arrays
    double *event_Qn_real = new double[order_max];
    double *event_Qn_real_err = new double[order_max];
    double *event_Qn_imag = new double[order_max];
    double *event_Qn_imag_err = new double[order_max];
    double **event_Qn_diff_real = new double* [order_max];
    double **event_Qn_diff_real_err = new double* [order_max];
    double **event_Qn_diff_imag = new double* [order_max];
    double **event_Qn_diff_imag_err = new double* [order_max];
    for (int i = 0; i < order_max; i++) {
        event_Qn_real[i] = 0.0;
        event_Qn_real_err[i] = 0.0;
        event_Qn_imag[i] = 0.0;
        event_Qn_imag_err[i] = 0.0;
        event_Qn_diff_real[i] = new double[npT];
        event_Qn_diff_real_err[i] = new double[npT];
        event_Qn_diff_imag[i] = new double[npT];
        event_Qn_diff_imag_err[i] = new double[npT];
        for (int j = 0; j < npT; j++) {
            event_Qn_diff_real[i][j] = 0.0;
            event_Qn_diff_real_err[i][j] = 0.0;
            event_Qn_diff_imag[i][j] = 0.0;
            event_Qn_diff_imag_err[i][j] = 0.0;
        }
    }
    double *event_Qn_p_real, *event_Qn_p_real_err;
    double *event_Qn_p_imag, *event_Qn_p_imag_err;
    double *event_Qn_m_real, *event_Qn_m_real_err;
    double *event_Qn_m_imag, *event_Qn_m_imag_err;
    double **event_Qn_p_diff_real, **event_Qn_p_diff_real_err;
    double **event_Qn_p_diff_imag, **event_Qn_p_diff_imag_err;
    double **event_Qn_m_diff_real, **event_Qn_m_diff_real_err;
    double **event_Qn_m_diff_imag, **event_Qn_m_diff_imag_err;
    if (flag_charge_dependence == 1) {
        event_Qn_p_real = new double[order_max];
        event_Qn_p_real_err = new double[order_max];
        event_Qn_p_imag = new double[order_max];
        event_Qn_p_imag_err = new double[order_max];
        event_Qn_m_real = new double[order_max];
        event_Qn_m_real_err = new double[order_max];
        event_Qn_m_imag = new double[order_max];
        event_Qn_m_imag_err = new double[order_max];
        event_Qn_p_diff_real = new double* [order_max];
        event_Qn_p_diff_real_err = new double* [order_max];
        event_Qn_p_diff_imag = new double* [order_max];
        event_Qn_p_diff_imag_err = new double* [order_max];
        event_Qn_m_diff_real = new double* [order_max];
        event_Qn_m_diff_real_err = new double* [order_max];
        event_Qn_m_diff_imag = new double* [order_max];
        event_Qn_m_diff_imag_err = new double* [order_max];
        for (int i = 0; i < order_max; i++) {
            event_Qn_p_diff_real[i] = new double[npT];
            event_Qn_p_diff_real_err[i] = new double[npT];
            event_Qn_p_diff_imag[i] = new double[npT];
            event_Qn_p_diff_imag_err[i] = new double[npT];
            event_Qn_m_diff_real[i] = new double[npT];
            event_Qn_m_diff_real_err[i] = new double[npT];
            event_Qn_m_diff_imag[i] = new double[npT];
            event_Qn_m_diff_imag_err[i] = new double[npT];
        }
    }
    double **event_Qn_rap_real = new double* [N_rap];
    double **event_Qn_rap_real_err = new double* [N_rap];
    double **event_Qn_rap_imag = new double* [N_rap];
    double **event_Qn_rap_imag_err = new double* [N_rap];
    for (int i = 0; i < N_rap; i++) {
        event_Qn_rap_real[i] = new double[order_max];
        event_Qn_rap_real_err[i] = new double[order_max];
        event_Qn_rap_imag[i] = new double[order_max];
        event_Qn_rap_imag_err[i] = new double[order_max];
        for (int j = 0; j < order_max; j++) {
            event_Qn_rap_real[i][j] = 0.0;
            event_Qn_rap_real_err[i][j] = 0.0;
            event_Qn_rap_imag[i][j] = 0.0;
            event_Qn_rap_imag_err[i][j] = 0.0;
        }
    }
    double **event_Qn_p_rap_real, **event_Qn_p_rap_real_err;
    double **event_Qn_p_rap_imag, **event_Qn_p_rap_imag_err;
    double **event_Qn_m_rap_real, **event_Qn_m_rap_real_err;
    double **event_Qn_m_rap_imag, **event_Qn_m_rap_imag_err;
    if (flag_charge_dependence == 1) {
        event_Qn_p_rap_real = new double* [N_rap];
        event_Qn_p_rap_real_err = new double* [N_rap];
        event_Qn_p_rap_imag = new double* [N_rap];
        event_Qn_p_rap_imag_err = new double* [N_rap];
        event_Qn_m_rap_real = new double* [N_rap];
        event_Qn_m_rap_real_err = new double* [N_rap];
        event_Qn_m_rap_imag = new double* [N_rap];
        event_Qn_m_rap_imag_err = new double* [N_rap];
        for (int i = 0; i < N_rap; i++) {
            event_Qn_p_rap_real[i] = new double[order_max];
            event_Qn_p_rap_real_err[i] = new double[order_max];
            event_Qn_p_rap_imag[i] = new double[order_max];
            event_Qn_p_rap_imag_err[i] = new double[order_max];
            event_Qn_m_rap_real[i] = new double[order_max];
            event_Qn_m_rap_real_err[i] = new double[order_max];
            event_Qn_m_rap_imag[i] = new double[order_max];
            event_Qn_m_rap_imag_err[i] = new double[order_max];
            for (int j = 0; j < order_max; j++) {
                event_Qn_p_rap_real[i][j] = 0.0;
                event_Qn_p_rap_real_err[i][j] = 0.0;
                event_Qn_p_rap_imag[i][j] = 0.0;
                event_Qn_p_rap_imag_err[i][j] = 0.0;
                event_Qn_m_rap_real[i][j] = 0.0;
                event_Qn_m_rap_real_err[i][j] = 0.0;
                event_Qn_m_rap_imag[i][j] = 0.0;
                event_Qn_m_rap_imag_err[i][j] = 0.0;
            }
        }
    }

    // start the loop
    while (!particle_list->end_of_file()) {
        cout << "Reading event: " << event_id + 1 
             << "-" << event_id + buffer_size << " ... " << flush;
        particle_list->read_in_particle_samples();
        cout << " processing ..." << flush;
        int nev = particle_list->get_number_of_events();
        for (int iev = 0; iev < nev; iev++) {
            event_id++;
            calculate_Qn_vector(iev,
                                event_Qn_real, event_Qn_real_err,
                                event_Qn_imag, event_Qn_imag_err,
                                event_Qn_diff_real, event_Qn_diff_real_err,
                                event_Qn_diff_imag, event_Qn_diff_imag_err);
            for (int i = 0; i < order_max; i++) {
                Qn_vector_real[i] += event_Qn_real[i];
                Qn_vector_real_err[i] += event_Qn_real_err[i];
                Qn_vector_imag[i] += event_Qn_imag[i];
                Qn_vector_imag_err[i] += event_Qn_imag_err[i];
                for (int j = 0; j < npT; j++) {
                    Qn_diff_vector_real[i][j] += event_Qn_diff_real[i][j];
                    Qn_diff_vector_real_err[i][j] += (
                                                 event_Qn_diff_real_err[i][j]);
                    Qn_diff_vector_imag[i][j] += event_Qn_diff_imag[i][j];
                    Qn_diff_vector_imag_err[i][j] += (
                                                 event_Qn_diff_imag_err[i][j]);
                }
            }

            if (flag_correlation == 1) {
                calculate_two_particle_correlation(
                        event_Qn_real, event_Qn_imag,
                        event_Qn_diff_real, event_Qn_diff_imag);
                calculate_three_particle_correlation(
                        event_Qn_real, event_Qn_imag,
                        event_Qn_real, event_Qn_imag,
                        event_Qn_real, event_Qn_imag, 0, C_nmk, C_nmk_err);
                calculate_four_particle_correlation_SC(
                        event_Qn_real, event_Qn_imag,
                        event_Qn_real, event_Qn_imag,
                        event_Qn_real, event_Qn_imag,
                        event_Qn_real, event_Qn_imag, 0, SC_mn, SC_mn_err);
                if (flag_charge_dependence == 1) {
                    calculate_Qn_vector_positive_charge(iev,
                        event_Qn_p_real, event_Qn_p_real_err,
                        event_Qn_p_imag, event_Qn_p_imag_err,
                        event_Qn_p_diff_real, event_Qn_p_diff_real_err,
                        event_Qn_p_diff_imag, event_Qn_p_diff_imag_err);
                    calculate_Qn_vector_negative_charge(iev,
                        event_Qn_m_real, event_Qn_m_real_err,
                        event_Qn_m_imag, event_Qn_m_imag_err,
                        event_Qn_m_diff_real, event_Qn_m_diff_real_err,
                        event_Qn_m_diff_imag, event_Qn_m_diff_imag_err);
                    calculate_three_particle_correlation(
                            event_Qn_p_real, event_Qn_p_imag,
                            event_Qn_p_real, event_Qn_p_imag,
                            event_Qn_real, event_Qn_imag, 0,
                            C_nmk_ss, C_nmk_ss_err);
                    calculate_three_particle_correlation(
                            event_Qn_m_real, event_Qn_m_imag,
                            event_Qn_m_real, event_Qn_m_imag,
                            event_Qn_real, event_Qn_imag, 0,
                            C_nmk_ss, C_nmk_ss_err);
                    calculate_three_particle_correlation(
                            event_Qn_p_real, event_Qn_p_imag,
                            event_Qn_m_real, event_Qn_m_imag,
                            event_Qn_real, event_Qn_imag, 1,
                            C_nmk_os, C_nmk_os_err);
                    calculate_three_particle_correlation(
                            event_Qn_m_real, event_Qn_m_imag,
                            event_Qn_p_real, event_Qn_p_imag,
                            event_Qn_real, event_Qn_imag, 1,
                            C_nmk_os, C_nmk_os_err);
                }
            }

            if (rapidity_distribution_flag == 1) {
                calculate_rapidity_distribution(iev,
                        event_Qn_rap_real, event_Qn_rap_real_err,
                        event_Qn_rap_imag, event_Qn_rap_imag_err, 0);
                for (int i = 0; i < N_rap; i++) {
                    for (int j = 0; j < order_max; j++) {
                        vn_real_rapidity_dis_array[i][j] += (
                                                event_Qn_rap_real[i][j]);
                        vn_real_rapidity_dis_array_err[i][j] += (
                                                event_Qn_rap_real_err[i][j]);
                        vn_imag_rapidity_dis_array[i][j] += (
                                                event_Qn_rap_imag[i][j]);
                        vn_imag_rapidity_dis_array_err[i][j] += (
                                                event_Qn_rap_imag_err[i][j]);
                    }
                }
                if (flag_correlation == 1) {
                    calculate_two_particle_correlation_deltaeta(
                            event_Qn_rap_real, event_Qn_rap_imag);
                    calculate_three_particle_correlation_deltaeta(
                            event_Qn_rap_real, event_Qn_rap_imag,
                            event_Qn_rap_real, event_Qn_rap_imag,
                            event_Qn_rap_real, event_Qn_rap_imag, 1, 1,
                            C_nmk_eta12, C_nmk_eta12_err);
                    calculate_three_particle_correlation_deltaeta(
                            event_Qn_rap_real, event_Qn_rap_imag,
                            event_Qn_rap_real, event_Qn_rap_imag,
                            event_Qn_rap_real, event_Qn_rap_imag, 2, 1,
                            C_nmk_eta13, C_nmk_eta13_err);
                    if (flag_charge_dependence == 1) {
                        calculate_rapidity_distribution(iev,
                                event_Qn_p_rap_real, event_Qn_p_rap_real_err,
                                event_Qn_p_rap_imag, event_Qn_p_rap_imag_err,
                                1);
                        calculate_rapidity_distribution(iev,
                                event_Qn_m_rap_real, event_Qn_m_rap_real_err,
                                event_Qn_m_rap_imag, event_Qn_m_rap_imag_err,
                                2);
                        calculate_three_particle_correlation_deltaeta(
                                event_Qn_p_rap_real, event_Qn_p_rap_imag,
                                event_Qn_m_rap_real, event_Qn_m_rap_imag,
                                event_Qn_rap_real, event_Qn_rap_imag, 1, 2,
                                C_nmk_eta12_os, C_nmk_eta12_os_err);
                        calculate_three_particle_correlation_deltaeta(
                                event_Qn_m_rap_real, event_Qn_m_rap_imag,
                                event_Qn_p_rap_real, event_Qn_p_rap_imag,
                                event_Qn_rap_real, event_Qn_rap_imag, 1, 2,
                                C_nmk_eta12_os, C_nmk_eta12_os_err);
                        calculate_three_particle_correlation_deltaeta(
                                event_Qn_p_rap_real, event_Qn_p_rap_imag,
                                event_Qn_m_rap_real, event_Qn_m_rap_imag,
                                event_Qn_rap_real, event_Qn_rap_imag, 2, 2,
                                C_nmk_eta13_os, C_nmk_eta13_os_err);
                        calculate_three_particle_correlation_deltaeta(
                                event_Qn_m_rap_real, event_Qn_m_rap_imag,
                                event_Qn_p_rap_real, event_Qn_p_rap_imag,
                                event_Qn_rap_real, event_Qn_rap_imag, 2, 2,
                                C_nmk_eta13_os, C_nmk_eta13_os_err);
                        calculate_three_particle_correlation_deltaeta(
                                event_Qn_p_rap_real, event_Qn_p_rap_imag,
                                event_Qn_p_rap_real, event_Qn_p_rap_imag,
                                event_Qn_rap_real, event_Qn_rap_imag, 1, 1,
                                C_nmk_eta12_ss, C_nmk_eta12_ss_err);
                        calculate_three_particle_correlation_deltaeta(
                                event_Qn_m_rap_real, event_Qn_m_rap_imag,
                                event_Qn_m_rap_real, event_Qn_m_rap_imag,
                                event_Qn_rap_real, event_Qn_rap_imag, 1, 1,
                                C_nmk_eta12_ss, C_nmk_eta12_ss_err);
                        calculate_three_particle_correlation_deltaeta(
                                event_Qn_p_rap_real, event_Qn_p_rap_imag,
                                event_Qn_p_rap_real, event_Qn_p_rap_imag,
                                event_Qn_rap_real, event_Qn_rap_imag, 2, 1,
                                C_nmk_eta13_ss, C_nmk_eta13_ss_err);
                        calculate_three_particle_correlation_deltaeta(
                                event_Qn_m_rap_real, event_Qn_m_rap_imag,
                                event_Qn_m_rap_real, event_Qn_m_rap_imag,
                                event_Qn_rap_real, event_Qn_rap_imag, 2, 1,
                                C_nmk_eta13_ss, C_nmk_eta13_ss_err);
                    }
                }
            }

            if (check_spatial_flag == 1)
                check_dNdSV(iev);
        }
        cout << " done!" << endl;
    }
    total_number_of_events = event_id;
    output_Qn_vectors();
    if (rapidity_distribution_flag == 1) {
        output_rapidity_distribution();
        if (flag_correlation == 1) {
            output_two_particle_correlation_rap();
            output_three_particle_correlation_rap();
        }
    }
    if (check_spatial_flag == 1)
        output_dNdSV();
    if (flag_correlation == 1) {
        output_two_particle_correlation();
        output_three_particle_correlation();
        output_four_particle_SC_correlation();
    }

    // clean up
    delete[] event_Qn_real;
    delete[] event_Qn_real_err;
    delete[] event_Qn_imag;
    delete[] event_Qn_imag_err;
    for (int i = 0; i < order_max; i++) {
        delete[] event_Qn_diff_real[i];
        delete[] event_Qn_diff_real_err[i];
        delete[] event_Qn_diff_imag[i];
        delete[] event_Qn_diff_imag_err[i];
    }
    for (int i = 0; i < N_rap; i++) {
        delete[] event_Qn_rap_real[i];
        delete[] event_Qn_rap_real_err[i];
        delete[] event_Qn_rap_imag[i];
        delete[] event_Qn_rap_imag_err[i];
    }
    delete[] event_Qn_diff_real;
    delete[] event_Qn_diff_real_err;
    delete[] event_Qn_diff_imag;
    delete[] event_Qn_diff_imag_err;
    delete[] event_Qn_rap_real;
    delete[] event_Qn_rap_real_err;
    delete[] event_Qn_rap_imag;
    delete[] event_Qn_rap_imag_err;
    if (flag_charge_dependence == 1) {
        delete[] event_Qn_p_real;
        delete[] event_Qn_p_real_err;
        delete[] event_Qn_p_imag;
        delete[] event_Qn_p_imag_err;
        delete[] event_Qn_m_real;
        delete[] event_Qn_m_real_err;
        delete[] event_Qn_m_imag;
        delete[] event_Qn_m_imag_err;
        for (int i = 0; i < order_max; i++) {
            delete[] event_Qn_p_diff_real[i];
            delete[] event_Qn_p_diff_real_err[i];
            delete[] event_Qn_p_diff_imag[i];
            delete[] event_Qn_p_diff_imag_err[i];
            delete[] event_Qn_m_diff_real[i];
            delete[] event_Qn_m_diff_real_err[i];
            delete[] event_Qn_m_diff_imag[i];
            delete[] event_Qn_m_diff_imag_err[i];
        }
        delete[] event_Qn_p_diff_real;
        delete[] event_Qn_p_diff_real_err;
        delete[] event_Qn_p_diff_imag;
        delete[] event_Qn_p_diff_imag_err;
        delete[] event_Qn_m_diff_real;
        delete[] event_Qn_m_diff_real_err;
        delete[] event_Qn_m_diff_imag;
        delete[] event_Qn_m_diff_imag_err;
        for (int i = 0; i < N_rap; i++) {
            delete[] event_Qn_p_rap_real[i];
            delete[] event_Qn_p_rap_real_err[i];
            delete[] event_Qn_p_rap_imag[i];
            delete[] event_Qn_p_rap_imag_err[i];
            delete[] event_Qn_m_rap_real[i];
            delete[] event_Qn_m_rap_real_err[i];
            delete[] event_Qn_m_rap_imag[i];
            delete[] event_Qn_m_rap_imag_err[i];
        }
        delete[] event_Qn_p_rap_real;
        delete[] event_Qn_p_rap_real_err;
        delete[] event_Qn_p_rap_imag;
        delete[] event_Qn_p_rap_imag_err;
        delete[] event_Qn_m_rap_real;
        delete[] event_Qn_m_rap_real_err;
        delete[] event_Qn_m_rap_imag;
        delete[] event_Qn_m_rap_imag_err;
    }
}

//! this function computes the pT-integrated and pT-differential Qn vector
//! within a given rapidity region in one event
void singleParticleSpectra::calculate_Qn_vector(int event_id,
        double *event_Qn_real, double *event_Qn_real_err,
        double *event_Qn_imag, double *event_Qn_imag_err,
        double **event_Qn_diff_real, double **event_Qn_diff_real_err,
        double **event_Qn_diff_imag, double **event_Qn_diff_imag_err) {
    
    // first clean the results arrays
    for (int i = 0; i < order_max; i++) {
        event_Qn_real[i] = 0.0;
        event_Qn_real_err[i] = 0.0;
        event_Qn_imag[i] = 0.0;
        event_Qn_imag_err[i] = 0.0;
        for (int j = 0; j < npT; j++) {
            event_Qn_diff_real[i][j] = 0.0;
            event_Qn_diff_real_err[i][j] = 0.0;
            event_Qn_diff_imag[i][j] = 0.0;
            event_Qn_diff_imag_err[i][j] = 0.0;
        }
    }

    int number_of_particles = particle_list->get_number_of_particles(event_id);
    for (int i = 0; i < number_of_particles; i++) {
        double pz_local = particle_list->get_particle(event_id, i).pz;
        double E_local = particle_list->get_particle(event_id, i).E;

        double rap_local;
        if (rap_type == 0) {
            double mass = particle_list->get_particle(event_id, i).mass;
            double pmag = sqrt(E_local*E_local - mass*mass);
            rap_local = 0.5*log((pmag + pz_local)/(pmag - pz_local));
        } else {
            rap_local = 0.5*log((E_local + pz_local)/(E_local - pz_local));
        }

        if (rap_local > rap_min && rap_local < rap_max) {
            double px_local = particle_list->get_particle(event_id, i).px;
            double py_local = particle_list->get_particle(event_id, i).py;
            double p_perp = sqrt(px_local*px_local + py_local*py_local);
            double p_phi = atan2(py_local, px_local);
            for (int iorder = 0; iorder < order_max; iorder++) {
                double cos_nphi = cos(iorder*p_phi);
                double sin_nphi = sin(iorder*p_phi);
                event_Qn_real[iorder] += cos_nphi;
                event_Qn_imag[iorder] += sin_nphi;
                event_Qn_real_err[iorder] += cos_nphi*cos_nphi;
                event_Qn_imag_err[iorder] += sin_nphi*sin_nphi;
            }
            if (p_perp > pT_min && p_perp < pT_max) {
                int p_idx = static_cast<int>((p_perp - pT_min)/dpT);
                pT_mean_array[p_idx] += p_perp;
                pT_mean_array_err[p_idx] += p_perp*p_perp;
                for (int iorder = 0; iorder < order_max; iorder++) {
                    double cos_nphi = cos(iorder*p_phi);
                    double sin_nphi = sin(iorder*p_phi);
                    event_Qn_diff_real[iorder][p_idx] += cos_nphi;
                    event_Qn_diff_imag[iorder][p_idx] += sin_nphi;
                    event_Qn_diff_real_err[iorder][p_idx] += cos_nphi*cos_nphi;
                    event_Qn_diff_imag_err[iorder][p_idx] += sin_nphi*sin_nphi;
                }
            }
        }
    }
}


//! this function computes the pT-integrated and pT-differential Qn vector
//! of postive charged hadrons within a given rapidity region in one event
void singleParticleSpectra::calculate_Qn_vector_positive_charge(int event_id,
        double *event_Qn_real, double *event_Qn_real_err,
        double *event_Qn_imag, double *event_Qn_imag_err,
        double **event_Qn_diff_real, double **event_Qn_diff_real_err,
        double **event_Qn_diff_imag, double **event_Qn_diff_imag_err) {
    
    // first clean the results arrays
    for (int i = 0; i < order_max; i++) {
        event_Qn_real[i] = 0.0;
        event_Qn_real_err[i] = 0.0;
        event_Qn_imag[i] = 0.0;
        event_Qn_imag_err[i] = 0.0;
        for (int j = 0; j < npT; j++) {
            event_Qn_diff_real[i][j] = 0.0;
            event_Qn_diff_real_err[i][j] = 0.0;
            event_Qn_diff_imag[i][j] = 0.0;
            event_Qn_diff_imag_err[i][j] = 0.0;
        }
    }

    int number_of_particles = (
                    particle_list->get_number_of_positive_particles(event_id));
    for (int i = 0; i < number_of_particles; i++) {
        double pz_local = particle_list->get_positive_particle(event_id, i).pz;
        double E_local = particle_list->get_positive_particle(event_id, i).E;

        double rap_local;
        if (rap_type == 0) {
            double mass = (
                    particle_list->get_positive_particle(event_id, i).mass);
            double pmag = sqrt(E_local*E_local - mass*mass);
            rap_local = 0.5*log((pmag + pz_local)/(pmag - pz_local));
        } else {
            rap_local = 0.5*log((E_local + pz_local)/(E_local - pz_local));
        }

        if (rap_local > rap_min && rap_local < rap_max) {
            double px_local = (
                    particle_list->get_positive_particle(event_id, i).px);
            double py_local = (
                    particle_list->get_positive_particle(event_id, i).py);
            double p_perp = sqrt(px_local*px_local + py_local*py_local);
            if (p_perp > pT_min && p_perp < pT_max) {
                double p_phi = atan2(py_local, px_local);
                int p_idx = (int)((p_perp - pT_min)/dpT);
                for (int iorder = 0; iorder < order_max; iorder++) {
                    double cos_nphi = cos(iorder*p_phi);
                    double sin_nphi = sin(iorder*p_phi);
                    event_Qn_real[iorder] += cos_nphi;
                    event_Qn_imag[iorder] += sin_nphi;
                    event_Qn_real_err[iorder] += cos_nphi*cos_nphi;
                    event_Qn_imag_err[iorder] += sin_nphi*sin_nphi;
                    event_Qn_diff_real[iorder][p_idx] += cos_nphi;
                    event_Qn_diff_imag[iorder][p_idx] += sin_nphi;
                    event_Qn_diff_real_err[iorder][p_idx] += cos_nphi*cos_nphi;
                    event_Qn_diff_imag_err[iorder][p_idx] += sin_nphi*sin_nphi;
                }
            }
        }
    }
}


//! this function computes the pT-integrated and pT-differential Qn vector
//! of negative charged hadrons within a given rapidity region in one event
void singleParticleSpectra::calculate_Qn_vector_negative_charge(int event_id,
        double *event_Qn_real, double *event_Qn_real_err,
        double *event_Qn_imag, double *event_Qn_imag_err,
        double **event_Qn_diff_real, double **event_Qn_diff_real_err,
        double **event_Qn_diff_imag, double **event_Qn_diff_imag_err) {
    
    // first clean the results arrays
    for (int i = 0; i < order_max; i++) {
        event_Qn_real[i] = 0.0;
        event_Qn_real_err[i] = 0.0;
        event_Qn_imag[i] = 0.0;
        event_Qn_imag_err[i] = 0.0;
        for (int j = 0; j < npT; j++) {
            event_Qn_diff_real[i][j] = 0.0;
            event_Qn_diff_real_err[i][j] = 0.0;
            event_Qn_diff_imag[i][j] = 0.0;
            event_Qn_diff_imag_err[i][j] = 0.0;
        }
    }

    int number_of_particles = (
                    particle_list->get_number_of_negative_particles(event_id));
    for (int i = 0; i < number_of_particles; i++) {
        double pz_local = particle_list->get_negative_particle(event_id, i).pz;
        double E_local = particle_list->get_negative_particle(event_id, i).E;

        double rap_local;
        if (rap_type == 0) {
            double mass = (
                    particle_list->get_negative_particle(event_id, i).mass);
            double pmag = sqrt(E_local*E_local - mass*mass);
            rap_local = 0.5*log((pmag + pz_local)/(pmag - pz_local));
        } else {
            rap_local = 0.5*log((E_local + pz_local)/(E_local - pz_local));
        }

        if (rap_local > rap_min && rap_local < rap_max) {
            double px_local = (
                    particle_list->get_negative_particle(event_id, i).px);
            double py_local = (
                    particle_list->get_negative_particle(event_id, i).py);
            double p_perp = sqrt(px_local*px_local + py_local*py_local);
            if (p_perp > pT_min && p_perp < pT_max) {
                double p_phi = atan2(py_local, px_local);
                int p_idx = (int)((p_perp - pT_min)/dpT);
                for (int iorder = 0; iorder < order_max; iorder++) {
                    double cos_nphi = cos(iorder*p_phi);
                    double sin_nphi = sin(iorder*p_phi);
                    event_Qn_real[iorder] += cos_nphi;
                    event_Qn_imag[iorder] += sin_nphi;
                    event_Qn_real_err[iorder] += cos_nphi*cos_nphi;
                    event_Qn_imag_err[iorder] += sin_nphi*sin_nphi;
                    event_Qn_diff_real[iorder][p_idx] += cos_nphi;
                    event_Qn_diff_imag[iorder][p_idx] += sin_nphi;
                    event_Qn_diff_real_err[iorder][p_idx] += cos_nphi*cos_nphi;
                    event_Qn_diff_imag_err[iorder][p_idx] += sin_nphi*sin_nphi;
                }
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
    if (rap_type == 0) {
        filename << path << "/particle_" << particle_monval << "_vndata"
                 << "_eta_" << rap_min << "_" << rap_max << ".dat";
    } else {
        filename << path << "/particle_" << particle_monval << "_vndata"
                 << "_y_" << rap_min << "_" << rap_max << ".dat";
    }
    ofstream output(filename.str().c_str());

    double total_N = Qn_vector_real[0];
    double dN_ev_avg = Qn_vector_real[0]/total_number_of_events/drapidity;
    double dN_ev_avg_err = sqrt(dN_ev_avg/total_number_of_events)/drapidity;
    if (particle_monval == 333) {
        // for phi(1020) need to rescale the yield by
        // reconstruction branching ratio
        total_N = total_N/reconst_branching_ratio;
        dN_ev_avg = dN_ev_avg/reconst_branching_ratio;
        dN_ev_avg_err = dN_ev_avg_err/reconst_branching_ratio;
    }
    output << scientific << setw(18) << setprecision(8) 
           << 0 << "   " << dN_ev_avg << "   " << dN_ev_avg_err << "   " 
           << 0.0 << "   " << 0.0 << endl;
    for (int iorder = 1; iorder < order_max; iorder++) {
        double vn_evavg_real = 0.0;
        double vn_evavg_imag = 0.0;
        double vn_real_err = 0.0;
        double vn_imag_err = 0.0;
        if (Qn_vector_real[0] > 0) {
            vn_evavg_real = Qn_vector_real[iorder]/Qn_vector_real[0];
            vn_evavg_imag = Qn_vector_imag[iorder]/Qn_vector_real[0];
            vn_real_err = (sqrt(Qn_vector_real_err[iorder]/Qn_vector_real[0] 
                                - vn_evavg_real*vn_evavg_real)
                           /sqrt(Qn_vector_real[0]));
            vn_imag_err = (sqrt(Qn_vector_imag_err[iorder]/Qn_vector_real[0] 
                                - vn_evavg_imag*vn_evavg_imag)
                           /sqrt(Qn_vector_real[0]));
            if (isnan(vn_real_err)) {
                vn_real_err = 0.0;
            }
            if (isnan(vn_imag_err)) {
                vn_imag_err = 0.0;
            }
        }
        output << scientific << setw(18) << setprecision(8) << iorder << "   " 
               << vn_evavg_real << "   " << vn_real_err << "   " 
               << vn_evavg_imag << "   " << vn_imag_err << endl;
    }
    // output total number of particles in the last row
    // this quantities is useful when one wants to reconst the Qn vectors
    output << scientific << setw(18) << setprecision(8) << 99 << "   "
           << total_N << "   " << 0.0 << "   " << 0.0 << "   " << 0.0 << endl;
    output.close();
    
    // pT-differential flow
    ostringstream filename_diff;
    if (rap_type == 0) {
        filename_diff << path 
                      << "/particle_" << particle_monval << "_vndata_diff"
                      << "_eta_" << rap_min << "_" << rap_max << ".dat";
    } else {
        filename_diff << path 
                      << "/particle_" << particle_monval << "_vndata_diff"
                      << "_y_" << rap_min << "_" << rap_max << ".dat";
    }
    ofstream output_diff(filename_diff.str().c_str());

    for (int ipT = 0; ipT < npT - 1; ipT++) {
        double total_NpT = Qn_diff_vector_real[0][ipT];
        double dNpT_ev_avg = Qn_diff_vector_real[0][ipT]/total_number_of_events;
        double dNpT_ev_avg_err = sqrt(dNpT_ev_avg/total_number_of_events);
        if (particle_monval == 333) {
            // for phi(1020) need to rescale the yield by
            // reconstruction branching ratio
            total_NpT = total_NpT/reconst_branching_ratio;
            dNpT_ev_avg = dNpT_ev_avg/reconst_branching_ratio;
            dNpT_ev_avg_err = dNpT_ev_avg_err/reconst_branching_ratio;
        }
        double mean_pT, mean_pT_err;
        if (dNpT_ev_avg > 0.) {
            mean_pT = pT_mean_array[ipT]/Qn_diff_vector_real[0][ipT];
            mean_pT_err = ((pT_mean_array_err[ipT]/Qn_diff_vector_real[0][ipT] 
                           - mean_pT*mean_pT)
                          /sqrt(Qn_diff_vector_real[0][ipT]));
        } else {
            mean_pT = pT_array[ipT];
            mean_pT_err = 0.0;
        }
        output_diff << scientific << setw(18) << setprecision(8) 
                    << mean_pT << "   " << mean_pT_err << "   " 
                    << dNpT_ev_avg/mean_pT/dpT/(2*M_PI)/drapidity << "   " 
                    << dNpT_ev_avg_err/mean_pT/dpT/(2*M_PI)/drapidity;
        for (int iorder = 1; iorder < order_max; iorder++) {
            if (dNpT_ev_avg > 0.) {
                double vn_evavg_real = (Qn_diff_vector_real[iorder][ipT]
                                        /Qn_diff_vector_real[0][ipT]);
                double vn_evavg_imag = (Qn_diff_vector_imag[iorder][ipT]
                                        /Qn_diff_vector_real[0][ipT]);
                double vn_evavg_real_err = (
                        sqrt(Qn_diff_vector_real_err[iorder][ipT]
                                /Qn_diff_vector_real[0][ipT] 
                             - vn_evavg_real*vn_evavg_real)
                        /sqrt(Qn_diff_vector_real[0][ipT]));
                double vn_evavg_imag_err = (
                        sqrt(Qn_diff_vector_imag_err[iorder][ipT]
                                /Qn_diff_vector_real[0][ipT] 
                             - vn_evavg_imag*vn_evavg_imag)
                        /sqrt(Qn_diff_vector_real[0][ipT]));
                if (isnan(vn_evavg_real_err)) {
                    vn_evavg_real_err = 0.0;
                }
                if (isnan(vn_evavg_imag_err)) {
                    vn_evavg_imag_err = 0.0;
                }
                output_diff << scientific << setw(18) << setprecision(8) 
                            << vn_evavg_real << "   " 
                            << vn_evavg_real_err << "   " 
                            << vn_evavg_imag << "   "
                            << vn_evavg_imag_err << "   ";
            } else {
                output_diff << scientific << setw(18) << setprecision(8) 
                            << 0.0e0 << "   " << 0.0e0 << "   " 
                            << 0.0e0 << "   " << 0.0e0 << "   ";
            }
        }
        // output total number of particles in the last column
        // this quantities is useful when one wants to reconst the Qn vectors
        output_diff << scientific << setw(18) << setprecision(8) << total_NpT;
        output_diff << endl;
    }
    output_diff.close();
}

//! This function computes the 2-particle correlation for Qn vectors
//! within one event
//!     Real(Qn*conj(Qn)) for n = 0, 1, ... , order_max
//!     Real(Qn(pT)*conj(Qn)) for n = 0, 1, ... , order_max
//! self correlation is subtracted assuming full overlap
void singleParticleSpectra::calculate_two_particle_correlation(
        double *event_Qn_real, double *event_Qn_imag,
        double **event_Qn_diff_real, double **event_Qn_diff_imag) {
    for (int i = 0; i < order_max; i++) {
        double Q2_local = (event_Qn_real[i]*event_Qn_real[i]
                            + event_Qn_imag[i]*event_Qn_imag[i]
                            - event_Qn_real[0]);
        Qn2_vector[i] += Q2_local;
        Qn2_vector_err[i] += Q2_local*Q2_local;
        for (int j = 0; j < npT; j++) {
            double QnSP_pT_local = (
                    event_Qn_diff_real[i][j]*event_Qn_real[i]
                    + event_Qn_diff_imag[i][j]*event_Qn_imag[i]
                    - event_Qn_diff_real[0][j]);
            QnSP_diff_vector[i][j] += QnSP_pT_local;
            QnSP_diff_vector_err[i][j] += QnSP_pT_local*QnSP_pT_local;
        }
    }
}


//! This function computes the 2-particle correlation for Qn vectors
//! as a function of \delta \eta within one event
//!     Real(Qn(eta_1)*conj(Qn(eta_2))) for n = 0, 1, ... , order_max
//! self correlation is subtracted when eta_1 = eta_2
void singleParticleSpectra::calculate_two_particle_correlation_deltaeta(
        double **event_Qn_diff_real, double **event_Qn_diff_imag) {
    double *temp_corr = new double[N_rap];
    for (int i = 0; i < order_max; i++) {
        for (int j = 0; j < N_rap; j++) {
            temp_corr[j] = 0.0;
        }
        for (int ii = 0; ii < N_rap; ii++) {
            for (int jj = 0; jj < N_rap; jj++) {
                int rap_idx = ((ii - jj) - static_cast<int>(
                                                    rapidity_dis_min/drap));
                if (rap_idx >= 0 && rap_idx < N_rap) {
                    double QnSP_rap_local = (
                        event_Qn_diff_real[ii][i]*event_Qn_diff_real[jj][i]
                        + event_Qn_diff_imag[ii][i]*event_Qn_diff_imag[jj][i]);
                    if (ii == jj) {
                        QnSP_rap_local -= event_Qn_diff_real[ii][0];
                    }
                    temp_corr[rap_idx] += QnSP_rap_local;
                }
            }
        }
        for (int j = 0; j < N_rap; j++) {
            QnSP_eta12[i][j] += temp_corr[j];
            QnSP_eta12_err[i][j] += temp_corr[j]*temp_corr[j];
        }
    }
    delete[] temp_corr;
}

//! This function outputs the event averaged two-particle flow correlation
void singleParticleSpectra::output_two_particle_correlation() {
    ostringstream filename;
    if (rap_type == 0) {
        filename << path << "/particle_" << particle_monval << "_vn2"
                 << "_eta_" << rap_min << "_" << rap_max << ".dat";
    } else {
        filename << path << "/particle_" << particle_monval << "_vn2"
                 << "_y_" << rap_min << "_" << rap_max << ".dat";
    }
    ofstream output(filename.str().c_str());
    output << "# n  vn{2}  vn{2}_err  <Qn*conj(Qn)>  <(Qn*conj(Qn))^2>"
           << endl;

    double num_pair = Qn2_vector[0]/total_number_of_events;
    double num_pair_stdsq = (
            Qn2_vector_err[0]/total_number_of_events - num_pair*num_pair);
    double num_pair_err = 0.0;
    if (num_pair_stdsq > 0) {
        num_pair_err = sqrt(num_pair_stdsq/total_number_of_events);
    }
    output << scientific << setw(18) << setprecision(8)
           << 0 << "  " << num_pair << "  " << num_pair_err << "  "
           << num_pair << "  " << Qn2_vector_err[0]/total_number_of_events
           << endl;
    for (int i = 1; i < order_max; i++) {
        double vn2_avg = 0.0;
        double vn2_err = 0.0;
        double Qn2_avg = Qn2_vector[i]/total_number_of_events;
        if (Qn2_avg > 0.) {
            vn2_avg = sqrt(Qn2_avg/num_pair);
            double Qn2_stdsq = (
                Qn2_vector_err[i]/total_number_of_events - Qn2_avg*Qn2_avg);
            if (Qn2_stdsq > 0) {
                double Qn2_err = sqrt(Qn2_stdsq/total_number_of_events);
                vn2_err = Qn2_err/num_pair;
            }
        }
        output << scientific << setw(18) << setprecision(8)
               << i << "  " << vn2_avg << "  " << vn2_err << "  "
               << Qn2_avg << "  " << Qn2_vector_err[i]/total_number_of_events
               << endl;
    }
}

//! This function outputs the two-particle flow correlation as a function of
//! delta eta between the two particles
void singleParticleSpectra::output_two_particle_correlation_rap() {
    ostringstream filename;
    filename << path << "/particle_" << particle_monval << "_vn2_eta12"
             << "_pT_" << vn_rapidity_dis_pT_min << "_"
             << vn_rapidity_dis_pT_max << ".dat";
    ofstream output(filename.str().c_str());
    if (rap_type == 0) {
        output << "# eta  vn{2}  vn{2}_err  <Qn*conj(Qn)>  <(Qn*conj(Qn))^2>"
               << endl;
    } else {
        output << "# y  vn{2}  vn{2}_err  <Qn*conj(Qn)>  <(Qn*conj(Qn))^2>"
               << endl;
    }

    for (int i = 0; i < N_rap; i++) {
        double eta_local = rapidity_dis_min + i*drap;
        double num_pair = QnSP_eta12[0][i]/total_number_of_events;
        double num_pair_stdsq = (
            QnSP_eta12_err[0][i]/total_number_of_events - num_pair*num_pair);
        double num_pair_err = 0.0;
        if (num_pair_stdsq > 0) {
            num_pair_err = sqrt(num_pair_stdsq/total_number_of_events);
        }
        output << scientific << setw(18) << setprecision(8)
               << eta_local << "  " << num_pair << "  " << num_pair_err
               << "  " << num_pair << "  "
               << Qn2_vector_err[0]/total_number_of_events << "  ";
        for (int iorder = 1; iorder < order_max; iorder++) {
            double vn2_avg = 0.0;
            double vn2_err = 0.0;
            double Qn2_avg = QnSP_eta12[iorder][i]/total_number_of_events;
            if (Qn2_avg > 0.) {
                vn2_avg = sqrt(Qn2_avg/num_pair);
                double Qn2_stdsq = (
                            QnSP_eta12_err[iorder][i]/total_number_of_events
                            - Qn2_avg*Qn2_avg);
                if (Qn2_stdsq > 0) {
                    double Qn2_err = sqrt(Qn2_stdsq/total_number_of_events);
                    vn2_err = Qn2_err/num_pair;
                }
            }
            output << scientific << setw(18) << setprecision(8)
                   << vn2_avg << "  " << vn2_err << "  "
                   << Qn2_avg << "  " 
                   << QnSP_eta12_err[iorder][i]/total_number_of_events << "  ";
        }
        output << endl;
    }
    output.close();
}

//! This function computes the 3-particle correlation for Qn vectors
//! within one event
//!     C_nmk = Real(Q_n*Q_m*conj(Q_k)) for (112), (123), (224), (235)
//! self correlation is subtracted assuming Qk's sample >= Qn's and Qm's
//! flag = 0: Qn = Qm <= Qk, flag = 1: Qn != Qm, Qn \in Qk, Qm \in Qk
//! flag = 2: no overlap
void singleParticleSpectra::calculate_three_particle_correlation(
        double *event_Q1_real, double *event_Q1_imag,
        double *event_Q2_real, double *event_Q2_imag,
        double *event_Q3_real, double *event_Q3_imag, int flag,
        double *corr, double *corr_err) {
    // C_nmk[0] = C_000 = N(N-1)(N-2) is the number of pairs
    // C_nmk[1] = C_112, C_nmk[2] = C_123, C_nmk[3] = C_224, C_nmk[4] = C_235
    // C_nmk[5] = C_134, C_nmk[6] = C_246, C_nmk[7] = C_336, C_nmk[8] = C_347
    int n[9] = {0, 1, 1, 2, 2, 1, 2, 3, 3};
    int m[9] = {0, 1, 2, 2, 3, 3, 4, 3, 4};
    for (int i = 0; i < num_corr; i++) {
        int k = n[i] + m[i];
        if (k > order_max) {
            cout << "Error: the vn needed for three-particle correlation "
                 << "is not computed, order_max = " << order_max
                 << "n+m = " << k << ". Make sure order_max > n+m!" << endl;
            exit(1);
        }
        double Qn_Qm_Qkstar = (
            (event_Q1_real[n[i]]*event_Q2_real[m[i]]
             - event_Q1_imag[n[i]]*event_Q2_imag[m[i]])*event_Q3_real[k]
            + (event_Q1_real[n[i]]*event_Q2_imag[m[i]]
               + event_Q1_imag[n[i]]*event_Q2_real[m[i]])*event_Q3_imag[k]);
        // full overlap between Q2 and Q3 = Q1[n]*conj(Q2[n])
        double Qn_Qnstar = (event_Q1_real[n[i]]*event_Q2_real[n[i]]
                            + event_Q1_imag[n[i]]*event_Q2_imag[n[i]]);
        // full overlap between Q1 and Q3 = Q2[m]*conj(Q1[m])
        double Qm_Qmstar = (event_Q2_real[m[i]]*event_Q1_real[m[i]]
                            + event_Q2_imag[m[i]]*event_Q1_imag[m[i]]);
        // full overlap between Q1 and Q2 = Q1[k]*conj(Q3[k])
        double Qk_Qkstar = (event_Q1_real[k]*event_Q3_real[k]
                            + event_Q1_imag[k]*event_Q3_imag[k]);
        double corr_local = 0.;
        if (flag == 0) {
            corr_local = (Qn_Qm_Qkstar - Qn_Qnstar - Qm_Qmstar - Qk_Qkstar
                             + 2.*event_Q1_real[0]);
        } else if (flag == 1) {
            corr_local = (Qn_Qm_Qkstar - Qn_Qnstar - Qm_Qmstar);
        } else if (flag == 2) {
            corr_local = Qn_Qm_Qkstar;
        }
        corr[i] += corr_local;
        corr_err[i] += corr_local*corr_local;
    }
}

//! This function computes the 3-particle correlation as a function of relative
//! rapidity between particle 1 and particle 2 within one event
//!     C_nmk(eta_12) = Real(Q_n(eta1)*Q_m(eta2)*conj(Q_k))
//! self correlation is subtracted assuming Qk's sample >= Qn's and Qm's
//! Qn != Qm, Qn \in Qk, Qm \in Qk
//! flag == 1 : eta_12, flag == 2: eta_13
//! flag_ch == 1: same charge, flag_ch == 2: opposite charge
void singleParticleSpectra::calculate_three_particle_correlation_deltaeta(
        double **event_Q1_real, double **event_Q1_imag,
        double **event_Q2_real, double **event_Q2_imag,
        double **event_Q3_real, double **event_Q3_imag, int flag, int flag_ch,
        double **corr_rap, double **corr_rap_err) {
    // C_nmk[0] = C_000 = N(N-1)(N-2) is the number of pairs
    // C_nmk[1] = C_112, C_nmk[2] = C_123, C_nmk[3] = C_224, C_nmk[4] = C_235
    // C_nmk[5] = C_134, C_nmk[6] = C_246, C_nmk[7] = C_336, C_nmk[8] = C_347
    int n[9] = {0, 1, 1, 2, 2, 1, 2, 3, 3};
    int m[9] = {0, 1, 2, 2, 3, 3, 4, 3, 4};
    double *temp_corr = new double[N_rap];
    for (int i = 0; i < num_corr; i++) {
        int k = n[i] + m[i];
        if (k > order_max) {
            cout << "Error: the vn needed for three-particle correlation "
                 << "is not computed, order_max = " << order_max
                 << "n+m = " << k << ". Make sure order_max > n+m!" << endl;
            exit(1);
        }

        double Q2_real = 0.0;
        double Q2_imag = 0.0;
        double Q3_real = 0.0;
        double Q3_imag = 0.0;
        if (flag == 1) {
            for (int j = 0; j < N_rap; j++) {
                Q3_real += event_Q3_real[j][k];
                Q3_imag += event_Q3_imag[j][k];
            }
        } else if (flag == 2) {
            for (int j = 0; j < N_rap; j++) {
                Q2_real += event_Q2_real[j][k];
                Q2_imag += event_Q2_imag[j][k];
            }
        }

        for (int kk = 0; kk < N_rap; kk++) {
            temp_corr[kk] = 0.0;
        }
        for (int ii = 0; ii < N_rap; ii++) {
            for (int jj = 0; jj < N_rap; jj++) {
                int rap_idx = ((ii - jj) - static_cast<int>(
                                                    rapidity_dis_min/drap));
                if (rap_idx >= 0 && rap_idx < N_rap) {
                    double Qn_Qm_Qkstar = 0.0;
                    double Qn_Qnstar = 0.0;
                    double Qm_Qmstar = 0.0;
                    double Qk_Qkstar = 0.0;
                    if (flag == 1) {
                        Qn_Qm_Qkstar = (
                            (event_Q1_real[ii][n[i]]*event_Q2_real[jj][m[i]]
                             - event_Q1_imag[ii][n[i]]*event_Q2_imag[jj][m[i]])
                            *Q3_real
                            + (event_Q1_real[ii][n[i]]*event_Q2_imag[jj][m[i]]
                               + event_Q1_imag[ii][n[i]]*event_Q2_real[jj][m[i]])
                            *Q3_imag);
                        // full overlap between Q2 and Q3 = Q1[n]*conj(Q2[n])
                        Qn_Qnstar = (
                            event_Q1_real[ii][n[i]]*event_Q2_real[jj][n[i]]
                            + event_Q1_imag[ii][n[i]]*event_Q2_imag[jj][n[i]]);
                        // full overlap between Q1 and Q3 = Q2[m]*conj(Q1[m])
                        Qm_Qmstar = (
                            event_Q2_real[ii][m[i]]*event_Q1_real[jj][m[i]]
                            + event_Q2_imag[ii][m[i]]*event_Q1_imag[jj][m[i]]);
                        // full overlap between Q1 and Q2 = Q1[k]*conj(Q3[k])
                        if (ii == jj && flag_ch == 1) {
                            Qk_Qkstar = (event_Q1_real[ii][k]*Q3_real
                                         + event_Q1_imag[ii][k]*Q3_imag);
                        }
                        double corr_local;
                        if (fabs(Qk_Qkstar) < 1e-15) {
                            corr_local = Qn_Qm_Qkstar - Qn_Qnstar - Qm_Qmstar;
                        } else {
                            corr_local = (
                                Qn_Qm_Qkstar - Qn_Qnstar - Qm_Qmstar
                                - Qk_Qkstar + 2.*event_Q1_real[ii][0]);
                        }
                        temp_corr[rap_idx] += corr_local;
                    } else if (flag == 2) {
                        Qn_Qm_Qkstar = (
                            (event_Q1_real[ii][n[i]]*Q2_real
                             - event_Q1_imag[ii][n[i]]*Q2_imag)
                            *event_Q3_real[jj][k]
                            + (event_Q1_real[ii][n[i]]*Q2_imag
                               + event_Q1_imag[ii][n[i]]*Q2_real)
                            *event_Q3_imag[jj][k]);
                        // full overlap between Q2 and Q3 = Q1[n]*conj(Q2[n])
                        Qn_Qnstar = (
                            event_Q1_real[ii][n[i]]*event_Q2_real[jj][n[i]]
                            + event_Q1_imag[ii][n[i]]*event_Q2_imag[jj][n[i]]);
                        // full overlap between Q1 and Q3 = Q2[m]*conj(Q1[m])
                        if (ii == jj) {
                            Qm_Qmstar = (Q2_real*event_Q1_real[jj][m[i]]
                                         + Q2_imag*event_Q1_imag[jj][m[i]]);
                        }
                        // full overlap between Q1 and Q2 = Q1[k]*conj(Q3[k])
                        if (flag_ch == 1) {
                            Qk_Qkstar = (
                                event_Q1_real[ii][k]*event_Q3_real[jj][k]
                                + event_Q1_imag[ii][k]*event_Q3_imag[jj][k]);
                        }

                        double corr_local;
                        if (fabs(Qk_Qkstar) > 1e-15
                                && fabs(Qm_Qmstar) > 1e-15) {
                            corr_local = (
                                Qn_Qm_Qkstar - Qn_Qnstar - Qm_Qmstar
                                - Qk_Qkstar + 2.*event_Q1_real[ii][0]);
                        } else {
                            corr_local = (Qn_Qm_Qkstar - Qn_Qnstar
                                          - Qm_Qmstar - Qk_Qkstar);
                        }
                        temp_corr[rap_idx] += corr_local;
                    }
                }
            }
        }
        for (int kk = 0; kk < N_rap; kk++) {
            corr_rap[i][kk] += temp_corr[kk];
            corr_rap_err[i][kk] += temp_corr[kk]*temp_corr[kk];
        }
    }
    delete[] temp_corr;
}

//! This function computes the 4-particle correlation for symmetric cumulants
//! using Qn vectors within one event
//!     SC_mn = <Q1_m*conj(Q2_m)*Q3_n*conj(Q4_n)>
//              - <Q1_m*conj(Q2_m)><Q3_n*conj(Q4_n)>
//! for (32), (42), (52), (43), (53)
//! self correlation is subtracted assuming Qk's sample >= Qn's and Qm's
void singleParticleSpectra::calculate_four_particle_correlation_SC(
        double *event_Q1_real, double *event_Q1_imag,
        double *event_Q2_real, double *event_Q2_imag,
        double *event_Q3_real, double *event_Q3_imag,
        double *event_Q4_real, double *event_Q4_imag, int flag,
        double *corr, double *corr_err) {
    // SC_mn[0] = SC_00 = N(N-1)(N-2)(N-3) is the number of pairs
    // SC_mn[1] = SC_32, SC_mn[2] = SC_42
    // SC_mn[3] = SC_52, SC_mn[4] = SC_43, SC_mn[5] = SC_53
    int m[6] = {0, 3, 4, 5, 4, 5};
    int n[6] = {0, 2, 2, 2, 3, 3};
    for (int i = 0; i < SC_num_corr; i++) {
        int k = m[i] + n[i];
        if (k > order_max) {
            cout << "Error: the vn needed for symmetric cumulants is not "
                 << "computed, order_max = " << order_max
                 << "m+n = " << k << ". Make sure order_max > m+n!" << endl;
            exit(1);
        }
        int l = m[i] - n[i];
        double Q1_conjQ2 = (event_Q1_real[m[i]]*event_Q2_real[m[i]]
                            + event_Q1_imag[m[i]]*event_Q2_imag[m[i]]);
        double Q3_conjQ4 = (event_Q3_real[n[i]]*event_Q4_real[n[i]]
                            + event_Q3_imag[n[i]]*event_Q4_imag[n[i]]);
        // full overlap of Q1 and Q3
        double Q1Q3_overlap = (
            event_Q1_real[k]*event_Q2_real[m[i]]*event_Q4_real[n[i]]
            - event_Q1_real[k]*event_Q2_imag[m[i]]*event_Q4_imag[n[i]]
            + event_Q1_imag[k]*event_Q2_real[m[i]]*event_Q4_imag[n[i]]
            + event_Q1_imag[k]*event_Q2_imag[m[i]]*event_Q4_real[n[i]]
        );
        // full overlap of Q1 and Q4
        double Q1Q4_overlap = (
            event_Q1_real[l]*event_Q2_real[m[i]]*event_Q3_real[n[i]]
            + event_Q1_real[l]*event_Q2_imag[m[i]]*event_Q3_imag[n[i]]
            - event_Q1_imag[l]*event_Q2_real[m[i]]*event_Q3_imag[n[i]]
            + event_Q1_imag[l]*event_Q2_imag[m[i]]*event_Q3_real[n[i]]
        );
        // full overlap of Q2 and Q4
        double Q2Q4_overlap = (
            event_Q2_real[k]*event_Q1_real[m[i]]*event_Q3_real[n[i]]
            - event_Q2_real[k]*event_Q1_imag[m[i]]*event_Q3_imag[n[i]]
            + event_Q2_imag[k]*event_Q1_real[m[i]]*event_Q3_imag[n[i]]
            + event_Q2_imag[k]*event_Q1_imag[m[i]]*event_Q3_real[n[i]]
        );
        // full overlap of Q2 and Q3
        double Q2Q3_overlap = (
            event_Q2_real[l]*event_Q1_real[m[i]]*event_Q4_real[n[i]]
            + event_Q2_real[l]*event_Q1_imag[m[i]]*event_Q4_imag[n[i]]
            - event_Q2_imag[l]*event_Q1_real[m[i]]*event_Q4_imag[n[i]]
            + event_Q2_imag[l]*event_Q1_imag[m[i]]*event_Q4_real[n[i]]
        );
        // full overlap of Q1 with Q3 and Q2 with Q4
        double Q1Q3_Q2Q4_overlap = (event_Q1_real[k]*event_Q2_real[k]
                                    + event_Q1_imag[k]*event_Q2_imag[k]);
        // full overlap of Q1 with Q4 and Q2 with Q3
        double Q1Q4_Q2Q3_overlap = (event_Q1_real[l]*event_Q2_real[l]
                                    + event_Q1_imag[l]*event_Q2_imag[l]);
        // full overlap of Q1 with Q2 and Q3
        double Q1Q2Q3_overlap = (event_Q1_real[n[i]]*event_Q4_real[n[i]]
                                 + event_Q1_imag[n[i]]*event_Q4_imag[n[i]]);
        // full overlap of Q1 with Q2 and Q4
        double Q1Q2Q4_overlap = (event_Q1_real[n[i]]*event_Q3_real[n[i]]
                                 + event_Q1_imag[n[i]]*event_Q3_imag[n[i]]);
        // full overlap of Q1 with Q3 and Q4
        double Q1Q3Q4_overlap = (event_Q1_real[m[i]]*event_Q2_real[m[i]]
                                 + event_Q1_imag[m[i]]*event_Q2_imag[m[i]]);
        // full overlap of Q2 with Q3 and Q4
        double Q2Q3Q4_overlap = (event_Q1_real[m[i]]*event_Q2_real[m[i]]
                                 + event_Q1_imag[m[i]]*event_Q2_imag[m[i]]);
        // full overlap of Q1 with Q2, Q3, and Q4
        double Q1Q2Q3Q4_overlap = event_Q1_real[0];

        double corr_local = 0.;
        if (flag == 0) {
            corr_local = (
                Q1_conjQ2*Q3_conjQ4 - Q1Q3_overlap - Q1Q4_overlap
                - Q2Q3_overlap - Q2Q4_overlap
                + Q1Q3_Q2Q4_overlap + Q1Q4_Q2Q3_overlap
                - (event_Q1_real[0] - 4.)/2.*Q1Q2Q3_overlap
                - (event_Q1_real[0] - 4.)/2.*Q1Q2Q4_overlap
                - (event_Q1_real[0] - 4.)/2.*Q1Q3Q4_overlap
                - (event_Q1_real[0] - 4.)/2.*Q2Q3Q4_overlap
                + Q1Q2Q3Q4_overlap*(Q1Q2Q3Q4_overlap - 6.)
            );
        }
        corr[i] += corr_local;
        corr_err[i] += corr_local*corr_local;
    }
}

//! This function outputs the event averaged three-particle correlation
void singleParticleSpectra::output_three_particle_correlation() {
    ostringstream filename;
    if (rap_type == 0) {
        filename << path << "/particle_" << particle_monval << "_Cmnk"
                 << "_eta_" << rap_min << "_" << rap_max << ".dat";
    } else {
        filename << path << "/particle_" << particle_monval << "_Cmnk"
                 << "_y_" << rap_min << "_" << rap_max << ".dat";
    }
    ofstream output(filename.str().c_str());
    output << "# n  C_nmk  C_nmk_err" << endl;

    double num_pair = C_nmk[0]/total_number_of_events;
    double num_pair_stdsq = (
            C_nmk_err[0]/total_number_of_events - num_pair*num_pair);
    double num_pair_err = 0.0;
    if (num_pair_stdsq > 0) {
        num_pair_err = sqrt(num_pair_stdsq/total_number_of_events);
    }
    output << scientific << setw(18) << setprecision(8)
           << 0 << "  " << num_pair << "  " << num_pair_err << endl;
    for (int i = 1; i < num_corr; i++) {
        double Cnmk_avg = C_nmk[i]/total_number_of_events;
        double Cnmk_stdsq = (
                C_nmk_err[i]/total_number_of_events - Cnmk_avg*Cnmk_avg);
        Cnmk_avg = Cnmk_avg/num_pair;
        double Cnmk_err = 0.0;
        if (Cnmk_stdsq > 0) {
            Cnmk_err = sqrt(Cnmk_stdsq/total_number_of_events);
            Cnmk_err = Cnmk_err/num_pair;
        }
        output << scientific << setw(18) << setprecision(8)
               << i << "  " << Cnmk_avg << "  " << Cnmk_err << endl;
    }
    if (flag_charge_dependence == 1) {
        ostringstream filename_ss;
        if (rap_type == 0) {
            filename_ss << path << "/particle_" << particle_monval
                        << "_Cmnk_ss_eta_" << rap_min << "_"
                        << rap_max << ".dat";
        } else {
            filename_ss << path << "/particle_" << particle_monval
                        << "_Cmnk_ss_y_" << rap_min << "_"
                        << rap_max << ".dat";
        }
        ofstream output_ss(filename_ss.str().c_str());
        output_ss << "# n  C_nmk_ss  C_nmk_ss_err" << endl;

        double num_pair_ss = C_nmk_ss[0]/total_number_of_events;
        double num_pair_ss_stdsq = (
                C_nmk_ss_err[0]/total_number_of_events
                - num_pair_ss*num_pair_ss);
        double num_pair_ss_err = 0.0;
        if (num_pair_ss_stdsq > 0) {
            num_pair_ss_err = sqrt(num_pair_ss_stdsq/total_number_of_events);
        }
        output_ss << scientific << setw(18) << setprecision(8)
                  << 0 << "  " << num_pair_ss << "  " << num_pair_ss_err
                  << endl;
        for (int i = 1; i < num_corr; i++) {
            double Cnmk_ss_avg = C_nmk_ss[i]/total_number_of_events;
            double Cnmk_ss_stdsq = (
                    C_nmk_ss_err[i]/total_number_of_events
                    - Cnmk_ss_avg*Cnmk_ss_avg);
            Cnmk_ss_avg = Cnmk_ss_avg/num_pair_ss;
            double Cnmk_ss_err = 0.0;
            if (Cnmk_ss_stdsq > 0) {
                Cnmk_ss_err = sqrt(Cnmk_ss_stdsq/total_number_of_events);
                Cnmk_ss_err = Cnmk_ss_err/num_pair;
            }
            output_ss << scientific << setw(18) << setprecision(8)
                      << i << "  " << Cnmk_ss_avg << "  " << Cnmk_ss_err
                      << endl;
        }
        ostringstream filename_os;
        if (rap_type == 0) {
            filename_os << path << "/particle_" << particle_monval
                        << "_Cmnk_os_eta_" << rap_min << "_"
                        << rap_max << ".dat";
        } else {
            filename_os << path << "/particle_" << particle_monval
                        << "_Cmnk_os_y_" << rap_min << "_"
                        << rap_max << ".dat";
        }
        ofstream output_os(filename_os.str().c_str());
        output_os << "# n  C_nmk_os  C_nmk_os_err" << endl;

        double num_pair_os = C_nmk_os[0]/total_number_of_events;
        double num_pair_os_stdsq = (
                C_nmk_os_err[0]/total_number_of_events
                - num_pair_os*num_pair_os);
        double num_pair_os_err = 0.0;
        if (num_pair_os_stdsq > 0) {
            num_pair_os_err = sqrt(num_pair_os_stdsq/total_number_of_events);
        }
        output_os << scientific << setw(18) << setprecision(8)
                  << 0 << "  " << num_pair_os << "  " << num_pair_os_err
                  << endl;
        for (int i = 1; i < num_corr; i++) {
            double Cnmk_os_avg = C_nmk_os[i]/total_number_of_events;
            double Cnmk_os_stdsq = (
                    C_nmk_os_err[i]/total_number_of_events
                    - Cnmk_os_avg*Cnmk_os_avg);
            Cnmk_os_avg = Cnmk_os_avg/num_pair_os;
            double Cnmk_os_err = 0.0;
            if (Cnmk_os_stdsq > 0) {
                Cnmk_os_err = sqrt(Cnmk_os_stdsq/total_number_of_events);
                Cnmk_os_err = Cnmk_os_err/num_pair;
            }
            output_os << scientific << setw(18) << setprecision(8)
                      << i << "  " << Cnmk_os_avg << "  " << Cnmk_os_err
                      << endl;
        }
    }
}

//! This function outputs the rapidity dependent three-particle correlation
void singleParticleSpectra::output_three_particle_correlation_rap() {
    ostringstream filename1, filename2;
    filename1 << path << "/particle_" << particle_monval << "_Cmnk_eta12"
              << "_pT_" << vn_rapidity_dis_pT_min << "_"
              << vn_rapidity_dis_pT_max << ".dat";
    filename2 << path << "/particle_" << particle_monval << "_Cmnk_eta13"
              << "_pT_" << vn_rapidity_dis_pT_min << "_"
              << vn_rapidity_dis_pT_max << ".dat";
    ofstream output1(filename1.str().c_str());
    ofstream output2(filename2.str().c_str());
    if (rap_type == 0) {
        output1 << "# eta12  C_nmk  C_nmk_err "
                << "(000, 112, 123, 224, 235, 134, 246, 336, 347)" << endl;
        output2 << "# eta13  C_nmk  C_nmk_err " 
                << "(000, 112, 123, 224, 235, 134, 246, 336, 347)" << endl;
    } else {
        output1 << "# y12  C_nmk  C_nmk_err "
                << "(000, 112, 123, 224, 235, 134, 246, 336, 347)" << endl;
        output2 << "# y13  C_nmk  C_nmk_err " 
                << "(000, 112, 123, 224, 235, 134, 246, 336, 347)" << endl;
    }

    for (int j = 0; j < N_rap; j++) {
        double eta_local = rapidity_dis_min + j*drap;
        double num_pair = C_nmk_eta12[0][j]/total_number_of_events;
        double num_pair_stdsq = (
            C_nmk_eta12_err[0][j]/total_number_of_events - num_pair*num_pair);
        double num_pair_err = 0.0;
        if (num_pair_stdsq > 0) {
            num_pair_err = sqrt(num_pair_stdsq/total_number_of_events);
        }
        output1 << scientific << setw(18) << setprecision(8)
                << eta_local << "  " << num_pair << "  " << num_pair_err
                << "  ";
        double num_pair_2 = C_nmk_eta13[0][j]/total_number_of_events;
        double num_pair_2_stdsq = (C_nmk_eta13_err[0][j]/total_number_of_events
                                   - num_pair_2*num_pair_2);
        num_pair_err = 0.0;
        if (num_pair_2_stdsq > 0) {
            num_pair_err = sqrt(num_pair_2_stdsq/total_number_of_events);
        }
        output2 << scientific << setw(18) << setprecision(8)
                << eta_local << "  " << num_pair_2 << "  " << num_pair_err
                << "  ";
        for (int i = 1; i < num_corr; i++) {
            double Cnmk_avg = C_nmk_eta12[i][j]/total_number_of_events;
            double Cnmk_stdsq = (C_nmk_eta12_err[i][j]/total_number_of_events
                                 - Cnmk_avg*Cnmk_avg);
            Cnmk_avg = Cnmk_avg/num_pair;
            double Cnmk_err = 0.0;
            if (Cnmk_stdsq > 0) {
                Cnmk_err = sqrt(Cnmk_stdsq/total_number_of_events);
                Cnmk_err = Cnmk_err/num_pair;
            }
            output1 << scientific << setw(18) << setprecision(8)
                    << Cnmk_avg << "  " << Cnmk_err << "  ";
            Cnmk_avg = C_nmk_eta13[i][j]/total_number_of_events;
            Cnmk_stdsq = (C_nmk_eta13_err[i][j]/total_number_of_events
                          - Cnmk_avg*Cnmk_avg);
            Cnmk_avg = Cnmk_avg/num_pair_2;
            Cnmk_err = 0.0;
            if (Cnmk_stdsq > 0) {
                Cnmk_err = sqrt(Cnmk_stdsq/total_number_of_events);
                Cnmk_err = Cnmk_err/num_pair_2;
            }
            Cnmk_stdsq = Cnmk_stdsq/num_pair;
            output2 << scientific << setw(18) << setprecision(8)
                    << Cnmk_avg << "  " << Cnmk_err << "  ";
        }
        output1 << endl;
        output2 << endl;
    }
    output1.close();
    output2.close();
    if (flag_charge_dependence == 1) {
        ostringstream filename_ss1, filename_ss2;
        ostringstream filename_os1, filename_os2;
        filename_ss1 << path << "/particle_" << particle_monval
                     << "_Cmnk_eta12_ss_pT_" << vn_rapidity_dis_pT_min
                     << "_" << vn_rapidity_dis_pT_max << ".dat";
        filename_os1 << path << "/particle_" << particle_monval
                     << "_Cmnk_eta12_os_pT_" << vn_rapidity_dis_pT_min
                     << "_" << vn_rapidity_dis_pT_max << ".dat";
        filename_ss2 << path << "/particle_" << particle_monval
                     << "_Cmnk_eta13_ss_pT_" << vn_rapidity_dis_pT_min
                     << "_" << vn_rapidity_dis_pT_max << ".dat";
        filename_os2 << path << "/particle_" << particle_monval
                     << "_Cmnk_eta13_os_pT_" << vn_rapidity_dis_pT_min
                     << "_" << vn_rapidity_dis_pT_max << ".dat";
        ofstream output_ss1(filename_ss1.str().c_str());
        ofstream output_os1(filename_os1.str().c_str());
        ofstream output_ss2(filename_ss2.str().c_str());
        ofstream output_os2(filename_os2.str().c_str());
        if (rap_type == 0) {
            output_ss1 << "# eta12  C_nmk_ss  C_nmk_ss_err" << endl;
            output_os1 << "# eta12  C_nmk_os  C_nmk_os_err" << endl;
            output_ss2 << "# eta13  C_nmk_ss  C_nmk_ss_err" << endl;
            output_os2 << "# eta13  C_nmk_os  C_nmk_os_err" << endl;
        } else {
            output_ss1 << "# y12  C_nmk_ss  C_nmk_ss_err" << endl;
            output_os1 << "# y12  C_nmk_os  C_nmk_os_err" << endl;
            output_ss2 << "# y13  C_nmk_ss  C_nmk_ss_err" << endl;
            output_os2 << "# y13  C_nmk_os  C_nmk_os_err" << endl;
        }

        for (int j = 0; j < N_rap; j++) {
            double eta_local = rapidity_dis_min + j*drap;
            double num_pair_ss = C_nmk_eta12_ss[0][j]/total_number_of_events;
            double num_pair_ss_stdsq = (
                    C_nmk_eta12_ss_err[0][j]/total_number_of_events
                    - num_pair_ss*num_pair_ss);
            double num_pair_ss_err = 0.0;
            if (num_pair_ss_stdsq > 0) {
                num_pair_ss_err = (
                        sqrt(num_pair_ss_stdsq/total_number_of_events));
            }
            output_ss1 << scientific << setw(18) << setprecision(8)
                       << eta_local << "  " << num_pair_ss << "  "
                       << num_pair_ss_err << "  ";
            double num_pair_os = C_nmk_eta12_os[0][j]/total_number_of_events;
            double num_pair_os_stdsq = (
                    C_nmk_eta12_os_err[0][j]/total_number_of_events
                    - num_pair_os*num_pair_os);
            double num_pair_os_err = 0.0;
            if (num_pair_os_stdsq > 0) {
                num_pair_os_err = (
                    sqrt(num_pair_os_stdsq/total_number_of_events));
            }
            output_os1 << scientific << setw(18) << setprecision(8)
                       << eta_local << "  " << num_pair_os << "  "
                       << num_pair_os_err << "  ";
            double num_pair_ss2 = C_nmk_eta13_ss[0][j]/total_number_of_events;
            double num_pair_ss2_stdsq = (
                    C_nmk_eta13_ss_err[0][j]/total_number_of_events
                    - num_pair_ss2*num_pair_ss2);
            double num_pair_ss2_err = 0.0;
            if (num_pair_ss2_stdsq > 0) {
                num_pair_ss2_err = (
                    sqrt(num_pair_ss2_stdsq/total_number_of_events));
            }
            output_ss2 << scientific << setw(18) << setprecision(8)
                       << eta_local << "  " << num_pair_ss2 << "  "
                       << num_pair_ss2_err << "  ";
            double num_pair_os2 = C_nmk_eta13_os[0][j]/total_number_of_events;
            double num_pair_os2_stdsq = (
                    C_nmk_eta13_os_err[0][j]/total_number_of_events
                    - num_pair_os2*num_pair_os2);
            double num_pair_os2_err = 0.0;
            if (num_pair_os2_stdsq > 0) {
                num_pair_os2_err = (
                    sqrt(num_pair_os2_stdsq/total_number_of_events));
            }
            output_os2 << scientific << setw(18) << setprecision(8)
                       << eta_local << "  " << num_pair_os2 << "  "
                       << num_pair_os2_err << "  ";
            for (int i = 1; i < num_corr; i++) {
                double Cnmk_ss_avg = (
                        C_nmk_eta12_ss[i][j]/total_number_of_events);
                double Cnmk_ss_stdsq = (
                        C_nmk_eta12_ss_err[i][j]/total_number_of_events
                        - Cnmk_ss_avg*Cnmk_ss_avg);
                Cnmk_ss_avg = Cnmk_ss_avg/num_pair_ss;
                double Cnmk_ss_err = 0.0;
                if (Cnmk_ss_stdsq > 0) {
                    Cnmk_ss_err = sqrt(Cnmk_ss_stdsq/total_number_of_events);
                    Cnmk_ss_err = Cnmk_ss_err/num_pair_ss;
                }
                output_ss1 << scientific << setw(18) << setprecision(8)
                           << Cnmk_ss_avg << "  " << Cnmk_ss_err << "  ";
                double Cnmk_os_avg = (
                        C_nmk_eta12_os[i][j]/total_number_of_events);
                double Cnmk_os_stdsq = (
                        C_nmk_eta12_os_err[i][j]/total_number_of_events
                        - Cnmk_os_avg*Cnmk_os_avg);
                Cnmk_os_avg = Cnmk_os_avg/num_pair_os;
                double Cnmk_os_err = 0.0;
                if (Cnmk_os_stdsq > 0) {
                    Cnmk_os_err = sqrt(Cnmk_os_stdsq/total_number_of_events);
                    Cnmk_os_err = Cnmk_os_err/num_pair_os;
                }
                output_os1 << scientific << setw(18) << setprecision(8)
                           << Cnmk_os_avg << "  " << Cnmk_os_err << "  ";
                Cnmk_ss_avg = (
                        C_nmk_eta13_ss[i][j]/total_number_of_events);
                Cnmk_ss_stdsq = (
                        C_nmk_eta13_ss_err[i][j]/total_number_of_events
                        - Cnmk_ss_avg*Cnmk_ss_avg);
                Cnmk_ss_avg = Cnmk_ss_avg/num_pair_ss2;
                Cnmk_ss_err = 0.0;
                if (Cnmk_ss_stdsq > 0) {
                    Cnmk_ss_err = sqrt(Cnmk_ss_stdsq/total_number_of_events);
                    Cnmk_ss_err = Cnmk_ss_err/num_pair_ss2;
                }
                output_ss2 << scientific << setw(18) << setprecision(8)
                           << Cnmk_ss_avg << "  " << Cnmk_ss_err << "  ";
                Cnmk_os_avg = (
                        C_nmk_eta13_os[i][j]/total_number_of_events);
                Cnmk_os_stdsq = (
                        C_nmk_eta13_os_err[i][j]/total_number_of_events
                        - Cnmk_os_avg*Cnmk_os_avg);
                Cnmk_os_avg = Cnmk_os_avg/num_pair_os2;
                Cnmk_os_err = 0.0;
                if (Cnmk_os_stdsq > 0) {
                    Cnmk_os_err = sqrt(Cnmk_os_stdsq/total_number_of_events);
                    Cnmk_os_err = Cnmk_os_err/num_pair_os2;
                }
                output_os2 << scientific << setw(18) << setprecision(8)
                           << Cnmk_os_avg << "  " << Cnmk_os_err << "  ";
            }
            output_ss1 << endl;
            output_os1 << endl;
            output_ss2 << endl;
            output_os2 << endl;
        }
        output_ss1.close();
        output_os1.close();
        output_ss2.close();
        output_os2.close();
    }
}

//! This function outputs the event averaged four-particle correlation
void singleParticleSpectra::output_four_particle_SC_correlation() {
    ostringstream filename;
    if (rap_type == 0) {
        filename << path << "/particle_" << particle_monval << "_SCmn"
                 << "_eta_" << rap_min << "_" << rap_max << ".dat";
    } else {
        filename << path << "/particle_" << particle_monval << "_SCmn"
                 << "_y_" << rap_min << "_" << rap_max << ".dat";
    }
    ofstream output(filename.str().c_str());
    output << "# n  SC_mn  SC_mn_err" << endl;
    double num_pair = SC_mn[0]/total_number_of_events;
    double num_pair_stdsq = (
            SC_mn_err[0]/total_number_of_events - num_pair*num_pair);
    double num_pair_err = 0.0;
    if (num_pair_stdsq > 0) {
        num_pair_err = sqrt(num_pair_stdsq/total_number_of_events);
    }
    output << scientific << setw(18) << setprecision(8)
           << 0 << "  " << num_pair << "  " << num_pair_err << endl;
    for (int i = 1; i < SC_num_corr; i++) {
        double SC_mn_avg = SC_mn[i]/total_number_of_events;
        double SC_mn_stdsq = (
                SC_mn_err[i]/total_number_of_events - SC_mn_avg*SC_mn_avg);
        SC_mn_avg = SC_mn_avg/num_pair;
        double SC_mn_avg_err = 0.0;
        if (SC_mn_stdsq > 0) {
            SC_mn_avg_err = sqrt(SC_mn_stdsq/total_number_of_events);
            SC_mn_avg_err = SC_mn_avg_err/num_pair;
        }
        output << scientific << setw(18) << setprecision(8)
               << i << "  " << SC_mn_avg << "  " << SC_mn_avg_err << endl;
    }
}

//! this function computes the pT-integrated Qn vector as a function of
//! rapidity in one event
void singleParticleSpectra::calculate_rapidity_distribution(int event_id,
        double **event_Qn_real, double **event_Qn_real_err,
        double **event_Qn_imag, double **event_Qn_imag_err, int flag) {
    // first clean the results arrays
    for (int i = 0; i < N_rap; i++) {
        for (int j = 0; j < order_max; j++) {
            event_Qn_real[i][j] = 0.0;
            event_Qn_real_err[i][j] = 0.0;
            event_Qn_imag[i][j] = 0.0;
            event_Qn_imag_err[i][j] = 0.0;
        }
    }

    int number_of_particles = 0;
    if (flag == 0) {
        number_of_particles = particle_list->get_number_of_particles(event_id);
    } else if (flag == 1) {
        number_of_particles = (
                particle_list->get_number_of_positive_particles(event_id));
    } else if (flag == 2) {
        number_of_particles = (
                particle_list->get_number_of_negative_particles(event_id));
    } else {
        cout << "[Error]:calculate_rapidity_distribution: "
             << "unrecogonized flag = " << flag << endl;
        exit(1);
    }
    for (int i = 0; i < number_of_particles; i++) {
        double pz_local = 0.0;
        double E_local = 0.0;
        if (flag == 0) {
            pz_local = particle_list->get_particle(event_id, i).pz;
            E_local = particle_list->get_particle(event_id, i).E;
        } else if (flag == 1) {
            pz_local = particle_list->get_positive_particle(event_id, i).pz;
            E_local = particle_list->get_positive_particle(event_id, i).E;
        } else if (flag == 2) {
            pz_local = particle_list->get_negative_particle(event_id, i).pz;
            E_local = particle_list->get_negative_particle(event_id, i).E;
        }

        double rap_local = 0.0;
        if (rap_type == 0) {
            double mass = 0.0;
            if (flag == 0) {
                mass = particle_list->get_particle(event_id, i).mass;
            } else if (flag == 1) {
                mass = particle_list->get_positive_particle(event_id, i).mass;
            } else if (flag == 2) {
                mass = particle_list->get_negative_particle(event_id, i).mass;
            }
            double pmag = sqrt(E_local*E_local - mass*mass);
            rap_local = 0.5*log((pmag + pz_local)/(pmag - pz_local));
        } else {
            rap_local = 0.5*log((E_local + pz_local)/(E_local - pz_local));
        }
        
        if (rap_local > rapidity_dis_min && rap_local < rapidity_dis_max) {
            int rap_idx = static_cast<int>(
                    (rap_local - rapidity_dis_min)/drap);
            if (flag == 0) {
                rapidity_array[rap_idx] += rap_local;
                dNdy_array[rap_idx]++;
            }
            // calcualte vn
            double px_local = 0.0;
            double py_local = 0.0;
            if (flag == 0) {
                px_local = particle_list->get_particle(event_id, i).px;
                py_local = particle_list->get_particle(event_id, i).py;
            } else if (flag == 1) {
                px_local = particle_list->get_positive_particle(event_id, i).px;
                py_local = particle_list->get_positive_particle(event_id, i).py;
            } else if (flag == 2) {
                px_local = particle_list->get_negative_particle(event_id, i).px;
                py_local = particle_list->get_negative_particle(event_id, i).py;
            }
            double p_perp = sqrt(px_local*px_local + py_local*py_local);
            if (p_perp > vn_rapidity_dis_pT_min 
                && p_perp < vn_rapidity_dis_pT_max) {
                double p_phi = atan2(py_local, px_local);
                for (int iorder = 0; iorder < order_max; iorder++) {
                    double cos_nphi = cos(iorder*p_phi);
                    double sin_nphi = sin(iorder*p_phi);
                    event_Qn_real[rap_idx][iorder] += cos_nphi;
                    event_Qn_imag[rap_idx][iorder] += sin_nphi;
                    event_Qn_real_err[rap_idx][iorder] += cos_nphi*cos_nphi;
                    event_Qn_imag_err[rap_idx][iorder] += sin_nphi*sin_nphi;
                }
            }
        }
    }
}

void singleParticleSpectra::output_rapidity_distribution() {
    ostringstream filename;
    if (rap_type == 0) {
        filename << path << "/particle_" << particle_monval << "_dNdeta"
                 << "_pT_" << vn_rapidity_dis_pT_min << "_"
                 << vn_rapidity_dis_pT_max << ".dat";
    } else {
        filename << path << "/particle_" << particle_monval << "_dNdy"
                 << "_pT_" << vn_rapidity_dis_pT_min << "_"
                 << vn_rapidity_dis_pT_max << ".dat";
    }
    ofstream output(filename.str().c_str());
    output << "#y  dN/dy  dN/dy_err  vn_real  vn_real_err  " 
           << "vn_imag  vn_imag_err  vn_rms   vn_rms_err"
           << endl;
    for (int i = 0; i < N_rap; i++) {
        rapidity_array[i] = rapidity_array[i]/(dNdy_array[i] + 1.);
        double total_Nrap = dNdy_array[i];
        dNdy_array[i] = dNdy_array[i]/total_number_of_events;
        double dNdy_err = sqrt(dNdy_array[i]/total_number_of_events);

        if (particle_monval == 333) {
            total_Nrap = total_Nrap/reconst_branching_ratio;
            dNdy_array[i] = dNdy_array[i]/reconst_branching_ratio;
            dNdy_err = dNdy_err/reconst_branching_ratio;
        }

        output << scientific << setw(18) << setprecision(8)
               << rapidity_array[i] << "   " << dNdy_array[i]/drap << "   " 
               << dNdy_err/drap << "   ";

        for (int iorder = 1; iorder < order_max; iorder++) {
            double vn_evavg_real = (vn_real_rapidity_dis_array[i][iorder]
                              /(vn_real_rapidity_dis_array[i][0] + 1e-30));
            double vn_evavg_imag = (vn_imag_rapidity_dis_array[i][iorder]
                              /(vn_real_rapidity_dis_array[i][0] + 1e-30));
            double vn_real_err = (
                sqrt(vn_real_rapidity_dis_array_err[i][iorder]
                     /(vn_real_rapidity_dis_array[i][0] + 1e-30)
                     - vn_evavg_real*vn_evavg_real)
                /sqrt(vn_real_rapidity_dis_array[i][0] + 1e-30));
            double vn_imag_err = (
                sqrt(vn_imag_rapidity_dis_array_err[i][iorder]
                     /(vn_real_rapidity_dis_array[i][0] + 1e-30)
                     - vn_evavg_imag*vn_evavg_imag)
                /sqrt(vn_real_rapidity_dis_array[i][0] + 1e-30));
            if (isnan(vn_real_err)) {
                vn_real_err = 0.0;
            }
            if (isnan(vn_imag_err)) {
                vn_imag_err = 0.0;
            }

            double vn_rms = sqrt(vn_evavg_real*vn_evavg_real 
                                 + vn_evavg_imag*vn_evavg_imag);
            double vn_rms_err = (sqrt(pow(vn_evavg_real*vn_real_err, 2)
                                      + pow(vn_evavg_imag*vn_imag_err, 2))
                                 /(vn_rms + 1e-30));
         
            output << scientific << setw(18) << setprecision(8) 
                   << vn_evavg_real << "   " << vn_real_err << "   " 
                   << vn_evavg_imag << "   " << vn_imag_err << "   "
                   << vn_rms << "   " << vn_rms_err << "   ";
        }
        // output total number of particles in the last column
        // this quantities is useful when one wants to reconst the Qn vectors
        output << scientific << setw(18) << setprecision(8) << total_Nrap;
        output << endl;
    }
    output.close();
}

void singleParticleSpectra::check_dNdSV(int event_id) {
    int number_of_particles = particle_list->get_number_of_particles(event_id);
    for (int i = 0; i < number_of_particles; i++) {
        double pz_local = particle_list->get_particle(event_id, i).pz;
        double E_local = particle_list->get_particle(event_id, i).E;

        double rap_local;
        if (rap_type == 0) {
            double mass = particle_list->get_particle(event_id, i).mass;
            double pmag = sqrt(E_local*E_local - mass*mass);
            rap_local = 0.5*log((pmag + pz_local)/(pmag - pz_local));
        } else {
            rap_local = 0.5*log((E_local + pz_local)/(E_local - pz_local));
        }

        if (rap_local > rap_min && rap_local < rap_max) {
            // first dN/dtau
            double t_local = particle_list->get_particle(event_id, i).t;
            double z_local = particle_list->get_particle(event_id, i).z;
            double tau_local = sqrt(t_local*t_local - z_local*z_local);
            if (tau_local > tau_min && tau_local < tau_max) {
                double random = drand48()*intrinsic_dtau;
                int idx = (int)((tau_local + random - tau_min)/dtau);
                tau_array[idx] += tau_local;
                dNdtau_array[idx]++;
            }

            // second dN/dx
            double x_local = particle_list->get_particle(event_id, i).x;
            double y_local = particle_list->get_particle(event_id, i).y;
            if (fabs(y_local) < 0.5) {
                if (x_local > spatial_x_min && x_local < spatial_x_max) {
                    double random = drand48()*intrinsic_dx;
                    int idx = (int)((x_local + random - spatial_x_min)
                                    /dspatial_x);
                    xpt_array[idx] += x_local;
                    dNdx1_array[idx]++;
                }
            }
            if (fabs(x_local) < 0.5) {
                if (y_local > spatial_x_min && y_local < spatial_x_max) {
                    double random = drand48()*intrinsic_dx;
                    int idx = (int)((y_local + random - spatial_x_min)
                                    /dspatial_x);
                    ypt_array[idx] += y_local;
                    dNdx2_array[idx]++;
                }
            }

            // dN/(dtau dx)
            if (fabs(y_local) < 0.5) {
                if (tau_local > tau_min && tau_local < tau_max) {
                    double random = drand48()*intrinsic_dtau;
                    int idx_tau = (int)((tau_local + random - tau_min)/dtau);
                    if (x_local > spatial_x_min && x_local < spatial_x_max) {
                        double random = drand48()*intrinsic_dx;
                        int idx_x = (int)((x_local + random - spatial_x_min)
                                          /dspatial_x);
                        dNdtaudx1_array[idx_tau][idx_x]++;
                    }
                }
            }
            if (fabs(x_local) < 0.5) {
                if (tau_local > tau_min && tau_local < tau_max) {
                    double random = drand48()*intrinsic_dtau;
                    int idx_tau = (int)((tau_local + random - tau_min)/dtau);
                    if (y_local > spatial_x_min && y_local < spatial_x_max) {
                        double random = drand48()*intrinsic_dx;
                        int idx_x = (int)((y_local + random - spatial_x_min)
                                          /dspatial_x);
                        dNdtaudx2_array[idx_tau][idx_x]++;
                    }
                }
            }

            // third dN/deta_s
            double etas_local = (0.5*log((t_local + z_local)
                                         /(t_local - z_local)));
            double y_minus_etas = rap_local - etas_local;
            if (y_minus_etas > eta_s_min && y_minus_etas < eta_s_max) {
                double random = drand48()*intrinsic_detas;
                int idx = (int)((y_minus_etas + random - eta_s_min)/deta_s);
                eta_s_array[idx] += y_minus_etas;
                dNdetas_array[idx]++;
            }
        }
    }
}

void singleParticleSpectra::output_dNdSV() {
    // first dN/dtau
    ostringstream filename;
    filename << path << "/check_" << particle_monval << "_dNdtau.dat";
    ofstream output(filename.str().c_str());
    for (int i = 0; i < N_tau; i++) {
        tau_array[i] = tau_array[i]/(dNdtau_array[i] + 1.);
        dNdtau_array[i] = dNdtau_array[i]/total_number_of_events;
        double dNdtau_err = sqrt(dNdtau_array[i]/total_number_of_events);
        if (particle_monval == 333) {
            dNdtau_array[i] = dNdtau_array[i]/reconst_branching_ratio;
            dNdtau_err =  dNdtau_err/reconst_branching_ratio;
        }
        output << scientific << setw(18) << setprecision(8)
               << tau_array[i] << "   " << dNdtau_array[i]/dtau << "   " 
               << dNdtau_err/dtau << endl;
    }
    output.close();
    
    // second dN/dx
    ostringstream filename2;
    filename2 << path << "/check_" << particle_monval << "_dNdx.dat";
    ofstream output2(filename2.str().c_str());
    for (int i = 0; i < N_xpt; i++) {
        xpt_array[i] = xpt_array[i]/(dNdx1_array[i] + 1.);
        ypt_array[i] = ypt_array[i]/(dNdx2_array[i] + 1.);
        dNdx1_array[i] = dNdx1_array[i]/total_number_of_events;
        dNdx2_array[i] = dNdx2_array[i]/total_number_of_events;
        double dNdx1_err = sqrt(dNdx1_array[i]/total_number_of_events);
        double dNdx2_err = sqrt(dNdx2_array[i]/total_number_of_events);
        if (particle_monval == 333) {
            dNdx1_array[i] = dNdx1_array[i]/reconst_branching_ratio;
            dNdx2_array[i] = dNdx2_array[i]/reconst_branching_ratio;
            dNdx1_err = dNdx1_err/reconst_branching_ratio;
            dNdx2_err = dNdx2_err/reconst_branching_ratio;
        }
        output2 << scientific << setw(18) << setprecision(8)
                << xpt_array[i] << "   " << dNdx1_array[i]/dspatial_x << "   " 
                << dNdx1_err/dspatial_x << "   " 
                << ypt_array[i] << "   " << dNdx2_array[i]/dspatial_x << "   " 
                << dNdx2_err/dspatial_x << endl;
    }
    output2.close();
    
    // dN/(dtau dx)
    ostringstream filename2_1, filename2_2;
    filename2_1 << path << "/check_" << particle_monval << "_dNdtaudx1.dat";
    filename2_2 << path << "/check_" << particle_monval << "_dNdtaudx2.dat";
    ofstream output2_1(filename2_1.str().c_str());
    ofstream output2_2(filename2_2.str().c_str());
    for (int i = 0; i < N_tau; i++) {
        for (int j = 0; j < N_xpt; j++) {
            dNdtaudx1_array[i][j] = (
                            dNdtaudx1_array[i][j]/total_number_of_events);
            dNdtaudx2_array[i][j] = (
                            dNdtaudx2_array[i][j]/total_number_of_events);
            if (particle_monval == 333) {
                dNdtaudx1_array[i][j] = (
                            dNdtaudx1_array[i][j]/reconst_branching_ratio);
                dNdtaudx2_array[i][j] = (
                            dNdtaudx2_array[i][j]/reconst_branching_ratio);
            }
            output2_1 << scientific << setw(18) << setprecision(8)
                      << dNdtaudx1_array[i][j]/dspatial_x/dtau << "   ";
            output2_2 << scientific << setw(18) << setprecision(8)
                      << dNdtaudx2_array[i][j]/dspatial_x/dtau << "   ";
        }
        output2_1 << endl;
        output2_2 << endl;
    }
    output2_1.close();
    output2_2.close();

    // third dN/detas
    ostringstream filename3;
    filename3 << path << "/check_" << particle_monval << "_dNdetas.dat";
    ofstream output3(filename3.str().c_str());
    for (int i = 0; i < N_eta_s; i++) {
        eta_s_array[i] = eta_s_array[i]/(dNdetas_array[i] + 1.);
        dNdetas_array[i] = dNdetas_array[i]/total_number_of_events;
        double dNdetas_err = sqrt(dNdetas_array[i]/total_number_of_events);
        if (particle_monval == 333) {
            dNdetas_array[i] = dNdetas_array[i]/reconst_branching_ratio;
            dNdetas_err = dNdetas_err/reconst_branching_ratio;
        }
        output3 << scientific << setw(18) << setprecision(8)
                << eta_s_array[i] << "   " << dNdetas_array[i]/deta_s << "   " 
                << dNdetas_err/deta_s << endl;
    }
    output3.close();
}

