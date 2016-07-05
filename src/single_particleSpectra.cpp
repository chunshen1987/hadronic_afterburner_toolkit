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
    Qn_vector_real = new double [order_max];
    Qn_vector_imag = new double [order_max];
    Qn_vector_real_err = new double [order_max];
    Qn_vector_imag_err = new double [order_max];
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
        Qn_diff_vector_real[i] = new double [npT];
        Qn_diff_vector_imag[i] = new double [npT];
        Qn_diff_vector_real_err[i] = new double [npT];
        Qn_diff_vector_imag_err[i] = new double [npT];
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
}

void singleParticleSpectra::calculate_Qn_vector_shell() {
    int event_id = 0;
    int buffer_size = particle_list->get_event_buffer_size();
    while (!particle_list->end_of_file()) {
        cout << "Reading event: " << event_id+1 
             << "-" << event_id + buffer_size << " ... " << flush;
        particle_list->read_in_particle_samples();
        cout << " processing ..." << flush;
        int nev = particle_list->get_number_of_events();
        for (int iev = 0; iev < nev; iev++) {
            event_id++;
            calculate_Qn_vector(iev);

            if (rapidity_distribution_flag == 1)
                calculate_rapidity_distribution(iev);

            if (check_spatial_flag == 1)
                check_dNdSV(iev);
        }
        cout << " done!" << endl;
    }
    total_number_of_events = event_id;
    output_Qn_vectors();
    if (rapidity_distribution_flag == 1)
        output_rapidity_distribution();
    if (check_spatial_flag == 1)
        output_dNdSV();
}

void singleParticleSpectra::calculate_Qn_vector(int event_id) {
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
            if (p_perp > pT_min && p_perp < pT_max) {
                double p_phi = atan2(py_local, px_local);
                int p_idx = (int)((p_perp - pT_min)/dpT);
                pT_mean_array[p_idx] += p_perp;
                pT_mean_array_err[p_idx] += p_perp*p_perp;
                for (int iorder = 0; iorder < order_max; iorder++) {
                    double cos_nphi = cos(iorder*p_phi);
                    double sin_nphi = sin(iorder*p_phi);
                    Qn_vector_real[iorder] += cos_nphi;
                    Qn_vector_imag[iorder] += sin_nphi;
                    Qn_vector_real_err[iorder] += cos_nphi*cos_nphi;
                    Qn_vector_imag_err[iorder] += sin_nphi*sin_nphi;
                    Qn_diff_vector_real[iorder][p_idx] += cos_nphi;
                    Qn_diff_vector_imag[iorder][p_idx] += sin_nphi;
                    Qn_diff_vector_real_err[iorder][p_idx] += cos_nphi*cos_nphi;
                    Qn_diff_vector_imag_err[iorder][p_idx] += sin_nphi*sin_nphi;
                }
            }
        }
    }
}

void singleParticleSpectra::output_Qn_vectors() {
    // pT-integrated flow
    ostringstream filename;
    filename << path << "/particle_" << particle_monval << "_vndata.dat";
    ofstream output(filename.str().c_str());

    double dN_ev_avg = Qn_vector_real[0]/total_number_of_events;
    double dN_ev_avg_err = sqrt(dN_ev_avg/total_number_of_events);
    if (particle_monval == 333) {
        // for phi(1020) need to rescale the yield by
        // reconstruction branching ratio
        dN_ev_avg = dN_ev_avg/reconst_branching_ratio;
        dN_ev_avg_err = dN_ev_avg_err/reconst_branching_ratio;
    }
    output << scientific << setw(18) << setprecision(8) 
           << 0 << "   " << dN_ev_avg << "   " << dN_ev_avg_err << "   " 
           << 0.0 << "   " << 0.0 << endl;
    for (int iorder = 1; iorder < order_max; iorder++) {
        double vn_evavg_real = Qn_vector_real[iorder]/Qn_vector_real[0];
        double vn_evavg_imag = Qn_vector_imag[iorder]/Qn_vector_real[0];
        double vn_real_err = sqrt(Qn_vector_real_err[iorder]/Qn_vector_real[0] 
                                  - vn_evavg_real*vn_evavg_real)
                             /sqrt(Qn_vector_real[0]);
        double vn_imag_err = sqrt(Qn_vector_imag_err[iorder]/Qn_vector_real[0] 
                                  - vn_evavg_imag*vn_evavg_imag)
                             /sqrt(Qn_vector_real[0]);

        output << scientific << setw(18) << setprecision(8) << iorder << "   " 
               << vn_evavg_real << "   " << vn_real_err << "   " 
               << vn_evavg_imag << "   " << vn_imag_err << endl;
    }
    output.close();
    
    // pT-differential flow
    ostringstream filename_diff;
    filename_diff << path 
                  << "/particle_" << particle_monval << "_vndata_diff.dat";
    ofstream output_diff(filename_diff.str().c_str());

    for (int ipT = 0; ipT < npT - 1; ipT++) {
        double dNpT_ev_avg = Qn_diff_vector_real[0][ipT]/total_number_of_events;
        double dNpT_ev_avg_err = sqrt(dNpT_ev_avg/total_number_of_events);
        if (particle_monval == 333) {
            // for phi(1020) need to rescale the yield by
            // reconstruction branching ratio
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
                    << dNpT_ev_avg/mean_pT/dpT/(2*M_PI) << "   " 
                    << dNpT_ev_avg_err/mean_pT/dpT/(2*M_PI);
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
                output_diff << scientific << setw(18) << setprecision(8) 
                            << vn_evavg_real << "   " 
                            << vn_evavg_real_err << "   " 
                            << vn_evavg_imag << "   " << vn_evavg_imag_err;
            } else {
                output_diff << scientific << setw(18) << setprecision(8) 
                            << 0.0e0 << "   " << 0.0e0 << "   " 
                            << 0.0e0 << "   " << 0.0e0;
            }
        }
        output_diff << endl;
    }
    output_diff.close();
}

void singleParticleSpectra::calculate_rapidity_distribution(int event_id) {
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
        
        if (rap_local > rapidity_dis_min && rap_local < rapidity_dis_max) {
            int rap_idx = (int)((rap_local - rapidity_dis_min)/drap);
            rapidity_array[rap_idx] += rap_local;
            dNdy_array[rap_idx]++;
            
            // calcualte vn
            double px_local = particle_list->get_particle(event_id, i).px;
            double py_local = particle_list->get_particle(event_id, i).py;
            double p_perp = sqrt(px_local*px_local + py_local*py_local);
            if(p_perp > vn_rapidity_dis_pT_min 
               && p_perp < vn_rapidity_dis_pT_max) {
                double p_phi = atan2(py_local, px_local);
                for (int iorder = 0; iorder < order_max; iorder++) {
                    double cos_nphi = cos(iorder*p_phi);
                    double sin_nphi = sin(iorder*p_phi);
                    vn_real_rapidity_dis_array[rap_idx][iorder] += cos_nphi;
                    vn_imag_rapidity_dis_array[rap_idx][iorder] += sin_nphi;
                    vn_real_rapidity_dis_array_err[rap_idx][iorder] += (
                                                            cos_nphi*cos_nphi);
                    vn_imag_rapidity_dis_array_err[rap_idx][iorder] += (
                                                            sin_nphi*sin_nphi);
                }
            }
        }
    }
}

void singleParticleSpectra::output_rapidity_distribution() {
    ostringstream filename;
    if (rap_type == 0)
        filename << path << "/particle_" << particle_monval << "_dNdeta.dat";
    else
        filename << path << "/particle_" << particle_monval << "_dNdy.dat";
    ofstream output(filename.str().c_str());
    output << "#y  dN/dy  dN/dy_err  vn_real  vn_real_err  " 
           << "vn_imag  vn_imag_err  vn_rms   vn_rms_err"
           << endl;
    for (int i = 0; i < N_rap; i++) {
        rapidity_array[i] = rapidity_array[i]/(dNdy_array[i] + 1.);
        dNdy_array[i] = dNdy_array[i]/total_number_of_events;
        double dNdy_err = sqrt(dNdy_array[i]/total_number_of_events);

        if (particle_monval == 333) {
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

