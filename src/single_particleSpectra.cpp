#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <iomanip>
#include "parameters.h"
#include "single_particleSpectra.h"

using namespace std;

singleParticleSpectra::singleParticleSpectra(ParameterReader *paraRdr_in, string path_in, particleSamples *particle_list_in)
{
    paraRdr = paraRdr_in;
    path = path_in;
    particle_list = particle_list_in;

    particle_monval = paraRdr->getVal("particle_monval");

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
    for(int i = 0; i < npT; i++)
    {
        pT_array[i] = pT_min + dpT*i;
        pT_mean_array[i] = 0.0;
        pT_mean_array_err[i] = 0.0;
    }
    for(int i = 0; i < order_max; i++)
    {
        Qn_vector_real[i] = 0.0;
        Qn_vector_imag[i] = 0.0;
        Qn_vector_real_err[i] = 0.0;
        Qn_vector_imag_err[i] = 0.0;
        Qn_diff_vector_real[i] = new double [npT];
        Qn_diff_vector_imag[i] = new double [npT];
        Qn_diff_vector_real_err[i] = new double [npT];
        Qn_diff_vector_imag_err[i] = new double [npT];
        for(int j = 0; j < npT; j++)
        {
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

    // check dN/dtau distribution
    check_spatial_flag = paraRdr->getVal("check_spatial_dis");
    if(check_spatial_flag == 1)
    {
        // dN/dtau
        intrinsic_dtau = paraRdr->getVal("intrinsic_dtau");
        N_tau = 50;
        tau_min = 0.6;
        tau_max = 15.0;
        dtau = (tau_max - tau_min)/(N_tau - 1);
        tau_array = new double [N_tau];
        dNdtau_array = new double [N_tau];
        for(int i = 0; i < N_tau; i++)
        {
            tau_array[i] = 0.0;
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
        for(int i = 0; i < N_xpt; i++)
        {
            xpt_array[i] = 0.0;
            ypt_array[i] = 0.0;
            dNdx1_array[i] = 0.0;
            dNdx2_array[i] = 0.0;
        }

        // dN/deta_s
        intrinsic_detas = paraRdr->getVal("intrinsic_detas");
        N_eta_s = 40;
        eta_s_min = - 3.0;
        eta_s_max = 3.0;
        deta_s = (eta_s_max - eta_s_min)/(N_eta_s - 1);
        eta_s_array = new double [N_eta_s];
        dNdetas_array = new double [N_eta_s];
        for(int i = 0; i < N_eta_s; i++)
        {
            eta_s_array[i] = 0.0;
            dNdetas_array[i] = 0.0;
        }

    }
}

singleParticleSpectra::~singleParticleSpectra()
{
    delete [] pT_array;
    delete [] pT_mean_array;
    delete [] pT_mean_array_err;
    delete [] Qn_vector_real;
    delete [] Qn_vector_imag;
    delete [] Qn_vector_real_err;
    delete [] Qn_vector_imag_err;
    for(int i = 0; i < order_max; i++)
    {
        delete [] Qn_diff_vector_real[i];
        delete [] Qn_diff_vector_imag[i];
        delete [] Qn_diff_vector_real_err[i];
        delete [] Qn_diff_vector_imag_err[i];
    }
    delete [] Qn_diff_vector_real;
    delete [] Qn_diff_vector_imag;
    delete [] Qn_diff_vector_real_err;
    delete [] Qn_diff_vector_imag_err;

    if(check_spatial_flag == 1)
    {
        delete [] tau_array;
        delete [] dNdtau_array;
        delete [] xpt_array;
        delete [] ypt_array;
        delete [] dNdx1_array;
        delete [] dNdx2_array;
        delete [] eta_s_array;
        delete [] dNdetas_array;
    }
}

void singleParticleSpectra::calculate_Qn_vector_shell()
{
    int event_id = 0;
    int buffer_size = particle_list->get_event_buffer_size();
    while(!particle_list->end_of_file())
    {
        cout << "Reading event: " << event_id+1 << "-" << event_id + buffer_size << " ... " << flush;
        particle_list->read_in_particle_samples();
        cout << " processing ..." << flush;
        int nev = particle_list->get_number_of_events();
        for(int iev = 0; iev < nev; iev++)
        {
            event_id++;
            int number_of_particles = particle_list->get_number_of_particles(iev);
            calculate_Qn_vector(iev);

            if(check_spatial_flag == 1)
                check_dNdSV(iev);
        }
        cout << " done!" << endl;
    }
    total_number_of_events = event_id;
    output_Qn_vectors();
    if(check_spatial_flag == 1)
        output_dNdSV();
}

void singleParticleSpectra::calculate_Qn_vector(int event_id)
{
    int number_of_particles = particle_list->get_number_of_particles(event_id);

    for(int i = 0; i < number_of_particles; i++)
    {
        double pz_local = particle_list->get_particle(event_id, i).pz;
        double E_local = particle_list->get_particle(event_id, i).E;

        double rap_local;
        if(rap_type == 0)
        {
            double mass = particle_list->get_particle(event_id, i).mass;
            double pmag = sqrt(E_local*E_local - mass*mass);
            rap_local = 0.5*log((pmag + pz_local)/(pmag - pz_local));
        }
        else
            rap_local = 0.5*log((E_local + pz_local)/(E_local - pz_local));

        if(rap_local > rap_min && rap_local < rap_max)
        {
            double px_local = particle_list->get_particle(event_id, i).px;
            double py_local = particle_list->get_particle(event_id, i).py;
            double p_perp = sqrt(px_local*px_local + py_local*py_local);
            if(p_perp > pT_min && p_perp < pT_max)
            {
                double p_phi = atan2(py_local, px_local);
                int p_idx = (int)((p_perp - pT_min)/dpT);
                pT_mean_array[p_idx] += p_perp;
                pT_mean_array_err[p_idx] += p_perp*p_perp;
                for(int iorder = 0; iorder < order_max; iorder++)
                {
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

void singleParticleSpectra::output_Qn_vectors()
{
    // pT-integrated flow
    ostringstream filename;
    filename << path << "/particle_" << particle_monval << "_vndata.dat";
    ofstream output(filename.str().c_str());

    double dN_ev_avg = Qn_vector_real[0]/total_number_of_events;
    double dN_ev_avg_err = sqrt(dN_ev_avg/total_number_of_events);
    output << scientific << setw(18) << setprecision(8) 
           << 0 << "   " << dN_ev_avg << "   " << dN_ev_avg_err << "   " 
           << 0.0 << "   " << 0.0 << endl;
    for(int iorder = 1; iorder < order_max; iorder++)
    {
        double vn_evavg_real = Qn_vector_real[iorder]/Qn_vector_real[0];
        double vn_evavg_imag = Qn_vector_imag[iorder]/Qn_vector_real[0];
        double vn_real_err = sqrt(Qn_vector_real_err[iorder]/Qn_vector_real[0] - vn_evavg_real*vn_evavg_real)/sqrt(Qn_vector_real[0]);
        double vn_imag_err = sqrt(Qn_vector_imag_err[iorder]/Qn_vector_real[0] - vn_evavg_imag*vn_evavg_imag)/sqrt(Qn_vector_real[0]);

        output << scientific << setw(18) << setprecision(8) << iorder << "   " 
               << vn_evavg_real << "   " << vn_real_err << "   " 
               << vn_evavg_imag << "   " << vn_imag_err << endl;
    }
    output.close();
    
    // pT-differential flow
    ostringstream filename_diff;
    filename_diff << path << "/particle_" << particle_monval << "_vndata_diff.dat";
    ofstream output_diff(filename_diff.str().c_str());

    for(int ipT = 0; ipT < npT - 1; ipT++)
    {
        double dNpT_ev_avg = Qn_diff_vector_real[0][ipT]/total_number_of_events;
        double dNpT_ev_avg_err = sqrt(dNpT_ev_avg/total_number_of_events);
        double mean_pT = pT_mean_array[ipT]/Qn_diff_vector_real[0][ipT];
        double mean_pT_err = (pT_mean_array_err[ipT]/Qn_diff_vector_real[0][ipT] - mean_pT*mean_pT)/sqrt(Qn_diff_vector_real[0][ipT]);
        output_diff << scientific << setw(18) << setprecision(8) 
                    << mean_pT << "   " << mean_pT_err << "   " 
                    << dNpT_ev_avg/mean_pT/dpT/(2*M_PI) << "   " 
                    << dNpT_ev_avg_err/mean_pT/dpT/(2*M_PI);
        for(int iorder = 1; iorder < order_max; iorder++)
        {
            double vn_evavg_real = Qn_diff_vector_real[iorder][ipT]/Qn_diff_vector_real[0][ipT];
            double vn_evavg_imag = Qn_diff_vector_imag[iorder][ipT]/Qn_diff_vector_real[0][ipT];
            double vn_evavg_real_err = sqrt(Qn_diff_vector_real_err[iorder][ipT]/Qn_diff_vector_real[0][ipT] - vn_evavg_real*vn_evavg_real)/sqrt(Qn_diff_vector_real[0][ipT]);
            double vn_evavg_imag_err = sqrt(Qn_diff_vector_imag_err[iorder][ipT]/Qn_diff_vector_real[0][ipT] - vn_evavg_imag*vn_evavg_imag)/sqrt(Qn_diff_vector_real[0][ipT]);
            output_diff << scientific << setw(18) << setprecision(8) 
                        << vn_evavg_real << "   " << vn_evavg_real_err << "   " 
                        << vn_evavg_imag << "   " << vn_evavg_imag_err;
        }
        output_diff << endl;
    }
    output_diff.close();
}

void singleParticleSpectra::check_dNdSV(int event_id)
{
    int number_of_particles = particle_list->get_number_of_particles(event_id);
    for(int i = 0; i < number_of_particles; i++)
    {
        double pz_local = particle_list->get_particle(event_id, i).pz;
        double E_local = particle_list->get_particle(event_id, i).E;

        double rap_local;
        if(rap_type == 0)
        {
            double mass = particle_list->get_particle(event_id, i).mass;
            double pmag = sqrt(E_local*E_local - mass*mass);
            rap_local = 0.5*log((pmag + pz_local)/(pmag - pz_local));
        }
        else
            rap_local = 0.5*log((E_local + pz_local)/(E_local - pz_local));

        if(rap_local > rap_min && rap_local < rap_max)
        {
            // first dN/dtau
            double t_local = particle_list->get_particle(event_id, i).t;
            double z_local = particle_list->get_particle(event_id, i).z;
            double tau_local = sqrt(t_local*t_local - z_local*z_local);
            if(tau_local > tau_min && tau_local < tau_max)
            {
                double random = drand48()*intrinsic_dtau;
                int idx = (int)((tau_local + random - tau_min)/dtau);
                tau_array[idx] += tau_local;
                dNdtau_array[idx]++;
            }
            // second dN/dx
            double x_local = particle_list->get_particle(event_id, i).x;
            double y_local = particle_list->get_particle(event_id, i).y;
            if(fabs(y_local) < 0.5)
            {
                if(x_local > spatial_x_min && x_local < spatial_x_max)
                {
                    double random = drand48()*intrinsic_dx;
                    int idx = (int)((x_local + random - spatial_x_min)/dspatial_x);
                    xpt_array[idx] += x_local;
                    dNdx1_array[idx]++;
                }
            }
            if(fabs(x_local) < 0.5)
            {
                if(y_local > spatial_x_min && y_local < spatial_x_max)
                {
                    double random = drand48()*intrinsic_dx;
                    int idx = (int)((y_local + random - spatial_x_min)/dspatial_x);
                    ypt_array[idx] += y_local;
                    dNdx2_array[idx]++;
                }
            }
            // third dN/deta_s
            double etas_local = 0.5*log((t_local + z_local)/(t_local - z_local));
            double y_minus_etas = rap_local - etas_local;
            if(y_minus_etas > eta_s_min && y_minus_etas < eta_s_max)
            {
                double random = drand48()*intrinsic_detas;
                int idx = (int)((y_minus_etas + random - eta_s_min)/deta_s);
                eta_s_array[idx] += y_minus_etas;
                dNdetas_array[idx]++;
            }
        }
    }
}

void singleParticleSpectra::output_dNdSV()
{
    // first dN/dtau
    ostringstream filename;
    filename << path << "/check_" << particle_monval << "_dNdtau.dat";
    ofstream output(filename.str().c_str());
    for(int i = 0; i < N_tau; i++)
    {
        tau_array[i] = tau_array[i]/(dNdtau_array[i] + 1e-15);
        dNdtau_array[i] = dNdtau_array[i]/total_number_of_events;
        double dNdtau_err = sqrt(dNdtau_array[i]/total_number_of_events);
        output << scientific << setw(18) << setprecision(8)
               << tau_array[i] << "   " << dNdtau_array[i]/dtau << "   " 
               << dNdtau_err/dtau << endl;
    }
    output.close();
    
    // second dN/dx
    ostringstream filename2;
    filename2 << path << "/check_" << particle_monval << "_dNdx.dat";
    ofstream output2(filename2.str().c_str());
    for(int i = 0; i < N_xpt; i++)
    {
        xpt_array[i] = xpt_array[i]/(dNdx1_array[i] + 1e-15);
        ypt_array[i] = ypt_array[i]/(dNdx2_array[i] + 1e-15);
        dNdx1_array[i] = dNdx1_array[i]/total_number_of_events;
        dNdx2_array[i] = dNdx2_array[i]/total_number_of_events;
        double dNdx1_err = sqrt(dNdx1_array[i]/total_number_of_events);
        double dNdx2_err = sqrt(dNdx2_array[i]/total_number_of_events);
        output2 << scientific << setw(18) << setprecision(8)
                << xpt_array[i] << "   " << dNdx1_array[i]/dspatial_x << "   " << dNdx1_err/dspatial_x << "   " 
                << ypt_array[i] << "   " << dNdx2_array[i]/dspatial_x << "   " << dNdx2_err/dspatial_x << endl;
    }
    output2.close();

    // third dN/detas
    ostringstream filename3;
    filename3 << path << "/check_" << particle_monval << "_dNdetas.dat";
    ofstream output3(filename3.str().c_str());
    for(int i = 0; i < N_eta_s; i++)
    {
        eta_s_array[i] = eta_s_array[i]/(dNdetas_array[i] + 1e-15);
        dNdetas_array[i] = dNdetas_array[i]/total_number_of_events;
        double dNdetas_err = sqrt(dNdetas_array[i]/total_number_of_events);
        output3 << scientific << setw(18) << setprecision(8)
                << eta_s_array[i] << "   " << dNdetas_array[i]/deta_s << "   " 
                << dNdetas_err/deta_s << endl;
    }
    output3.close();
}

