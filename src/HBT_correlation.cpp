#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <iomanip>
#include "parameters.h"
#include "HBT_correlation.h"
using namespace std;

HBT_correlation::HBT_correlation(ParameterReader* paraRdr_in, string path_in, particleSamples *particle_list_in)
{
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
    for(int i = 0; i < qnpts; i++)
    {
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

    dKT = (KT_max - KT_min)/(n_KT - 1);
    dKphi = 2*M_PI/(n_Kphi - 1);
    
    KT_array = new double [n_KT];
    number_of_pairs_numerator_KTdiff = new unsigned long long int [n_KT];
    number_of_pairs_denormenator_KTdiff = new unsigned long long int [n_KT];
    for(int i = 0; i < n_KT; i++)
    {
        KT_array[i] = KT_min + i*dKT;
        number_of_pairs_numerator_KTdiff[i] = 0;
        number_of_pairs_denormenator_KTdiff[i] = 0;
    }

    Kphi_array = new double [n_Kphi];
    for(int i = 0; i < n_Kphi; i++)
        Kphi_array[i] = 0.0 + i*dKphi;

    number_of_pairs_numerator_KTKphidiff = new unsigned long long int* [n_KT];
    number_of_pairs_denormenator_KTKphidiff = new unsigned long long int* [n_KT];
    for(int i = 0; i < n_KT; i++)
    {
        number_of_pairs_numerator_KTKphidiff[i] = new unsigned long long int [n_Kphi];
        number_of_pairs_denormenator_KTKphidiff[i] = new unsigned long long int [n_Kphi];
        for(int j = 0; j < n_Kphi; j++)
        {
            number_of_pairs_numerator_KTKphidiff[i] = 0;
            number_of_pairs_denormenator_KTKphidiff[i] = 0;
        }
    }
    

    number_of_mixed_events = paraRdr->getVal("number_of_mixed_events");
    number_of_oversample_events = paraRdr->getVal("number_of_oversample_events");
    number_pairs_num = 0;
    number_pairs_denorm = 0;
    if(azimuthal_flag == 0)
    {
        q_out_mean = new double ***[n_KT];
        q_side_mean = new double ***[n_KT];
        q_long_mean = new double ***[n_KT];
        correl_3d_num = new double *** [n_KT];
        correl_3d_num_count = new double *** [n_KT];
        correl_3d_denorm = new double *** [n_KT];
        for(int iK = 0; iK < n_KT; iK++)
        {
            q_out_mean[iK] = new double ** [qnpts];
            q_side_mean[iK] = new double ** [qnpts];
            q_long_mean[iK] = new double ** [qnpts];
            correl_3d_num[iK] = new double ** [qnpts];
            correl_3d_num_count[iK] = new double ** [qnpts];
            correl_3d_denorm[iK] = new double ** [qnpts];
            for(int i = 0; i < qnpts; i++)
            {
                q_out_mean[iK][i] = new double * [qnpts];
                q_side_mean[iK][i] = new double * [qnpts];
                q_long_mean[iK][i] = new double * [qnpts];
                correl_3d_num[iK][i] = new double * [qnpts];
                correl_3d_num_count[iK][i] = new double * [qnpts];
                correl_3d_denorm[iK][i] = new double * [qnpts];
                for(int j = 0; j < qnpts; j++)
                {
                    q_out_mean[iK][i][j] = new double [qnpts];
                    q_side_mean[iK][i][j] = new double [qnpts];
                    q_long_mean[iK][i][j] = new double [qnpts];
                    correl_3d_num[iK][i][j] = new double [qnpts];
                    correl_3d_num_count[iK][i][j] = new double [qnpts];
                    correl_3d_denorm[iK][i][j] = new double [qnpts];
                    for(int k = 0; k < qnpts; k++)
                    {
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
    }
    else
    {
        q_out_diff_mean = new double **** [n_KT];
        q_side_diff_mean = new double **** [n_KT];
        q_long_diff_mean = new double **** [n_KT];
        correl_3d_Kphi_diff_num = new double **** [n_KT];
        correl_3d_Kphi_diff_num_count = new double **** [n_KT];
        correl_3d_Kphi_diff_denorm = new double **** [n_KT];
        for(int iK = 0; iK < n_KT; iK++)
        {
            q_out_diff_mean[iK] = new double *** [n_Kphi];
            q_side_diff_mean[iK] = new double *** [n_Kphi];
            q_long_diff_mean[iK] = new double *** [n_Kphi];
            correl_3d_Kphi_diff_num[iK] = new double *** [n_Kphi];
            correl_3d_Kphi_diff_num_count[iK] = new double *** [n_Kphi];
            correl_3d_Kphi_diff_denorm[iK] = new double *** [n_Kphi];
            for(int iphi = 0; iphi < n_Kphi; iphi++)
            {
                q_out_diff_mean[iK][iphi] = new double ** [qnpts];
                q_side_diff_mean[iK][iphi] = new double ** [qnpts];
                q_long_diff_mean[iK][iphi] = new double ** [qnpts];
                correl_3d_Kphi_diff_num[iK][iphi] = new double ** [qnpts];
                correl_3d_Kphi_diff_num_count[iK][iphi] = new double ** [qnpts];
                correl_3d_Kphi_diff_denorm[iK][iphi] = new double ** [qnpts];
                for(int i = 0; i < qnpts; i++)
                {
                    q_out_diff_mean[iK][iphi][i] = new double * [qnpts];
                    q_side_diff_mean[iK][iphi][i] = new double * [qnpts];
                    q_long_diff_mean[iK][iphi][i] = new double * [qnpts];
                    correl_3d_Kphi_diff_num[iK][iphi][i] = new double * [qnpts];
                    correl_3d_Kphi_diff_num_count[iK][iphi][i] = new double * [qnpts];
                    correl_3d_Kphi_diff_denorm[iK][iphi][i] = new double * [qnpts];
                    for(int j = 0; j < qnpts; j++)
                    {
                        q_out_diff_mean[iK][iphi][i][j] = new double [qnpts];
                        q_side_diff_mean[iK][iphi][i][j] = new double [qnpts];
                        q_long_diff_mean[iK][iphi][i][j] = new double [qnpts];
                        correl_3d_Kphi_diff_num[iK][iphi][i][j] = new double [qnpts];
                        correl_3d_Kphi_diff_num_count[iK][iphi][i][j] = new double [qnpts];
                        correl_3d_Kphi_diff_denorm[iK][iphi][i][j] = new double [qnpts];
                        for(int k = 0; k < qnpts; k++)
                        {
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

HBT_correlation::~HBT_correlation()
{
    delete [] q_out;
    delete [] q_side;
    delete [] q_long;
    if(azimuthal_flag == 0)
    {
        for(int iK = 0; iK < n_KT; iK++)
        {
            for(int i = 0; i < qnpts; i++)
            {
                for(int j = 0; j < qnpts; j++)
                {
                    delete [] q_out_mean[iK][i][j];
                    delete [] q_side_mean[iK][i][j];
                    delete [] q_long_mean[iK][i][j];
                    delete [] correl_3d_num[iK][i][j];
                    delete [] correl_3d_num_count[iK][i][j];
                    delete [] correl_3d_denorm[iK][i][j];
                }
                delete [] q_out_mean[iK][i];
                delete [] q_side_mean[iK][i];
                delete [] q_long_mean[iK][i];
                delete [] correl_3d_num[iK][i];
                delete [] correl_3d_num_count[iK][i];
                delete [] correl_3d_denorm[iK][i];
            }
            delete [] q_out_mean[iK];
            delete [] q_side_mean[iK];
            delete [] q_long_mean[iK];
            delete [] correl_3d_num[iK];
            delete [] correl_3d_num_count[iK];
            delete [] correl_3d_denorm[iK];
        }
        delete [] q_out_mean;
        delete [] q_side_mean;
        delete [] q_long_mean;
        delete [] correl_3d_num;
        delete [] correl_3d_num_count;
        delete [] correl_3d_denorm;
    }
    else
    {
        for(int iK = 0; iK < n_KT; iK++)
        {
            for(int iphi = 0; iphi < n_Kphi; iphi++)
            {
                for(int i = 0; i < qnpts; i++)
                {
                    for(int j = 0; j < qnpts; j++)
                    {
                        delete [] q_out_diff_mean[iK][iphi][i][j];
                        delete [] q_side_diff_mean[iK][iphi][i][j];
                        delete [] q_long_diff_mean[iK][iphi][i][j];
                        delete [] correl_3d_Kphi_diff_num[iK][iphi][i][j];
                        delete [] correl_3d_Kphi_diff_num_count[iK][iphi][i][j];
                        delete [] correl_3d_Kphi_diff_denorm[iK][iphi][i][j];
                    }
                    delete [] q_out_diff_mean[iK][iphi][i];
                    delete [] q_side_diff_mean[iK][iphi][i];
                    delete [] q_long_diff_mean[iK][iphi][i];
                    delete [] correl_3d_Kphi_diff_num[iK][iphi][i];
                    delete [] correl_3d_Kphi_diff_num_count[iK][iphi][i];
                    delete [] correl_3d_Kphi_diff_denorm[iK][iphi][i];
                }
                delete [] q_out_diff_mean[iK][iphi];
                delete [] q_side_diff_mean[iK][iphi];
                delete [] q_long_diff_mean[iK][iphi];
                delete [] correl_3d_Kphi_diff_num[iK][iphi];
                delete [] correl_3d_Kphi_diff_num_count[iK][iphi];
                delete [] correl_3d_Kphi_diff_denorm[iK][iphi];
            }
            delete [] q_out_diff_mean[iK];
            delete [] q_side_diff_mean[iK];
            delete [] q_long_diff_mean[iK];
            delete [] correl_3d_Kphi_diff_num[iK];
            delete [] correl_3d_Kphi_diff_num_count[iK];
            delete [] correl_3d_Kphi_diff_denorm[iK];
        }
        delete [] q_out_diff_mean;
        delete [] q_side_diff_mean;
        delete [] q_long_diff_mean;
        delete [] correl_3d_Kphi_diff_num;
        delete [] correl_3d_Kphi_diff_num_count;
        delete [] correl_3d_Kphi_diff_denorm;
    }
    delete [] KT_array;
    for(int i = 0; i < n_KT; i++)
    {
        delete [] number_of_pairs_numerator_KTKphidiff[i];
        delete [] number_of_pairs_denormenator_KTKphidiff[i];
    }
    delete [] number_of_pairs_numerator_KTKphidiff;
    delete [] number_of_pairs_denormenator_KTKphidiff;
    delete [] number_of_pairs_numerator_KTdiff;
    delete [] number_of_pairs_denormenator_KTdiff;
}

void HBT_correlation::calculate_HBT_correlation_function()
{
    int event_id = 0;
    int buffer_size = particle_list->get_event_buffer_size();
    while(!particle_list->end_of_file())
    {
        cout << "Reading event: " << event_id+1 << "-" << event_id + buffer_size << " ... " << flush;
        particle_list->read_in_particle_samples();
        particle_list->read_in_particle_samples_mixed_event();
        cout << " processing ..." << flush;
        int nev = particle_list->get_number_of_events();
        int n_skip_ev = (int)(nev/number_of_oversample_events);
        // first pairs from the same event
        for(int iev = 0; iev < n_skip_ev; iev++)
        {
            int *event_list = new int [number_of_oversample_events];
            for(int isample = 0; isample < number_of_oversample_events; isample++)
                event_list[isample] = iev + isample*n_skip_ev;
            combine_and_bin_particle_pairs(event_list);
            delete [] event_list;
        }

        // then pairs from mixed events
        int nev_mixed_event = particle_list->get_number_of_mixed_events();
        for(int iev = 0; iev < nev_mixed_event; iev++)
        {
            int *mixed_event_list = new int [number_of_mixed_events];
            int count = 0;
            while(count < number_of_mixed_events)
            {
                int mixed_event_id = rand() % nev;
                mixed_event_list[count] = mixed_event_id;
                count++;
            }
            combine_and_bin_particle_pairs_mixed_events(iev, mixed_event_list);
            event_id++;
            delete [] mixed_event_list;
        }
        cout << " done!" << endl;
    }
    if(azimuthal_flag == 0)
        output_correlation_function();
    else
        output_correlation_function_Kphi_differential();
}

void HBT_correlation::combine_and_bin_particle_pairs(int* event_list)
{
    double hbarC_inv = 1./hbarC;
    int number_of_particles = 0;
    for(int i = 0; i < number_of_oversample_events; i++)
        number_of_particles += particle_list->get_number_of_particles(event_list[i]);
    particle_info* temp_particle_list = new particle_info [number_of_particles];  // local cache
    long int idx = 0;
    for(int j = 0; j < number_of_oversample_events; j++)
    {
        int event_id = event_list[j];
        int event_number_particle = particle_list->get_number_of_particles(event_list[j]);
        for(int i = 0; i < event_number_particle; i++)
        {
            temp_particle_list[idx].px = particle_list->get_particle(event_id, i).px;
            temp_particle_list[idx].py = particle_list->get_particle(event_id, i).py;
            temp_particle_list[idx].pz = particle_list->get_particle(event_id, i).pz;
            temp_particle_list[idx].E = particle_list->get_particle(event_id, i).E;
            temp_particle_list[idx].mass = particle_list->get_particle(event_id, i).mass;
            temp_particle_list[idx].x = particle_list->get_particle(event_id, i).x;
            temp_particle_list[idx].y = particle_list->get_particle(event_id, i).y;
            temp_particle_list[idx].z = particle_list->get_particle(event_id, i).z;
            temp_particle_list[idx].t = particle_list->get_particle(event_id, i).t;
            idx++;
        }
    }

    // nested pair loop
    for(int i = 0; i < number_of_particles; i++)
    {
        for(int j = i+1; j < number_of_particles; j++)
        {
            double K_z = 0.5*(temp_particle_list[i].pz + temp_particle_list[j].pz);
            double K_E = 0.5*(temp_particle_list[i].E + temp_particle_list[j].E);
            double K_rap = 0.5*log((K_E + K_z)/(K_E - K_z));
            if(K_rap > Krap_min && K_rap < Krap_max)  // check rapidity cut
            {
                double K_x = 0.5*(temp_particle_list[i].px + temp_particle_list[j].px);
                double K_y = 0.5*(temp_particle_list[i].py + temp_particle_list[j].py);
                double K_perp = sqrt(K_x*K_x + K_y*K_y);
                if(K_perp > KT_min && K_perp < KT_max)  // check K_T cut
                {
                    int Kperp_idx = (int)((K_perp - KT_min)/dKT);
                    double local_K_phi;
                    int Kphi_idx;
                    if(azimuthal_flag == 0)
                    {
                        if(number_of_pairs_numerator_KTdiff[Kperp_idx] > needed_number_of_pairs)
                            continue;
                        number_of_pairs_numerator_KTdiff[Kperp_idx]++;
                    }
                    else
                    {
                        local_K_phi = atan2(K_y, K_x);
                        Kphi_idx = (int)((local_K_phi)/dKphi);
                        if(number_of_pairs_numerator_KTKphidiff[Kperp_idx][Kphi_idx] > needed_number_of_pairs)
                            continue;
                        number_of_pairs_numerator_KTKphidiff[Kperp_idx][Kphi_idx]++;
                    }

                    double cos_K_phi = K_x/K_perp;
                    double sin_K_phi = K_y/K_perp;

                    double q_z = temp_particle_list[i].pz - temp_particle_list[j].pz;
                    double local_q_long = q_z;
                    if(local_q_long < (q_min - delta_q/2.) || local_q_long >= (q_max + delta_q/2.)) continue;

                    double q_x = temp_particle_list[i].px - temp_particle_list[j].px;
                    double q_y = temp_particle_list[i].py - temp_particle_list[j].py;
                    double local_q_out = q_x*cos_K_phi + q_y*sin_K_phi;
                    if(local_q_out < (q_min - delta_q/2.) || local_q_out >= (q_max + delta_q/2.)) continue;
                    double local_q_side = q_y*cos_K_phi - q_x*sin_K_phi;
                    if(local_q_side < (q_min - delta_q/2.) || local_q_side >= (q_max + delta_q/2.)) continue;
                    
                    int qout_idx = (int)((local_q_out - (q_min - delta_q/2.))/delta_q);
                    int qside_idx = (int)((local_q_side - (q_min - delta_q/2.))/delta_q);
                    int qlong_idx = (int)((local_q_long - (q_min - delta_q/2.))/delta_q);

                    double q_E = temp_particle_list[i].E - temp_particle_list[j].E;

                    double t_diff = temp_particle_list[i].t - temp_particle_list[j].t;
                    double x_diff = temp_particle_list[i].x - temp_particle_list[j].x;
                    double y_diff = temp_particle_list[i].y - temp_particle_list[j].y;
                    double z_diff = temp_particle_list[i].z - temp_particle_list[j].z;

                    double cos_qx = cos((q_E*t_diff - q_x*x_diff - q_y*y_diff - q_z*z_diff)*hbarC_inv);

                    // bin results
                    if(azimuthal_flag == 0)
                    {
                        correl_3d_num_count[Kperp_idx][qout_idx][qside_idx][qlong_idx]++;
                        q_out_mean[Kperp_idx][qout_idx][qside_idx][qlong_idx] += local_q_out;
                        q_side_mean[Kperp_idx][qout_idx][qside_idx][qlong_idx] += local_q_side;
                        q_long_mean[Kperp_idx][qout_idx][qside_idx][qlong_idx] += local_q_long;
                        correl_3d_num[Kperp_idx][qout_idx][qside_idx][qlong_idx] += cos_qx;
                    }
                    else
                    {
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
    delete [] temp_particle_list;
}

void HBT_correlation::combine_and_bin_particle_pairs_mixed_events(int event_id1, int* mixed_event_list)
{
    int number_of_particles_1 = particle_list->get_number_of_particles(event_id1);
    particle_info* temp_particle_list_1 = new particle_info [number_of_particles_1];
    for(int i = 0; i < number_of_particles_1; i++)
    {
        temp_particle_list_1[i].px = particle_list->get_particle(event_id1, i).px;
        temp_particle_list_1[i].py = particle_list->get_particle(event_id1, i).py;
        temp_particle_list_1[i].pz = particle_list->get_particle(event_id1, i).pz;
        temp_particle_list_1[i].E = particle_list->get_particle(event_id1, i).E;
        temp_particle_list_1[i].mass = particle_list->get_particle(event_id1, i).mass;
        temp_particle_list_1[i].x = particle_list->get_particle(event_id1, i).x;
        temp_particle_list_1[i].y = particle_list->get_particle(event_id1, i).y;
        temp_particle_list_1[i].z = particle_list->get_particle(event_id1, i).z;
        temp_particle_list_1[i].t = particle_list->get_particle(event_id1, i).t;
    }

    // prepare for the mixed events
    int number_of_particles_2 = 0;
    for(int iev = 0; iev < number_of_mixed_events; iev++)
        number_of_particles_2 += particle_list->get_number_of_particles_mixed_event(mixed_event_list[iev]);
    particle_info* temp_particle_list_2 = new particle_info [number_of_particles_2];
    long int idx = 0;
    for(int iev = 0; iev < number_of_mixed_events; iev++)
    {
        // introduce a random rotation for the mixed event
        double random_rotation = drand48()*2*M_PI;
        double cos_phi = cos(random_rotation);
        double sin_phi = sin(random_rotation);
        int mixed_event_id = mixed_event_list[iev];
        int event_number_particle = particle_list->get_number_of_particles_mixed_event(mixed_event_id);
        for(int i = 0; i < event_number_particle; i++)
        {
            double px_temp = particle_list->get_particle_from_mixed_event(mixed_event_id, i).px;
            double py_temp = particle_list->get_particle_from_mixed_event(mixed_event_id, i).py;
            temp_particle_list_2[idx].px = px_temp*cos_phi - py_temp*sin_phi;
            temp_particle_list_2[idx].py = px_temp*sin_phi + py_temp*cos_phi;
            double x_temp = particle_list->get_particle_from_mixed_event(mixed_event_id, i).x;
            double y_temp = particle_list->get_particle_from_mixed_event(mixed_event_id, i).y;
            temp_particle_list_2[idx].x = x_temp*cos_phi - y_temp*sin_phi;
            temp_particle_list_2[idx].y = x_temp*sin_phi + y_temp*cos_phi;

            temp_particle_list_2[idx].pz = particle_list->get_particle_from_mixed_event(mixed_event_id, i).pz;
            temp_particle_list_2[idx].E = particle_list->get_particle_from_mixed_event(mixed_event_id, i).E;
            temp_particle_list_2[idx].mass = particle_list->get_particle_from_mixed_event(mixed_event_id, i).mass;
            temp_particle_list_2[idx].z = particle_list->get_particle_from_mixed_event(mixed_event_id, i).z;
            temp_particle_list_2[idx].t = particle_list->get_particle_from_mixed_event(mixed_event_id, i).t;
            idx++;
        }
    }

    // nested pair loop
    for(int i = 0; i < number_of_particles_1; i++)
    {
        for(int j = 0; j < number_of_particles_2; j++)
        {
            double K_z = 0.5*(temp_particle_list_1[i].pz + temp_particle_list_2[j].pz);
            double K_E = 0.5*(temp_particle_list_1[i].E + temp_particle_list_2[j].E);
            double K_rap = 0.5*log((K_E + K_z)/(K_E - K_z));
            if(K_rap > Krap_min && K_rap < Krap_max)
            {
                double K_x = 0.5*(temp_particle_list_1[i].px + temp_particle_list_2[j].px);
                double K_y = 0.5*(temp_particle_list_1[i].py + temp_particle_list_2[j].py);
                double K_perp = sqrt(K_x*K_x + K_y*K_y);
                if(K_perp > KT_min && K_perp < KT_max)
                {
                    int Kperp_idx = (int)((K_perp - KT_min)/dKT);
                    double local_K_phi;
                    int Kphi_idx;
                    if(azimuthal_flag == 0)
                    {
                        if(number_of_pairs_denormenator_KTdiff[Kperp_idx] > needed_number_of_pairs)
                            continue;
                        number_of_pairs_denormenator_KTdiff[Kperp_idx]++;
                    }
                    else
                    {
                        local_K_phi = atan2(K_y, K_x);
                        Kphi_idx = (int)((local_K_phi)/dKphi);
                        if(number_of_pairs_denormenator_KTKphidiff[Kperp_idx][Kphi_idx] > needed_number_of_pairs)
                            continue;
                        number_of_pairs_denormenator_KTKphidiff[Kperp_idx][Kphi_idx]++;
                    }
                    
                    double cos_K_phi = K_x/K_perp;
                    double sin_K_phi = K_y/K_perp;
                    
                    double q_z = temp_particle_list_1[i].pz - temp_particle_list_2[j].pz;
                    double local_q_long = q_z;
                    if(local_q_long < (q_min - delta_q/2.) || local_q_long >= (q_max + delta_q/2.)) continue;
                    
                    double q_x = temp_particle_list_1[i].px - temp_particle_list_2[j].px;
                    double q_y = temp_particle_list_1[i].py - temp_particle_list_2[j].py;
                    double local_q_out = q_x*cos_K_phi + q_y*sin_K_phi;
                    if(local_q_out < (q_min - delta_q/2.) || local_q_out >= (q_max + delta_q/2.)) continue;
                    double local_q_side = q_y*cos_K_phi - q_x*sin_K_phi;
                    if(local_q_side < (q_min - delta_q/2.) || local_q_side >= (q_max + delta_q/2.)) continue;
                    
                    int qout_idx = (int)((local_q_out - (q_min - delta_q/2.))/delta_q); int qside_idx = (int)((local_q_side - (q_min - delta_q/2.))/delta_q);
                    int qlong_idx = (int)((local_q_long - (q_min - delta_q/2.))/delta_q);
// bin results
                    if(azimuthal_flag == 0)
                        correl_3d_denorm[Kperp_idx][qout_idx][qside_idx][qlong_idx] += 1.0;
                    else
                        correl_3d_Kphi_diff_denorm[Kperp_idx][Kphi_idx][qout_idx][qside_idx][qlong_idx] += 1.0;
                }
            }
        }
    }
    delete [] temp_particle_list_1;
    delete [] temp_particle_list_2;
}

void HBT_correlation::bin_into_correlation_function(int type, int num_pair, particle_pair* pairlist)
{
    if(type == 0)  // pairs from the same event
        number_pairs_num += num_pair;
    else   // pairs from mixed events
        number_pairs_denorm += num_pair;

    for(int ipair = 0; ipair < num_pair; ipair++)
    {
        double rap = pairlist[ipair].K_rap;
        if(rap > Krap_min && rap < Krap_max)
        {
            double KT_local = pairlist[ipair].K_perp;
            if(KT_local > KT_min && KT_local < KT_max)
            {
                int Kperp_idx = (int)((KT_local - KT_min)/dKT);
                double q_out_local = pairlist[ipair].q_out;
                if(q_out_local < (q_min - delta_q/2.) || q_out_local >= (q_max + delta_q/2.)) continue;
                double q_side_local = pairlist[ipair].q_side;
                if(q_side_local < (q_min - delta_q/2.) || q_side_local >= (q_max + delta_q/2.)) continue;
                double q_long_local = pairlist[ipair].q_long;
                if(q_long_local < (q_min - delta_q/2.) || q_long_local >= (q_max + delta_q/2.)) continue;

                int qout_idx = (int)((pairlist[ipair].q_out - (q_min - delta_q/2.))/delta_q);
                int qside_idx = (int)((pairlist[ipair].q_side - (q_min - delta_q/2.))/delta_q);
                int qlong_idx = (int)((pairlist[ipair].q_long - (q_min - delta_q/2.))/delta_q);

                if(azimuthal_flag == 0)
                {
                    if(type == 0)  // pairs from same event
                    {
                        correl_3d_num_count[Kperp_idx][qout_idx][qside_idx][qlong_idx]++;
                        q_out_mean[Kperp_idx][qout_idx][qside_idx][qlong_idx] += pairlist[ipair].q_out;
                        q_side_mean[Kperp_idx][qout_idx][qside_idx][qlong_idx] += pairlist[ipair].q_side;
                        q_long_mean[Kperp_idx][qout_idx][qside_idx][qlong_idx] += pairlist[ipair].q_long;
                        correl_3d_num[Kperp_idx][qout_idx][qside_idx][qlong_idx] += pairlist[ipair].cos_qx;
                    }
                    else   // pairs from mixed events
                        correl_3d_denorm[Kperp_idx][qout_idx][qside_idx][qlong_idx] += 1.0;
                }
                else
                {
                    int Kphi_idx = (int)((pairlist[ipair].K_phi)/dKphi);
                    if(type == 0)  // pairs from same event
                    {
                        correl_3d_Kphi_diff_num_count[Kperp_idx][Kphi_idx][qout_idx][qside_idx][qlong_idx]++;
                        q_out_diff_mean[Kperp_idx][Kphi_idx][qout_idx][qside_idx][qlong_idx] += pairlist[ipair].q_out;
                        q_side_diff_mean[Kperp_idx][Kphi_idx][qout_idx][qside_idx][qlong_idx] += pairlist[ipair].q_side;
                        q_long_diff_mean[Kperp_idx][Kphi_idx][qout_idx][qside_idx][qlong_idx] += pairlist[ipair].q_long;
                        correl_3d_Kphi_diff_num[Kperp_idx][Kphi_idx][qout_idx][qside_idx][qlong_idx] += pairlist[ipair].cos_qx;
                    }
                    else     // pairs from mixed events
                        correl_3d_Kphi_diff_denorm[Kperp_idx][Kphi_idx][qout_idx][qside_idx][qlong_idx] += 1.0;
                }
            }
        }
    }
}

void HBT_correlation::output_correlation_function()
{
    //double npair_ratio = (double)number_pairs_num/(double)number_pairs_denorm;
    double npair_ratio = 1.0;
    for(int iK = 0; iK < n_KT - 1; iK++)
    {
        ostringstream filename;
        filename << path << "/HBT_correlation_function_KT_" << KT_array[iK] << "_" << KT_array[iK+1] << ".dat";
        ofstream output(filename.str().c_str());
        for(int iqlong = 0; iqlong < qnpts; iqlong++)
        {
            for(int iqout = 0; iqout < qnpts; iqout++)
            {
                for(int iqside = 0; iqside < qnpts; iqside++)
                {
                    int npart_num = correl_3d_num_count[iK][iqout][iqside][iqlong];
                    int npart_denorm = correl_3d_denorm[iK][iqout][iqside][iqlong];
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
                        q_out_local = q_out_mean[iK][iqout][iqside][iqlong]/npart_num;
                        q_side_local = q_side_mean[iK][iqout][iqside][iqlong]/npart_num;
                        q_long_local = q_long_mean[iK][iqout][iqside][iqlong]/npart_num;

                        correl_fun_num = correl_3d_num[iK][iqout][iqside][iqlong];
                        correl_fun_denorm = correl_3d_denorm[iK][iqout][iqside][iqlong];
                        correl_fun_val = correl_fun_num/(correl_fun_denorm*npair_ratio);
                    }

                    output << scientific << setw(18) << setprecision(8) 
                           << q_out_local << "    " << q_side_local << "    " << q_long_local << "    "
                           << npart_num << "    " << correl_fun_num << "    "  << correl_fun_denorm << "    "
                           << correl_fun_val << "    " << 0.0 << endl;
                }
            }
        }
        output.close();
    }
}

void HBT_correlation::output_correlation_function_Kphi_differential()
{
    //double npair_ratio = (double)number_pairs_num/(double)number_pairs_denorm;
    double npair_ratio = 1.0;
    for(int iK = 0; iK < n_KT - 1; iK++)
    {
        for(int iKphi = 0; iKphi < n_Kphi; iKphi++)
        {
            ostringstream filename;
            filename << path << "/HBT_correlation_function_KT_" << KT_array[iK] << "_" << KT_array[iK+1] << "_Kphi_" << Kphi_array[iKphi] << ".dat";
            ofstream output(filename.str().c_str());
            for(int iqlong = 0; iqlong < qnpts; iqlong++)
            {
                for(int iqout = 0; iqout < qnpts; iqout++)
                {
                    for(int iqside = 0; iqside < qnpts; iqside++)
                    {
                        int npart_num = correl_3d_Kphi_diff_num_count[iK][iKphi][iqout][iqside][iqlong];
                        int npart_denorm = correl_3d_Kphi_diff_denorm[iK][iKphi][iqout][iqside][iqlong];
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
                            q_out_local = q_out_diff_mean[iK][iKphi][iqout][iqside][iqlong]/npart_num;
                            q_side_local = q_side_diff_mean[iK][iKphi][iqout][iqside][iqlong]/npart_num;
                            q_long_local = q_long_diff_mean[iK][iKphi][iqout][iqside][iqlong]/npart_num;

                            correl_fun_num = correl_3d_Kphi_diff_num[iK][iKphi][iqout][iqside][iqlong];
                            correl_fun_denorm = correl_3d_Kphi_diff_denorm[iK][iKphi][iqout][iqside][iqlong];
                            correl_fun_val = correl_fun_num/correl_fun_denorm/npair_ratio;
                        }

                        output << scientific << setw(18) << setprecision(8) 
                               << q_out_local << "    " << q_side_local << "    "
                               << q_long_local << "    "
                               << npart_num << "    "
                               << correl_fun_num << "    " << correl_fun_denorm << "    "
                               << correl_fun_val << "    " << 0.0 << endl;
                    }
                }
            }
            output.close();
        }
    }
}
