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
    for(int i = 0; i < n_KT; i++)
        KT_array[i] = KT_min + i*dKT;

    Kphi_array = new double [n_Kphi];
    for(int i = 0; i < n_Kphi; i++)
        Kphi_array[i] = 0.0 + i*dKphi;

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
        correl_3d_num_err = new double *** [n_KT];
        correl_3d_num_count = new double *** [n_KT];
        correl_3d_denorm = new double *** [n_KT];
        for(int iK = 0; iK < n_KT; iK++)
        {
            q_out_mean[iK] = new double ** [qnpts];
            q_side_mean[iK] = new double ** [qnpts];
            q_long_mean[iK] = new double ** [qnpts];
            correl_3d_num[iK] = new double ** [qnpts];
            correl_3d_num_err[iK] = new double ** [qnpts];
            correl_3d_num_count[iK] = new double ** [qnpts];
            correl_3d_denorm[iK] = new double ** [qnpts];
            for(int i = 0; i < qnpts; i++)
            {
                q_out_mean[iK][i] = new double * [qnpts];
                q_side_mean[iK][i] = new double * [qnpts];
                q_long_mean[iK][i] = new double * [qnpts];
                correl_3d_num[iK][i] = new double * [qnpts];
                correl_3d_num_err[iK][i] = new double * [qnpts];
                correl_3d_num_count[iK][i] = new double * [qnpts];
                correl_3d_denorm[iK][i] = new double * [qnpts];
                for(int j = 0; j < qnpts; j++)
                {
                    q_out_mean[iK][i][j] = new double [qnpts];
                    q_side_mean[iK][i][j] = new double [qnpts];
                    q_long_mean[iK][i][j] = new double [qnpts];
                    correl_3d_num[iK][i][j] = new double [qnpts];
                    correl_3d_num_err[iK][i][j] = new double [qnpts];
                    correl_3d_num_count[iK][i][j] = new double [qnpts];
                    correl_3d_denorm[iK][i][j] = new double [qnpts];
                    for(int k = 0; k < qnpts; k++)
                    {
                        q_out_mean[iK][i][j][k] = 0.0;
                        q_side_mean[iK][i][j][k] = 0.0;
                        q_long_mean[iK][i][j][k] = 0.0;
                        correl_3d_num[iK][i][j][k] = 0.0;
                        correl_3d_num_err[iK][i][j][k] = 0.0;
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
        correl_3d_Kphi_diff_num_err = new double **** [n_KT];
        correl_3d_Kphi_diff_num_count = new double **** [n_KT];
        correl_3d_Kphi_diff_denorm = new double **** [n_KT];
        for(int iK = 0; iK < n_KT; iK++)
        {
            q_out_diff_mean[iK] = new double *** [n_Kphi];
            q_side_diff_mean[iK] = new double *** [n_Kphi];
            q_long_diff_mean[iK] = new double *** [n_Kphi];
            correl_3d_Kphi_diff_num[iK] = new double *** [n_Kphi];
            correl_3d_Kphi_diff_num_err[iK] = new double *** [n_Kphi];
            correl_3d_Kphi_diff_num_count[iK] = new double *** [n_Kphi];
            correl_3d_Kphi_diff_denorm[iK] = new double *** [n_Kphi];
            for(int iphi = 0; iphi < n_Kphi; iphi++)
            {
                q_out_diff_mean[iK][iphi] = new double ** [qnpts];
                q_side_diff_mean[iK][iphi] = new double ** [qnpts];
                q_long_diff_mean[iK][iphi] = new double ** [qnpts];
                correl_3d_Kphi_diff_num[iK][iphi] = new double ** [qnpts];
                correl_3d_Kphi_diff_num_err[iK][iphi] = new double ** [qnpts];
                correl_3d_Kphi_diff_num_count[iK][iphi] = new double ** [qnpts];
                correl_3d_Kphi_diff_denorm[iK][iphi] = new double ** [qnpts];
                for(int i = 0; i < qnpts; i++)
                {
                    q_out_diff_mean[iK][iphi][i] = new double * [qnpts];
                    q_side_diff_mean[iK][iphi][i] = new double * [qnpts];
                    q_long_diff_mean[iK][iphi][i] = new double * [qnpts];
                    correl_3d_Kphi_diff_num[iK][iphi][i] = new double * [qnpts];
                    correl_3d_Kphi_diff_num_err[iK][iphi][i] = new double * [qnpts];
                    correl_3d_Kphi_diff_num_count[iK][iphi][i] = new double * [qnpts];
                    correl_3d_Kphi_diff_denorm[iK][iphi][i] = new double * [qnpts];
                    for(int j = 0; j < qnpts; j++)
                    {
                        q_out_diff_mean[iK][iphi][i][j] = new double [qnpts];
                        q_side_diff_mean[iK][iphi][i][j] = new double [qnpts];
                        q_long_diff_mean[iK][iphi][i][j] = new double [qnpts];
                        correl_3d_Kphi_diff_num[iK][iphi][i][j] = new double [qnpts];
                        correl_3d_Kphi_diff_num_err[iK][iphi][i][j] = new double [qnpts];
                        correl_3d_Kphi_diff_num_count[iK][iphi][i][j] = new double [qnpts];
                        correl_3d_Kphi_diff_denorm[iK][iphi][i][j] = new double [qnpts];
                        for(int k = 0; k < qnpts; k++)
                        {
                            q_out_diff_mean[iK][iphi][i][j][k] = 0.0;
                            q_side_diff_mean[iK][iphi][i][j][k] = 0.0;
                            q_long_diff_mean[iK][iphi][i][j][k] = 0.0;
                            correl_3d_Kphi_diff_num[iK][iphi][i][j][k] = 0.0;
                            correl_3d_Kphi_diff_num_err[iK][iphi][i][j][k] = 0.0;
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
                    delete [] correl_3d_num_err[iK][i][j];
                    delete [] correl_3d_num_count[iK][i][j];
                    delete [] correl_3d_denorm[iK][i][j];
                }
                delete [] q_out_mean[iK][i];
                delete [] q_side_mean[iK][i];
                delete [] q_long_mean[iK][i];
                delete [] correl_3d_num[iK][i];
                delete [] correl_3d_num_err[iK][i];
                delete [] correl_3d_num_count[iK][i];
                delete [] correl_3d_denorm[iK][i];
            }
            delete [] q_out_mean[iK];
            delete [] q_side_mean[iK];
            delete [] q_long_mean[iK];
            delete [] correl_3d_num[iK];
            delete [] correl_3d_num_err[iK];
            delete [] correl_3d_num_count[iK];
            delete [] correl_3d_denorm[iK];
        }
        delete [] q_out_mean;
        delete [] q_side_mean;
        delete [] q_long_mean;
        delete [] correl_3d_num;
        delete [] correl_3d_num_err;
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
                        delete [] correl_3d_Kphi_diff_num_err[iK][iphi][i][j];
                        delete [] correl_3d_Kphi_diff_num_count[iK][iphi][i][j];
                        delete [] correl_3d_Kphi_diff_denorm[iK][iphi][i][j];
                    }
                    delete [] q_out_diff_mean[iK][iphi][i];
                    delete [] q_side_diff_mean[iK][iphi][i];
                    delete [] q_long_diff_mean[iK][iphi][i];
                    delete [] correl_3d_Kphi_diff_num[iK][iphi][i];
                    delete [] correl_3d_Kphi_diff_num_err[iK][iphi][i];
                    delete [] correl_3d_Kphi_diff_num_count[iK][iphi][i];
                    delete [] correl_3d_Kphi_diff_denorm[iK][iphi][i];
                }
                delete [] q_out_diff_mean[iK][iphi];
                delete [] q_side_diff_mean[iK][iphi];
                delete [] q_long_diff_mean[iK][iphi];
                delete [] correl_3d_Kphi_diff_num[iK][iphi];
                delete [] correl_3d_Kphi_diff_num_err[iK][iphi];
                delete [] correl_3d_Kphi_diff_num_count[iK][iphi];
                delete [] correl_3d_Kphi_diff_denorm[iK][iphi];
            }
            delete [] q_out_diff_mean[iK];
            delete [] q_side_diff_mean[iK];
            delete [] q_long_diff_mean[iK];
            delete [] correl_3d_Kphi_diff_num[iK];
            delete [] correl_3d_Kphi_diff_num_err[iK];
            delete [] correl_3d_Kphi_diff_num_count[iK];
            delete [] correl_3d_Kphi_diff_denorm[iK];
        }
        delete [] q_out_diff_mean;
        delete [] q_side_diff_mean;
        delete [] q_long_diff_mean;
        delete [] correl_3d_Kphi_diff_num;
        delete [] correl_3d_Kphi_diff_num_err;
        delete [] correl_3d_Kphi_diff_num_count;
        delete [] correl_3d_Kphi_diff_denorm;
    }
    delete [] KT_array;
}

void HBT_correlation::calculate_HBT_correlation_function()
{
    int event_id = 0;
    int buffer_size = particle_list->get_event_buffer_size();
    while(!particle_list->end_of_file())
    {
        cout << "Reading event: " << event_id+1 << "-" << event_id + buffer_size << " ... " << flush;
        particle_list->read_in_particle_samples();
        cout << " processing ..." << flush;
        int nev = particle_list->get_number_of_events();
        // first pairs from the same event
        for(int iev = 0; iev < nev/number_of_oversample_events; iev++)
        {
            int n_skip_ev = nev/number_of_oversample_events;
            int number_of_particles = 0;
            int *event_list = new int [number_of_oversample_events];
            for(int isample = 0; isample < number_of_oversample_events; isample++)
            {
                event_list[isample] = iev + isample*n_skip_ev;
                number_of_particles += particle_list->get_number_of_particles(iev+isample*n_skip_ev);
            }
            long int num_of_pairs = number_of_particles*(number_of_particles - 1)/2;
            particle_pair *particle_pairs_list = new particle_pair [num_of_pairs];
            combine_particle_pairs(event_list, particle_pairs_list);
            bin_into_correlation_function(0, num_of_pairs, particle_pairs_list);
            delete [] particle_pairs_list;
        }

        // then pairs from mixed events
        for(int iev = 0; iev < nev; iev++)
        {
            event_id++;
            int number_of_particles = particle_list->get_number_of_particles(iev);
            int count = 0;
            while(1)
            {
                int mixed_event_id = rand() % buffer_size;
                if(mixed_event_id != iev)
                {
                    int number_of_particles_2 = particle_list->get_number_of_particles(mixed_event_id);
                    int num_of_mixed_pairs = number_of_particles*number_of_particles_2;
                    particle_pair *particle_mixed_pairs_list = new particle_pair [num_of_mixed_pairs];
                    combine_particle_pairs_mixed_events(iev, mixed_event_id, particle_mixed_pairs_list);
                    bin_into_correlation_function(1, num_of_mixed_pairs, particle_mixed_pairs_list);
                    delete [] particle_mixed_pairs_list;
                    count++;
                }
                if(count == number_of_mixed_events) break;
            }
        }
        cout << " done!" << endl;
    }
    nevent = event_id;
    if(azimuthal_flag == 0)
        output_correlation_function();
    else
        output_correlation_function_Kphi_differential();
}

void HBT_correlation::combine_particle_pairs(int* event_list, particle_pair* list)
{
    int number_of_particles = 0;
    for(int i = 0; i < number_of_oversample_events; i++)
        number_of_particles += particle_list->get_number_of_particles(event_list[i]);
    particle_info* temp_particle_list = new particle_info [number_of_particles];
    int idx = 0;
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
    int pair_id = 0;
    for(int i = 0; i < number_of_particles; i++)
    {
        for(int j = i+1; j < number_of_particles; j++)
        {
            double K_x = 0.5*(temp_particle_list[i].px + temp_particle_list[j].px);
            double K_y = 0.5*(temp_particle_list[i].py + temp_particle_list[j].py);
            double K_z = 0.5*(temp_particle_list[i].pz + temp_particle_list[j].pz);
            double K_E = 0.5*(temp_particle_list[i].E + temp_particle_list[j].E);

            double q_x = temp_particle_list[i].px - temp_particle_list[j].px;
            double q_y = temp_particle_list[i].py - temp_particle_list[j].py;
            double q_z = temp_particle_list[i].pz - temp_particle_list[j].pz;
            double q_E = temp_particle_list[i].E - temp_particle_list[j].E;

            double t_diff = temp_particle_list[i].t - temp_particle_list[j].t;
            double x_diff = temp_particle_list[i].x - temp_particle_list[j].x;
            double y_diff = temp_particle_list[i].y - temp_particle_list[j].y;
            double z_diff = temp_particle_list[i].z - temp_particle_list[j].z;

            double K_perp = sqrt(K_x*K_x + K_y*K_y);
            double cos_K_phi = K_x/K_perp;
            double sin_K_phi = K_y/K_perp;
            list[pair_id].K_perp = K_perp;
            list[pair_id].K_phi = atan2(K_y, K_x);
            list[pair_id].K_rap = 0.5*log((K_E + K_z)/(K_E - K_z));

            list[pair_id].q_out = q_x*cos_K_phi + q_y*sin_K_phi;
            list[pair_id].q_side = q_y*cos_K_phi - q_x*sin_K_phi;
            list[pair_id].q_long = q_z;
            list[pair_id].cos_qx = cos((q_E*t_diff - q_x*x_diff - q_y*y_diff - q_z*z_diff)/hbarC);
            pair_id++;
        }
    }
    delete [] temp_particle_list;
}

void HBT_correlation::combine_particle_pairs_mixed_events(int event_id1, int event_id2, particle_pair* list)
{
    int number_of_particles_1 = particle_list->get_number_of_particles(event_id1);
    int number_of_particles_2 = particle_list->get_number_of_particles(event_id2);
    particle_info* temp_particle_list_1 = new particle_info [number_of_particles_1];
    particle_info* temp_particle_list_2 = new particle_info [number_of_particles_2];
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
    // introduce a random rotation for the mixed event
    double random_rotation = drand48()*2*M_PI;
    double cos_phi = cos(random_rotation);
    double sin_phi = sin(random_rotation);
    for(int i = 0; i < number_of_particles_2; i++)
    {
        double px_temp = particle_list->get_particle(event_id2, i).px;
        double py_temp = particle_list->get_particle(event_id2, i).py;
        temp_particle_list_2[i].px = px_temp*cos_phi - py_temp*sin_phi;
        temp_particle_list_2[i].py = px_temp*sin_phi + py_temp*cos_phi;
        double x_temp = particle_list->get_particle(event_id2, i).x;
        double y_temp = particle_list->get_particle(event_id2, i).y;
        temp_particle_list_2[i].x = x_temp*cos_phi - y_temp*sin_phi;
        temp_particle_list_2[i].y = x_temp*sin_phi + y_temp*cos_phi;

        temp_particle_list_2[i].pz = particle_list->get_particle(event_id2, i).pz;
        temp_particle_list_2[i].E = particle_list->get_particle(event_id2, i).E;
        temp_particle_list_2[i].mass = particle_list->get_particle(event_id2, i).mass;
        temp_particle_list_2[i].z = particle_list->get_particle(event_id2, i).z;
        temp_particle_list_2[i].t = particle_list->get_particle(event_id2, i).t;
    }

    // nested pair loop
    int pair_id = 0;
    for(int i = 0; i < number_of_particles_1; i++)
    {
        for(int j = 0; j < number_of_particles_2; j++)
        {
            double K_x = 0.5*(temp_particle_list_1[i].px + temp_particle_list_2[j].px);
            double K_y = 0.5*(temp_particle_list_1[i].py + temp_particle_list_2[j].py);
            double K_z = 0.5*(temp_particle_list_1[i].pz + temp_particle_list_2[j].pz);
            double K_E = 0.5*(temp_particle_list_1[i].E + temp_particle_list_2[j].E);

            double q_x = temp_particle_list_1[i].px - temp_particle_list_2[j].px;
            double q_y = temp_particle_list_1[i].py - temp_particle_list_2[j].py;
            double q_z = temp_particle_list_1[i].pz - temp_particle_list_2[j].pz;

            double K_perp = sqrt(K_x*K_x + K_y*K_y);
            double cos_K_phi = K_x/K_perp;
            double sin_K_phi = K_y/K_perp;
            list[pair_id].K_perp = K_perp;
            list[pair_id].K_phi = atan2(K_y, K_x);
            list[pair_id].K_rap = 0.5*log((K_E + K_z)/(K_E - K_z));

            list[pair_id].q_out = q_x*cos_K_phi + q_y*sin_K_phi;
            list[pair_id].q_side = q_y*cos_K_phi - q_x*sin_K_phi;
            list[pair_id].q_long = q_z;
            list[pair_id].cos_qx = 1.0;
            pair_id++;
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
                if(q_out_local < (q_min - delta_q/2.) || q_out_local > (q_max + delta_q/2.)) continue;
                double q_side_local = pairlist[ipair].q_side;
                if(q_side_local < (q_min - delta_q/2.) || q_side_local > (q_max + delta_q/2.)) continue;
                double q_long_local = pairlist[ipair].q_long;
                if(q_long_local < (q_min - delta_q/2.) || q_long_local > (q_max + delta_q/2.)) continue;

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
                        correl_3d_num_err[Kperp_idx][qout_idx][qside_idx][qlong_idx] += pairlist[ipair].cos_qx*pairlist[ipair].cos_qx;
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
                        correl_3d_Kphi_diff_num_err[Kperp_idx][Kphi_idx][qout_idx][qside_idx][qlong_idx] += pairlist[ipair].cos_qx*pairlist[ipair].cos_qx;
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
    double npair_ratio = (double)number_pairs_num/(double)number_pairs_denorm;
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
                    double q_out_local, q_side_local, q_long_local;
                    double correl_fun_num, correl_fun_denorm, correl_fun_val;
                    double num_err, denorm_err, error;
                    if(npart_num < 2)
                    {
                        q_out_local = q_out[iqout];
                        q_side_local = q_side[iqside];
                        q_long_local = q_long[iqlong];
                        correl_fun_num = 0.0;
                        correl_fun_denorm = correl_3d_denorm[iK][iqout][iqside][iqlong];
                        correl_fun_val = 0.0;
                        error = 0.0;
                    }
                    else
                    {
                        q_out_local = q_out_mean[iK][iqout][iqside][iqlong]/npart_num;
                        q_side_local = q_side_mean[iK][iqout][iqside][iqlong]/npart_num;
                        q_long_local = q_long_mean[iK][iqout][iqside][iqlong]/npart_num;

                        correl_fun_num = correl_3d_num[iK][iqout][iqside][iqlong];
                        correl_fun_denorm = correl_3d_denorm[iK][iqout][iqside][iqlong];
                        correl_fun_val = correl_fun_num/correl_fun_denorm/npair_ratio;
                        num_err = sqrt((correl_3d_num_err[iK][iqout][iqside][iqlong]/npart_num - correl_fun_num/npart_num*correl_fun_num/npart_num)/npart_num);
                        denorm_err = 0.0;
                        error = sqrt(pow(num_err/correl_fun_denorm, 2) + pow(denorm_err*correl_fun_num/correl_fun_denorm/correl_fun_denorm, 2))/npair_ratio;
                    }

                    output << scientific << setw(18) << setprecision(8) 
                           << q_out_local << "    " << q_side_local << "    " 
                           << q_long_local << "    "
                           << correl_fun_num << "    "  << correl_fun_denorm << "    "
                           << correl_fun_val << "    "  << error 
                           << endl;
                }
            }
        }
        output.close();
    }
}

void HBT_correlation::output_correlation_function_Kphi_differential()
{
    double npair_ratio = (double)number_pairs_num/(double)number_pairs_denorm;
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
                        double q_out_local, q_side_local, q_long_local;
                        double correl_fun_num, correl_fun_denorm, correl_fun_val;
                        double num_err, denorm_err, error;
                        if(npart_num < 2)
                        {
                            q_out_local = q_out[iqout];
                            q_side_local = q_side[iqside];
                            q_long_local = q_long[iqlong];
                            correl_fun_num = 0.0;
                            correl_fun_denorm = correl_3d_Kphi_diff_denorm[iK][iKphi][iqout][iqside][iqlong];
                            correl_fun_val = 0.0;
                            error = 0.0;
                        }
                        else
                        {
                            q_out_local = q_out_diff_mean[iK][iKphi][iqout][iqside][iqlong]/npart_num;
                            q_side_local = q_side_diff_mean[iK][iKphi][iqout][iqside][iqlong]/npart_num;
                            q_long_local = q_long_diff_mean[iK][iKphi][iqout][iqside][iqlong]/npart_num;

                            correl_fun_num = correl_3d_Kphi_diff_num[iK][iKphi][iqout][iqside][iqlong];
                            correl_fun_denorm = correl_3d_Kphi_diff_denorm[iK][iKphi][iqout][iqside][iqlong];
                            correl_fun_val = correl_fun_num/correl_fun_denorm/npair_ratio;
                            num_err = sqrt((correl_3d_Kphi_diff_num[iK][iKphi][iqout][iqside][iqlong]/npart_num - correl_fun_num/npart_num*correl_fun_num/npart_num)/npart_num);
                            denorm_err = 0.0;
                            error = sqrt(pow(num_err/correl_fun_denorm, 2) + pow(correl_fun_num*denorm_err/correl_fun_denorm/correl_fun_denorm, 2))/npair_ratio;
                        }

                        output << scientific << setw(18) << setprecision(8) 
                               << q_out_local << "    " << q_side_local << "    "
                               << q_long_local << "    "
                               << correl_fun_num << "    " << correl_fun_denorm << "    "
                               << correl_fun_val<< "    " << error 
                               << endl;
                    }
                }
            }
            output.close();
        }
    }
}
