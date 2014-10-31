#include <cmath>
#include <iostream>
#include "parameters.h"
#include "HBT_correlation.h"
using namespace std;

HBT_correlation::HBT_correlation(ParameterReader* paraRdr_in, string path_in, particleSamples *particle_list_in)
{
    paraRdr = paraRdr_in;
    path = path_in;
    particle_list = particle_list_in;

    qnpts = paraRdr->getVal("qnpts");
    init_q = paraRdr->getVal("init_q");
    delta_q = paraRdr->getVal("delta_q");
    
    q_out = new double [qnpts];
    q_side = new double [qnpts];
    q_long = new double [qnpts];
    for(int i = 0; i < qnpts; i++)
    {
        q_out[i] = init_q + i*delta_q;
        q_side[i] = init_q + i*delta_q;
        q_long[i] = init_q + i*delta_q;
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

    if(azimuthal_flag == 0)
    {
        correl_3d = new double *** [n_KT];
        for(int iK = 0; iK < n_KT; iK++)
        {
            correl_3d[iK] = new double ** [qnpts];
            for(int i = 0; i < qnpts; i++)
            {
                correl_3d[iK][i] = new double * [qnpts];
                for(int j = 0; j < qnpts; j++)
                {
                    correl_3d[iK][i][j] = new double [qnpts];
                    for(int k = 0; k < qnpts; k++)
                        correl_3d[iK][i][j][k] = 0.0;
                }
            }
        }
    }
    else
    {
        correl_3d_Kphi_diff = new double **** [n_KT];
        for(int iK = 0; iK < n_KT; iK++)
        {
            correl_3d_Kphi_diff[iK] = new double *** [n_Kphi];
            for(int iphi = 0; iphi < n_Kphi; iphi++)
            {
                correl_3d_Kphi_diff[iK][iphi] = new double ** [qnpts];
                for(int i = 0; i < qnpts; i++)
                {
                    correl_3d_Kphi_diff[iK][iphi][i] = new double * [qnpts];
                    for(int j = 0; j < qnpts; j++)
                    {
                        correl_3d_Kphi_diff[iK][iphi][i][j] = new double [qnpts];
                        for(int k = 0; k < qnpts; k++)
                            correl_3d_Kphi_diff[iK][iphi][i][j][k] = 0.0;
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
                    delete [] correl_3d[iK][i][j];
                delete [] correl_3d[iK][i];
            }
            delete [] correl_3d[iK];
        }
        delete [] correl_3d;
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
                        delete [] correl_3d_Kphi_diff[iK][iphi][i][j];
                    delete [] correl_3d_Kphi_diff[iK][iphi][i];
                }
                delete [] correl_3d_Kphi_diff[iK][iphi];
            }
            delete [] correl_3d_Kphi_diff[iK];
        }
        delete [] correl_3d_Kphi_diff;
    }
    delete [] KT_array;
}

void HBT_correlation::calculate_HBT_correlation_function()
{
    int nev = particle_list->get_number_of_events();
    for(int iev = 0; iev < nev; iev++)
    {
        cout << "Processing event: " << iev << endl;
        int number_of_particles = particle_list->get_number_of_particles(iev);
        int num_of_pairs = number_of_particles*(number_of_particles - 1);

        particle_pair *particle_pairs_list = new particle_pair [num_of_pairs];

        combine_particle_pairs(iev, particle_pairs_list);
        bin_into_correlation_function(num_of_pairs, particle_pairs_list);

        delete [] particle_pairs_list;
    }

}

void HBT_correlation::combine_particle_pairs(int event_id, particle_pair* list)
{
    int number_of_particles = particle_list->get_number_of_particles(event_id);
    particle_info* temp_particle_list = new particle_info [number_of_particles];
    for(int i = 0; i < number_of_particles; i++)
    {
        temp_particle_list[i].px = particle_list->get_particle(event_id, i).px;
        temp_particle_list[i].py = particle_list->get_particle(event_id, i).py;
        temp_particle_list[i].pz = particle_list->get_particle(event_id, i).pz;
        temp_particle_list[i].E = particle_list->get_particle(event_id, i).E;
        temp_particle_list[i].mass = particle_list->get_particle(event_id, i).mass;
        temp_particle_list[i].x = particle_list->get_particle(event_id, i).x;
        temp_particle_list[i].y = particle_list->get_particle(event_id, i).y;
        temp_particle_list[i].z = particle_list->get_particle(event_id, i).z;
        temp_particle_list[i].t = particle_list->get_particle(event_id, i).t;
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
}

void HBT_correlation::bin_into_correlation_function(int num_pair, particle_pair* pairlist)
{
    for(int ipair = 0; ipair < num_pair; ipair++)
    {
        double rap = pairlist[ipair].K_rap;
        if(rap > Krap_min && rap < Krap_max)
        {
            double KT_local = pairlist[ipair].K_perp;
            if(KT_local > KT_min && KT_local < KT_max)
            {
                int Kperp_idx = (int)((KT_local - KT_min)/dKT);
                int qout_idx = abs((int)((pairlist[ipair].q_out - init_q)/delta_q));
                if(qout_idx >= qnpts) continue;
                int qside_idx = abs((int)((pairlist[ipair].q_side - init_q)/delta_q));
                if(qside_idx >= qnpts) continue;
                int qlong_idx = abs((int)((pairlist[ipair].q_long - init_q)/delta_q));
                if(qlong_idx >= qnpts) continue;
                if(azimuthal_flag == 0)
                    correl_3d[Kperp_idx][qout_idx][qside_idx][qlong_idx] += pairlist[ipair].cos_qx;
                else
                {
                    int Kphi_idx = (int)((pairlist[ipair].K_phi)/dKphi);
                    correl_3d_Kphi_diff[Kperp_idx][Kphi_idx][qout_idx][qside_idx][qlong_idx] += pairlist[ipair].cos_qx;
                }
            }

        }
    }
}
