#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
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

    number_of_mixed_events = paraRdr->getVal("number_of_mixed_events");
    number_pairs_num = 0;
    number_pairs_denorm = 0;
    if(azimuthal_flag == 0)
    {
        correl_3d_num = new double *** [n_KT];
        correl_3d_denorm = new double *** [n_KT];
        for(int iK = 0; iK < n_KT; iK++)
        {
            correl_3d_num[iK] = new double ** [qnpts];
            correl_3d_denorm[iK] = new double ** [qnpts];
            for(int i = 0; i < qnpts; i++)
            {
                correl_3d_num[iK][i] = new double * [qnpts];
                correl_3d_denorm[iK][i] = new double * [qnpts];
                for(int j = 0; j < qnpts; j++)
                {
                    correl_3d_num[iK][i][j] = new double [qnpts];
                    correl_3d_denorm[iK][i][j] = new double [qnpts];
                    for(int k = 0; k < qnpts; k++)
                    {
                        correl_3d_num[iK][i][j][k] = 0.0;
                        correl_3d_denorm[iK][i][j][k] = 0.0;
                    }
                }
            }
        }
    }
    else
    {
        correl_3d_Kphi_diff_num = new double **** [n_KT];
        correl_3d_Kphi_diff_denorm = new double **** [n_KT];
        for(int iK = 0; iK < n_KT; iK++)
        {
            correl_3d_Kphi_diff_num[iK] = new double *** [n_Kphi];
            correl_3d_Kphi_diff_denorm[iK] = new double *** [n_Kphi];
            for(int iphi = 0; iphi < n_Kphi; iphi++)
            {
                correl_3d_Kphi_diff_num[iK][iphi] = new double ** [qnpts];
                correl_3d_Kphi_diff_denorm[iK][iphi] = new double ** [qnpts];
                for(int i = 0; i < qnpts; i++)
                {
                    correl_3d_Kphi_diff_num[iK][iphi][i] = new double * [qnpts];
                    correl_3d_Kphi_diff_denorm[iK][iphi][i] = new double * [qnpts];
                    for(int j = 0; j < qnpts; j++)
                    {
                        correl_3d_Kphi_diff_num[iK][iphi][i][j] = new double [qnpts];
                        correl_3d_Kphi_diff_denorm[iK][iphi][i][j] = new double [qnpts];
                        for(int k = 0; k < qnpts; k++)
                        {
                            correl_3d_Kphi_diff_num[iK][iphi][i][j][k] = 0.0;
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
                    delete [] correl_3d_num[iK][i][j];
                    delete [] correl_3d_denorm[iK][i][j];
                }
                delete [] correl_3d_num[iK][i];
                delete [] correl_3d_denorm[iK][i];
            }
            delete [] correl_3d_num[iK];
            delete [] correl_3d_denorm[iK];
        }
        delete [] correl_3d_num;
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
                        delete [] correl_3d_Kphi_diff_num[iK][iphi][i][j];
                        delete [] correl_3d_Kphi_diff_denorm[iK][iphi][i][j];
                    }
                    delete [] correl_3d_Kphi_diff_num[iK][iphi][i];
                    delete [] correl_3d_Kphi_diff_denorm[iK][iphi][i];
                }
                delete [] correl_3d_Kphi_diff_num[iK][iphi];
                delete [] correl_3d_Kphi_diff_denorm[iK][iphi];
            }
            delete [] correl_3d_Kphi_diff_num[iK];
            delete [] correl_3d_Kphi_diff_denorm[iK];
        }
        delete [] correl_3d_Kphi_diff_num;
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
        for(int iev = 0; iev < nev; iev++)
        {
            event_id++;
            // first pairs from the same event
            int number_of_particles = particle_list->get_number_of_particles(iev);
            int num_of_pairs = number_of_particles*(number_of_particles - 1);
            particle_pair *particle_pairs_list = new particle_pair [num_of_pairs];
            combine_particle_pairs(iev, particle_pairs_list);
            bin_into_correlation_function(0, num_of_pairs, particle_pairs_list);
            delete [] particle_pairs_list;

            // then pairs from mixed events
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
    if(azimuthal_flag == 0)
        output_correlation_function();
    else
        output_correlation_function_Kphi_differential();
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
            double q_E = temp_particle_list_1[i].E - temp_particle_list_2[j].E;

            double t_diff = temp_particle_list_1[i].t - temp_particle_list_2[j].t;
            double x_diff = temp_particle_list_1[i].x - temp_particle_list_2[j].x;
            double y_diff = temp_particle_list_1[i].y - temp_particle_list_2[j].y;
            double z_diff = temp_particle_list_1[i].z - temp_particle_list_2[j].z;

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
                int qout_idx = abs((int)((pairlist[ipair].q_out - init_q)/delta_q));
                if(qout_idx >= qnpts) continue;
                int qside_idx = abs((int)((pairlist[ipair].q_side - init_q)/delta_q));
                if(qside_idx >= qnpts) continue;
                int qlong_idx = abs((int)((pairlist[ipair].q_long - init_q)/delta_q));
                if(qlong_idx >= qnpts) continue;
                if(azimuthal_flag == 0)
                {
                    if(type == 0)  // pairs from same event
                        correl_3d_num[Kperp_idx][qout_idx][qside_idx][qlong_idx] += pairlist[ipair].cos_qx;
                    else   // pairs from mixed events
                        correl_3d_denorm[Kperp_idx][qout_idx][qside_idx][qlong_idx] += 1.0;
                }
                else
                {
                    int Kphi_idx = (int)((pairlist[ipair].K_phi)/dKphi);
                    if(type == 0)  // pairs from same event
                        correl_3d_Kphi_diff_num[Kperp_idx][Kphi_idx][qout_idx][qside_idx][qlong_idx] += pairlist[ipair].cos_qx;
                    else     // pairs from mixed events
                        correl_3d_Kphi_diff_denorm[Kperp_idx][Kphi_idx][qout_idx][qside_idx][qlong_idx] += 1.0;
                }
            }
        }
    }
}

void HBT_correlation::output_correlation_function()
{
    for(int iK = 0; iK < n_KT; iK++)
    {
        ostringstream filename_num;
        ostringstream filename_denorm;
        filename_num << path << "/HBT_correlation_function_num_KT_" << KT_array[iK] << ".dat";
        filename_denorm << path << "/HBT_correlation_function_denorm_KT_" << KT_array[iK] << ".dat";
        ofstream output_num(filename_num.str().c_str());
        ofstream output_denorm(filename_denorm.str().c_str());
        for(int iqlong = 0; iqlong < qnpts; iqlong++)
        {
            for(int iqout = 0; iqout < qnpts; iqout++)
            {
                for(int iqside = 0; iqside < qnpts; iqside++)
                {
                    output_num << scientific << setw(18) << setprecision(8) 
                               << correl_3d_num[iK][iqout][iqside][iqlong] << "    ";
                    output_denorm << scientific << setw(18) << setprecision(8) 
                                  << correl_3d_denorm[iK][iqout][iqside][iqlong] << "    ";
                }
                output_num << endl;
                output_denorm << endl;
            }
            output_num << endl;
            output_denorm << endl;
        }
        output_num.close();
        output_denorm.close();
    }
}

void HBT_correlation::output_correlation_function_Kphi_differential()
{
    for(int iK = 0; iK < n_KT; iK++)
    {
        for(int iKphi = 0; iKphi < n_Kphi; iKphi++)
        {
            ostringstream filename_num;
            ostringstream filename_denorm;
            filename_num << path << "/HBT_correlation_function_num_KT_" << KT_array[iK] << "_Kphi_" << Kphi_array[iKphi] << ".dat";
            filename_denorm << path << "/HBT_correlation_function_denorm_KT_" << KT_array[iK] << "_Kphi_" << Kphi_array[iKphi] << ".dat";
            ofstream output_num(filename_num.str().c_str());
            ofstream output_denorm(filename_denorm.str().c_str());
            for(int iqlong = 0; iqlong < qnpts; iqlong++)
            {
                for(int iqout = 0; iqout < qnpts; iqout++)
                {
                    for(int iqside = 0; iqside < qnpts; iqside++)
                    {
                        output_num << scientific << setw(18) << setprecision(8) 
                                   << correl_3d_Kphi_diff_num[iK][iKphi][iqout][iqside][iqlong] << "    ";
                        output_denorm << scientific << setw(18) << setprecision(8) 
                                      << correl_3d_Kphi_diff_denorm[iK][iKphi][iqout][iqside][iqlong] << "    ";
                    }
                    output_num << endl;
                    output_denorm << endl;
                }
                output_num << endl;
                output_denorm << endl;
            }
            output_num.close();
            output_denorm.close();
        }
    }
}
