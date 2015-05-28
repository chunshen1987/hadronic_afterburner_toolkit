#include<iostream>
#include<sstream>
#include<string>
#include<fstream>
#include<cmath>
#include<iomanip>
#include<stdlib.h>
#include<vector>

#include "particleSamples.h"

using namespace std;

particleSamples::particleSamples(ParameterReader* paraRdr_in, string path_in)
{
    paraRdr = paraRdr_in;
    path = path_in;

    event_buffer_size = paraRdr->getVal("event_buffer_size");
    read_in_mode = paraRdr->getVal("read_in_mode");
    
    // read in particle Monte-Carlo number
    particle_monval = paraRdr->getVal("particle_monval");
    flag_isospin = paraRdr->getVal("distinguish_isospin");
    if(read_in_mode == 1 || read_in_mode == 3)
        get_UrQMD_id(particle_monval);

    particle_list = new vector< vector<particle_info>* >;

    ostringstream filename;
    if(read_in_mode == 0)
        filename << path << "/OSCAR.DAT";
    else if(read_in_mode == 1 || read_in_mode == 3)
        filename << path << "/particle_list.dat";

    inputfile.open(filename.str().c_str());
    // skip the header file for OSCAR
    string temp;
    if(read_in_mode == 0)
    {
        getline(inputfile, temp);
        getline(inputfile, temp);
        getline(inputfile, temp);
    }

    initialize_charged_hadron_urqmd_id_list();
}

particleSamples::~particleSamples()
{
    inputfile.close();
    for(int i = 0; i < particle_list->size(); i++)
        delete (*particle_list)[i];
    delete particle_list;
}

void particleSamples::initialize_charged_hadron_urqmd_id_list()
{
    charged_hadron_urqmd_id_list[0] = 101;   // pion
    charged_hadron_urqmd_id_list[1] = 106;   // kaon
    charged_hadron_urqmd_id_list[2] = 1;     // proton
    charged_hadron_urqmd_id_list[3] = 40;    // Sigma^+ and Sigma^-
    charged_hadron_urqmd_id_list[4] = 49;    // Xi^-
}

void particleSamples::get_UrQMD_id(int monval)
{
    // find the corresponding UrQMD id number
    if(monval == 211)  // pion^+
    {
        particle_urqmd_id = 101;
        particle_urqmd_isospin = 2;
    }
    else if(monval == -211)  // pion^-
    {
        particle_urqmd_id = 101;
        particle_urqmd_isospin = -2;
    }
    else if(monval == 321)  // Kaon^+
    {
        particle_urqmd_id = 106;
        particle_urqmd_isospin = 1;
    }
    else if(monval == -321)  // Kaon^+
    {
        particle_urqmd_id = -106;
        particle_urqmd_isospin = 1;
    }
    else if(monval == 2212)  // proton
    {
        particle_urqmd_id = 1;
        particle_urqmd_isospin = 1;
    }
    else if(monval == -2212)  // anti-proton
    {
        particle_urqmd_id = -1;
        particle_urqmd_isospin = 1;
    }
    else if(monval == 3222)  // Sigma^+
    {
        particle_urqmd_id = 40;
        particle_urqmd_isospin = 2;
    }
    else if(monval == -3222)  // anti-Sigma^+
    {
        particle_urqmd_id = -40;
        particle_urqmd_isospin = 2;
    }
    else if(monval == 3112)  // Sigma^-
    {
        particle_urqmd_id = 40;
        particle_urqmd_isospin = -2;
    }
    else if(monval == -3112)  // anti-Sigma^-
    {
        particle_urqmd_id = -40;
        particle_urqmd_isospin = -2;
    }
    else if(monval == 3312)  // Xi^-
    {
        particle_urqmd_id = 49;
        particle_urqmd_isospin = -1;
    }
    else if(monval == -3312)  // anti-Xi^-
    {
        particle_urqmd_id = -49;
        particle_urqmd_isospin = -1;
    }
    else if(monval == 3122)  // Lambda
    {
        particle_urqmd_id = 27;
        particle_urqmd_isospin = 0;
    }
    else if(monval == -3122)  // anti-Lambda
    {
        particle_urqmd_id = -27;
        particle_urqmd_isospin = 0;
    }
    else if(monval == 3334)  // Omega
    {
        particle_urqmd_id = 55;
        particle_urqmd_isospin = 0;
    }
    else if(monval == -3334)  // anti-Omega
    {
        particle_urqmd_id = -55;
        particle_urqmd_isospin = 0;
    }
    else if(monval == 9999)  // charged hadrons
    {
        particle_urqmd_id = 9999;
        particle_urqmd_isospin = 0;
    }
}

int particleSamples::read_in_particle_samples()
{
    if(read_in_mode == 0)
        read_in_particle_samples_OSCAR();
    else if(read_in_mode == 1)
        read_in_particle_samples_UrQMD();
    else if(read_in_mode == 3)
        read_in_particle_samples_Sangwook();

    return(0);
}

int particleSamples::decide_to_pick_UrQMD(int pid, int iso3, int charge)
{
    int pick_flag = 0;
    if(particle_urqmd_id == 9999)  // charged hadrons
    {
        int in_flag = 0;
        for(int i = 0; i < 5; i++)
        {
            if(pid == charged_hadron_urqmd_id_list[i])
            {
                in_flag = 1;
                break;
            }
        }
        if(in_flag == 1 && charge != 0)
            pick_flag = 1;
    }
    else
    {
        if(flag_isospin == 0)
        {
            if(pid == particle_urqmd_id 
               && abs(iso3) == abs(particle_urqmd_isospin))
                pick_flag = 1;
        }
        else
        {
            if(pid == particle_urqmd_id && iso3 == particle_urqmd_isospin)
                pick_flag = 1;
        }

    }
    return(pick_flag);
}

int particleSamples::read_in_particle_samples_OSCAR()
{
    // clean out the previous record
    for(int i = 0; i < particle_list->size(); i++)
        (*particle_list)[i]->clear();
    particle_list->clear();
    
    string temp_string;
    int event_id, n_particle, dummy;
    int ievent;
    int temp_monval;
    for(ievent = 0; ievent < event_buffer_size; ievent++)
    {
        getline(inputfile, temp_string);
        stringstream temp1(temp_string);
        temp1 >> event_id >> n_particle;
        if(!inputfile.eof())
        {
            particle_list->push_back(new vector<particle_info> );
            int idx = ievent;

            int pick_flag = 0;
            for(int ipart = 0; ipart < n_particle; ipart++)
            {
                getline(inputfile, temp_string);
                stringstream temp2(temp_string);
                temp2 >> dummy >> temp_monval;
                if(flag_isospin == 0)
                {
                    if(abs(temp_monval) == particle_monval)
                        pick_flag = 1;
                    else
                        pick_flag = 0;
                }
                else
                {
                    if(temp_monval == particle_monval)
                        pick_flag = 1;
                    else
                        pick_flag = 0;
                }
                if(pick_flag == 1)
                {
                     particle_info *temp_particle_info = new particle_info;
                     temp2 >> temp_particle_info->px 
                           >> temp_particle_info->py
                           >> temp_particle_info->pz 
                           >> temp_particle_info->E
                           >> temp_particle_info->mass 
                           >> temp_particle_info->x 
                           >> temp_particle_info->y
                           >> temp_particle_info->z 
                           >> temp_particle_info->t;
                     (*particle_list)[idx]->push_back(*temp_particle_info);
                }
            }
        }
        else
            break;
    }
    return(0);
}

int particleSamples::read_in_particle_samples_UrQMD()
{
    // clean out the previous record
    for(int i = 0; i < particle_list->size(); i++)
        (*particle_list)[i]->clear();
    particle_list->clear();

    string temp_string;
    int n_particle;
    double dummy;
    int ievent;
    int urqmd_pid, urqmd_iso3, urqmd_charge;
    double temp_mass;
    for(ievent = 0; ievent < event_buffer_size; ievent++)
    {
        getline(inputfile, temp_string);
        if(!inputfile.eof())
        {
            particle_list->push_back(new vector<particle_info> );

            // first skip the header
            for(int i = 0; i < 13; i++)
                getline(inputfile, temp_string);
            // then get number of particles within the event
            getline(inputfile, temp_string);
            stringstream temp1(temp_string);
            temp1 >> n_particle;
            getline(inputfile, temp_string);  // then get one useless line

            int idx = ievent;
            (*particle_list)[idx]->clear(); // clean out the previous record

            int pick_flag = 0;
            for(int ipart = 0; ipart < n_particle; ipart++)
            {
                getline(inputfile, temp_string);
                stringstream temp2(temp_string);
                temp2 >> dummy >> dummy >> dummy >> dummy
                      >> dummy >> dummy >> dummy >> dummy
                      >> temp_mass >> urqmd_pid >> urqmd_iso3 >> urqmd_charge;
                pick_flag = decide_to_pick_UrQMD(urqmd_pid, urqmd_iso3, urqmd_charge);
                if(pick_flag == 1)
                {
                     particle_info *temp_particle_info = new particle_info;
                     temp2 >> dummy >> dummy >> dummy;
                     temp2 >> temp_particle_info->t
                           >> temp_particle_info->x 
                           >> temp_particle_info->y
                           >> temp_particle_info->z 
                           >> temp_particle_info->E
                           >> temp_particle_info->px 
                           >> temp_particle_info->py
                           >> temp_particle_info->pz ;
                     temp_particle_info->mass = temp_mass;
                     (*particle_list)[idx]->push_back(*temp_particle_info);
                }
            }
        }
        else
            break;
    }
    return(0);
}

int particleSamples::read_in_particle_samples_Sangwook()
{
    // clean out the previous record
    for(int i = 0; i < particle_list->size(); i++)
        (*particle_list)[i]->clear();
    particle_list->clear();

    string temp_string;
    int n_particle;
    double dummy;
    int ievent;
    int urqmd_pid, urqmd_iso3, urqmd_charge;
    double temp_mass;
    for(ievent = 0; ievent < event_buffer_size; ievent++)
    {
        getline(inputfile, temp_string);
        if(!inputfile.eof())
        {
            particle_list->push_back(new vector<particle_info> );

            // get number of particles within the event
            stringstream temp1(temp_string);
            temp1 >> n_particle;
            getline(inputfile, temp_string);  // then get one useless line

            int idx = ievent;
            (*particle_list)[idx]->clear(); // clean out the previous record

            int pick_flag = 0;
            for(int ipart = 0; ipart < n_particle; ipart++)
            {
                getline(inputfile, temp_string);
                stringstream temp2(temp_string);
                temp2 >> dummy >> dummy >> dummy >> dummy
                      >> dummy >> dummy >> dummy >> dummy
                      >> temp_mass >> urqmd_pid >> urqmd_iso3 >> urqmd_charge;
                pick_flag = decide_to_pick_UrQMD(urqmd_pid, urqmd_iso3, urqmd_charge);
                if(pick_flag == 1)
                {
                     particle_info *temp_particle_info = new particle_info;
                     temp2 >> dummy >> dummy >> dummy;
                     temp2 >> temp_particle_info->t
                           >> temp_particle_info->x 
                           >> temp_particle_info->y
                           >> temp_particle_info->z 
                           >> temp_particle_info->E
                           >> temp_particle_info->px 
                           >> temp_particle_info->py
                           >> temp_particle_info->pz ;
                     temp_particle_info->mass = temp_mass;
                     (*particle_list)[idx]->push_back(*temp_particle_info);
                }
            }
        }
        else
            break;
    }
    return(0);
}
