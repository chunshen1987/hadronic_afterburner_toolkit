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
    if(read_in_mode == 1)
        get_UrQMD_id(particle_monval);

    particle_list = new vector< vector<particle_info>* >;
//    for(int i = 0; i < event_buffer_size; i++)
//        particle_list->push_back(new vector<particle_info> );

    ostringstream filename;
    if(read_in_mode == 0)
        filename << path << "/OSCAR.DAT";
    else if(read_in_mode == 1)
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
}

particleSamples::~particleSamples()
{
    inputfile.close();
    for(int i = 0; i < particle_list->size(); i++)
        delete (*particle_list)[i];
    delete particle_list;
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
    else if(monval == 2212)  // proton
    {
        particle_urqmd_id = 1;
        particle_urqmd_isospin = 1;
    }
}

int particleSamples::read_in_particle_samples()
{
    if(read_in_mode == 0)
        read_in_particle_samples_OSCAR();
    else if(read_in_mode == 1)
        read_in_particle_samples_UrQMD();

    return(0);
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

            for(int ipart = 0; ipart < n_particle; ipart++)
            {
                getline(inputfile, temp_string);
                stringstream temp2(temp_string);
                temp2 >> dummy >> temp_monval;
                if(temp_monval == particle_monval)
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
    int urqmd_pid, urqmd_iso3;
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

            for(int ipart = 0; ipart < n_particle; ipart++)
            {
                getline(inputfile, temp_string);
                stringstream temp2(temp_string);
                temp2 >> dummy >> dummy >> dummy >> dummy
                      >> dummy >> dummy >> dummy >> dummy
                      >> temp_mass >> urqmd_pid >> urqmd_iso3;
                if(urqmd_pid == particle_urqmd_id && urqmd_iso3 == particle_urqmd_isospin)
                {
                     particle_info *temp_particle_info = new particle_info;
                     temp2 >> dummy >> dummy >> dummy >> dummy;
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
