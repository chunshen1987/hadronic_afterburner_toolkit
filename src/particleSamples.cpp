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

    particle_monval = paraRdr->getVal("particle_monval");
    if(particle_monval == 211)
    {
        particle_urqmd_id = 101;
        particle_urqmd_isospin = 2;
    }
    else if(particle_monval == -211)
    {
        particle_urqmd_id = 101;
        particle_urqmd_isospin = -2;
    }
    else if(particle_monval == 321)
    {
        particle_urqmd_id = 106;
        particle_urqmd_isospin = 1;
    }
    else if(particle_monval == 2212)
    {
        particle_urqmd_id = 1;
        particle_urqmd_isospin = 1;
    }
    event_buffer_size = paraRdr->getVal("event_buffer_size");
    read_in_mode = paraRdr->getVal("read_in_mode");

    num_of_particles = new int [event_buffer_size];
    for(int i = 0; i < event_buffer_size; i++)
        num_of_particles[i] = 0;

    particle_list = new particle_info* [event_buffer_size];
    for(int i = 0; i < event_buffer_size; i++)
        particle_list[i] = new particle_info [2];

    ostringstream filename;
    if(read_in_mode == 0)
        filename << path << "/OSCAR.DAT";
    else
        filename << path << "/OSCAR.DAT";   // to be changed

    inputfile.open(filename.str().c_str());
    // skip the header file
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
    delete [] num_of_particles;

    for(int i = 0; i < event_buffer_size; i++)
        delete [] particle_list[i];
    delete [] particle_list;
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
            int idx = ievent;
            delete [] particle_list[idx]; // clean out the previous record
            particle_list[idx] = new particle_info [n_particle];

            int num_of_chosen_particle = 0;
            for(int ipart = 0; ipart < n_particle; ipart++)
            {
                getline(inputfile, temp_string);
                stringstream temp2(temp_string);
                temp2 >> dummy >> temp_monval;
                if(temp_monval == particle_monval)
                {
                     num_of_chosen_particle++;
                     temp2 >> particle_list[idx][ipart].px 
                           >> particle_list[idx][ipart].py
                           >> particle_list[idx][ipart].pz 
                           >> particle_list[idx][ipart].E
                           >> particle_list[idx][ipart].mass 
                           >> particle_list[idx][ipart].x 
                           >> particle_list[idx][ipart].y
                           >> particle_list[idx][ipart].z 
                           >> particle_list[idx][ipart].t;
                }
            }
            num_of_particles[idx] = num_of_chosen_particle;
        }
        else
            break;
    }
    end_event_idx = ievent - 1;
    return(0);
}

int particleSamples::read_in_particle_samples_UrQMD()
{
    string temp_string;
    
    int n_particle, dummy;
    int ievent;
    int urqmd_pid, urqmd_iso3;
    double temp_mass;
    for(ievent = 0; ievent < event_buffer_size; ievent++)
    {
        getline(inputfile, temp_string);
        if(!inputfile.eof())
        {
            // first skip the header
            for(int i = 0; i < 13; i++)
                getline(inputfile, temp_string);
            // then get number of particles within the event
            getline(inputfile, temp_string);
            stringstream temp1(temp_string);
            temp1 >> n_particle;

            int idx = ievent;
            delete [] particle_list[idx]; // clean out the previous record
            particle_list[idx] = new particle_info [n_particle];

            int num_of_chosen_particle = 0;
            for(int ipart = 0; ipart < n_particle; ipart++)
            {
                getline(inputfile, temp_string);
                stringstream temp2(temp_string);
                temp2 >> dummy >> dummy >> dummy >> dummy
                      >> dummy >> dummy >> dummy >> dummy
                      >> temp_mass >> urqmd_pid >> urqmd_iso3;
                if(urqmd_pid == particle_urqmd_id && urqmd_iso3 == particle_urqmd_isospin)
                {
                     num_of_chosen_particle++;
                     temp2 >> dummy >> dummy >> dummy >> dummy;
                     temp2 >> particle_list[idx][ipart].t
                           >> particle_list[idx][ipart].x 
                           >> particle_list[idx][ipart].y
                           >> particle_list[idx][ipart].z 
                           >> particle_list[idx][ipart].E
                           >> particle_list[idx][ipart].px 
                           >> particle_list[idx][ipart].py
                           >> particle_list[idx][ipart].pz ;
                     particle_list[idx][ipart].mass = temp_mass;
                }
            }
            num_of_particles[idx] = num_of_chosen_particle;
        }
        else
            break;
    }
    end_event_idx = ievent - 1;
    return(0);
}
