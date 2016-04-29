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
    run_mode = paraRdr->getVal("run_mode");
    if(run_mode == 1)
    {
        reject_decay_flag = paraRdr->getVal("reject_decay_flag");
        tau_reject = paraRdr->getVal("tau_reject");
    }
    else
    {
        reject_decay_flag = 0;
        tau_reject = 10000.;
    }
    
    // read in particle Monte-Carlo number
    particle_monval = paraRdr->getVal("particle_monval");
    flag_isospin = paraRdr->getVal("distinguish_isospin");
    if(read_in_mode == 1 || read_in_mode == 3 || read_in_mode == 4)
        get_UrQMD_id(particle_monval);

    particle_list = new vector< vector<particle_info>* >;
    particle_list_mixed_event = new vector< vector<particle_info>* >;

    ostringstream filename;
    ostringstream filename_mixed_event;
    if(read_in_mode == 0) {
        filename << path << "/OSCAR.DAT";
        filename_mixed_event << path << "/OSCAR_mixed_event.DAT";
    } else if (read_in_mode == 1 || read_in_mode == 3 || read_in_mode == 4
               || read_in_mode == 5) {
        filename << path << "/particle_list.dat";
        filename_mixed_event << path << "/particle_list_mixed_event.dat";
    }

    inputfile.open(filename.str().c_str());
    if(!inputfile.is_open())
    {
        cout << "particleSamples:: Error: input file: " << filename.str() 
             << " can not open!" << endl;
        exit(1);
    }
    if(run_mode == 1)
    {
        inputfile_mixed_event.open(filename_mixed_event.str().c_str());
        if(!inputfile_mixed_event.is_open())
        {
            cout << "particleSamples:: Error: input file: " 
                 << filename_mixed_event.str() << " can not open!" << endl;
            exit(1);
        }
    }

    // skip the header file for OSCAR
    string temp;
    if (read_in_mode == 0) {
        getline(inputfile, temp);
        getline(inputfile, temp);
        getline(inputfile, temp);
        if (run_mode == 1) {
            getline(inputfile_mixed_event, temp);
            getline(inputfile_mixed_event, temp);
            getline(inputfile_mixed_event, temp);
        }
    }
    
    // skip the header in JAM
    if (read_in_mode == 5) {
        getline(inputfile, temp);
    }

    initialize_charged_hadron_urqmd_id_list();
}

particleSamples::~particleSamples()
{
    inputfile.close();
    if(run_mode == 1)
        inputfile_mixed_event.close();
    for(int i = 0; i < particle_list->size(); i++)
        delete (*particle_list)[i];
    delete particle_list;
    for(int i = 0; i < particle_list_mixed_event->size(); i++)
        delete (*particle_list_mixed_event)[i];
    delete particle_list_mixed_event;
}

void particleSamples::initialize_charged_hadron_pdg_list()
{
    charged_hadron_pdg_list[0] = 211;       // pion
    charged_hadron_pdg_list[1] = 321;       // kaon
    charged_hadron_pdg_list[2] = 2212;      // proton
    charged_hadron_pdg_list[3] = 3222;      // Sigma^+
    charged_hadron_pdg_list[4] = 3112;      // Sigma^-
    charged_hadron_pdg_list[5] = 3312;      // Xi^-
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
    else if(monval == -321)  // Kaon^-
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
    else if(read_in_mode == 4)
        read_in_particle_samples_UrQMD_3p3();
    else if(read_in_mode == 5)
        read_in_particle_samples_JAM();

    return(0);
}

int particleSamples::read_in_particle_samples_mixed_event()
{
    if(read_in_mode == 0)
        read_in_particle_samples_OSCAR_mixed_event();
    else if(read_in_mode == 1)
        read_in_particle_samples_UrQMD_mixed_event();
    else if(read_in_mode == 3)
        read_in_particle_samples_mixed_event_Sangwook();
    else if(read_in_mode == 4)
        read_in_particle_samples_UrQMD_3p3_mixed_event();
    else if(read_in_mode == 5)
        read_in_particle_samples_JAM_mixed_event();

    return(0);
}

int particleSamples::decide_to_pick_UrQMD(int pid, int iso3, int charge,
                                          int parent_proc_type) {
    int pick_flag = 0;
    if (particle_urqmd_id == 9999) {  // charged hadrons
        int in_flag = 0;
        for (int i = 0; i < 5; i++) {
            if (abs(pid) == charged_hadron_urqmd_id_list[i]) {
                in_flag = 1;
                break;
            }
        }
        if (in_flag == 1 && charge != 0)
            pick_flag = 1;
    } else {
        if (flag_isospin == 0) {
            if (pid == particle_urqmd_id) {
                if (reject_decay_flag == 1 && parent_proc_type == 20) {
                    pick_flag = 0;
                } else {
                    pick_flag = 1;
                }
            }
        } else {
            if (pid == particle_urqmd_id && iso3 == particle_urqmd_isospin) {
                if (reject_decay_flag == 1 && parent_proc_type == 20) {
                    pick_flag = 0;
                } else {
                    pick_flag = 1;
                }
            }
        }
    }
    return(pick_flag);
}

int particleSamples::decide_to_pick_JAM(int pid) {
    int pick_flag = 0;
    if (particle_urqmd_id == 9999) {  // charged hadrons
        int in_flag = 0;
        for (int i = 0; i < 6; i++) {
            if (abs(pid) == charged_hadron_pdg_list[i]) {
                in_flag = 1;
                break;
            }
        }
    }
    return(pick_flag);
}

int particleSamples::read_in_particle_samples_OSCAR() {
    // clean out the previous record
    for (int i = 0; i < particle_list->size(); i++)
        (*particle_list)[i]->clear();
    particle_list->clear();
    
    string temp_string;
    int event_id, n_particle, dummy;
    int ievent;
    int temp_monval;
    for (ievent = 0; ievent < event_buffer_size; ievent++) {
        getline(inputfile, temp_string);
        stringstream temp1(temp_string);
        temp1 >> event_id >> n_particle;
        if (!inputfile.eof()) {
            particle_list->push_back(new vector<particle_info> );
            int idx = ievent;

            int pick_flag = 0;
            for (int ipart = 0; ipart < n_particle; ipart++) {
                getline(inputfile, temp_string);
                stringstream temp2(temp_string);
                temp2 >> dummy >> temp_monval;
                if (flag_isospin == 0) {
                    if (abs(temp_monval) == particle_monval)
                        pick_flag = 1;
                    else
                        pick_flag = 0;
                } else {
                    if (temp_monval == particle_monval)
                        pick_flag = 1;
                    else
                        pick_flag = 0;
                }
                if (pick_flag == 1) {
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
        } else {
            break;
        }
    }
    return(0);
}

int particleSamples::read_in_particle_samples_JAM() {
    // clean out the previous record
    for (int i = 0; i < particle_list->size(); i++)
        (*particle_list)[i]->clear();
    particle_list->clear();
    
    string temp_string;
    int event_id, n_particle, dummy;
    char cdummy;
    int ievent;
    int temp_monval;
    for (ievent = 0; ievent < event_buffer_size; ievent++) {
        getline(inputfile, temp_string);
        stringstream temp1(temp_string);
        temp1 >> cdummy >> event_id >> n_particle;
        if (!inputfile.eof()) {
            particle_list->push_back(new vector<particle_info> );
            int idx = ievent;

            int pick_flag = 0;
            for (int ipart = 0; ipart < n_particle; ipart++) {
                getline(inputfile, temp_string);
                stringstream temp2(temp_string);
                temp2 >> temp_monval;
                if (flag_isospin == 0) {
                    if (abs(temp_monval) == particle_monval)
                        pick_flag = 1;
                    else
                        pick_flag = 0;
                } else {
                    if (temp_monval == particle_monval)
                        pick_flag = 1;
                    else
                        pick_flag = 0;
                }
                if (particle_monval == 9999)
                    pick_flag = decide_to_pick_JAM(temp_monval);
                if (pick_flag == 1) {  // pick up
                     particle_info *temp_particle_info = new particle_info;
                     temp2 >> temp_particle_info->mass
                           >> temp_particle_info->px
                           >> temp_particle_info->py
                           >> temp_particle_info->pz
                           >> temp_particle_info->x 
                           >> temp_particle_info->y
                           >> temp_particle_info->z 
                           >> temp_particle_info->t;
                     temp_particle_info->E = sqrt(
                         temp_particle_info->mass*temp_particle_info->mass
                         + temp_particle_info->px*temp_particle_info->px
                         + temp_particle_info->py*temp_particle_info->py
                         + temp_particle_info->pz*temp_particle_info->pz);
                     (*particle_list)[idx]->push_back(*temp_particle_info);
                }
            }
        } else {
            break;
        }
    }
    return(0);
}

int particleSamples::read_in_particle_samples_OSCAR_mixed_event() {
    // clean out the previous record
    for (int i = 0; i < particle_list_mixed_event->size(); i++)
        (*particle_list_mixed_event)[i]->clear();
    particle_list_mixed_event->clear();
    
    string temp_string;
    int event_id, n_particle, dummy;
    int ievent;
    int temp_monval;
    for (ievent = 0; ievent < event_buffer_size; ievent++) {
        getline(inputfile_mixed_event, temp_string);
        stringstream temp1(temp_string);
        temp1 >> event_id >> n_particle;
        if (!inputfile_mixed_event.eof()) {
            particle_list_mixed_event->push_back(new vector<particle_info> );
            int idx = ievent;

            int pick_flag = 0;
            for (int ipart = 0; ipart < n_particle; ipart++) {
                getline(inputfile_mixed_event, temp_string);
                stringstream temp2(temp_string);
                temp2 >> dummy >> temp_monval;
                if (flag_isospin == 0) {
                    if (abs(temp_monval) == particle_monval)
                        pick_flag = 1;
                    else
                        pick_flag = 0;
                } else {
                    if (temp_monval == particle_monval)
                        pick_flag = 1;
                    else
                        pick_flag = 0;
                }
                if (pick_flag == 1) {
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
                     (*particle_list_mixed_event)[idx]->push_back(
                                                        *temp_particle_info);
                }
            }
        } else {
            break;
        }
    }
    return(0);
}

int particleSamples::read_in_particle_samples_JAM_mixed_event() {
    // clean out the previous record
    for (int i = 0; i < particle_list_mixed_event->size(); i++)
        (*particle_list_mixed_event)[i]->clear();
    particle_list_mixed_event->clear();
    
    string temp_string;
    int event_id, n_particle, dummy;
    char cdummy;
    int ievent;
    int temp_monval;
    for (ievent = 0; ievent < event_buffer_size; ievent++) {
        getline(inputfile_mixed_event, temp_string);
        stringstream temp1(temp_string);
        temp1 >> cdummy >> event_id >> n_particle;
        if (!inputfile_mixed_event.eof()) {
            particle_list_mixed_event->push_back(new vector<particle_info> );
            int idx = ievent;

            int pick_flag = 0;
            for (int ipart = 0; ipart < n_particle; ipart++) {
                getline(inputfile_mixed_event, temp_string);
                stringstream temp2(temp_string);
                temp2 >> temp_monval;
                if (flag_isospin == 0) {
                    if (abs(temp_monval) == particle_monval)
                        pick_flag = 1;
                    else
                        pick_flag = 0;
                } else {
                    if (temp_monval == particle_monval)
                        pick_flag = 1;
                    else
                        pick_flag = 0;
                }
                if (particle_monval == 9999)
                    pick_flag = decide_to_pick_JAM(temp_monval);
                if (pick_flag == 1) {
                     particle_info *temp_particle_info = new particle_info;
                     temp2 >> temp_particle_info->mass
                           >> temp_particle_info->px
                           >> temp_particle_info->py
                           >> temp_particle_info->pz
                           >> temp_particle_info->x
                           >> temp_particle_info->y
                           >> temp_particle_info->z
                           >> temp_particle_info->t;
                     temp_particle_info->E = sqrt(
                         temp_particle_info->mass*temp_particle_info->mass
                         + temp_particle_info->px*temp_particle_info->px
                         + temp_particle_info->py*temp_particle_info->py
                         + temp_particle_info->pz*temp_particle_info->pz);
                     (*particle_list_mixed_event)[idx]->push_back(
                                                        *temp_particle_info);
                }
            }
        } else {
            break;
        }
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
    int parent_proc_type;
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
            for(int i = 0; i < 16; i++)
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
                      >> temp_mass >> urqmd_pid >> urqmd_iso3 >> urqmd_charge
                      >> dummy >> dummy >> parent_proc_type;
                pick_flag = decide_to_pick_UrQMD(
                        urqmd_pid, urqmd_iso3, urqmd_charge, parent_proc_type);
                if(pick_flag == 1)
                {
                     particle_info *temp_particle_info = new particle_info;
                     temp2 >> temp_particle_info->t
                           >> temp_particle_info->x 
                           >> temp_particle_info->y
                           >> temp_particle_info->z 
                           >> temp_particle_info->E
                           >> temp_particle_info->px 
                           >> temp_particle_info->py
                           >> temp_particle_info->pz;
                     temp_particle_info->mass = temp_mass;
                     if(reject_decay_flag == 2 && parent_proc_type == 20)
                     {
                         double tau = sqrt(
                             temp_particle_info->t*temp_particle_info->t
                             - temp_particle_info->z*temp_particle_info->z);
                         if(tau > tau_reject)
                         {
                             pick_flag = 0;
                         }
                         else
                         {
                             pick_flag = 1;
                         }
                     }
                     if(pick_flag == 1)
                     {
                         (*particle_list)[idx]->push_back(*temp_particle_info);
                     }
                     else
                     {
                         delete temp_particle_info;
                     }
                }
            }
        }
        else
            break;
    }
    return(0);
}

int particleSamples::read_in_particle_samples_UrQMD_3p3()
{
    // clean out the previous record
    for(int i = 0; i < particle_list->size(); i++)
        (*particle_list)[i]->clear();
    particle_list->clear();

    string temp_string;
    int n_particle;
    double dummy;
    int parent_proc_type;
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
                      >> temp_mass >> urqmd_pid >> urqmd_iso3 >> urqmd_charge
                      >> dummy >> dummy >> parent_proc_type;
                pick_flag = decide_to_pick_UrQMD(
                        urqmd_pid, urqmd_iso3, urqmd_charge, parent_proc_type);
                if(pick_flag == 1)
                {
                     particle_info *temp_particle_info = new particle_info;
                     temp2 >> temp_particle_info->t
                           >> temp_particle_info->x 
                           >> temp_particle_info->y
                           >> temp_particle_info->z 
                           >> temp_particle_info->E
                           >> temp_particle_info->px 
                           >> temp_particle_info->py
                           >> temp_particle_info->pz;
                     temp_particle_info->mass = temp_mass;
                     if(reject_decay_flag == 2 && parent_proc_type == 20)
                     {
                         double tau = sqrt(
                             temp_particle_info->t*temp_particle_info->t
                             - temp_particle_info->z*temp_particle_info->z);
                         if(tau > tau_reject)
                         {
                             pick_flag = 0;
                         }
                         else
                         {
                             pick_flag = 1;
                         }
                     }
                     if(pick_flag == 1)
                     {
                         (*particle_list)[idx]->push_back(*temp_particle_info);
                     }
                     else
                     {
                         delete temp_particle_info;
                     }
                }
            }
        }
        else
            break;
    }
    return(0);
}

int particleSamples::read_in_particle_samples_UrQMD_mixed_event()
{
    // clean out the previous record
    for(int i = 0; i < particle_list_mixed_event->size(); i++)
        (*particle_list_mixed_event)[i]->clear();
    particle_list_mixed_event->clear();

    string temp_string;
    int n_particle;
    double dummy;
    int parent_proc_type;
    int ievent;
    int urqmd_pid, urqmd_iso3, urqmd_charge;
    double temp_mass;
    for(ievent = 0; ievent < event_buffer_size; ievent++)
    {
        getline(inputfile_mixed_event, temp_string);
        if(!inputfile_mixed_event.eof())
        {
            particle_list_mixed_event->push_back(new vector<particle_info> );

            // first skip the header
            for(int i = 0; i < 16; i++)
                getline(inputfile_mixed_event, temp_string);
            // then get number of particles within the event
            getline(inputfile_mixed_event, temp_string);
            stringstream temp1(temp_string);
            temp1 >> n_particle;
            // then get one useless line
            getline(inputfile_mixed_event, temp_string);  

            // clean out the previous record
            int idx = ievent;
            (*particle_list_mixed_event)[idx]->clear(); 

            int pick_flag = 0;
            for(int ipart = 0; ipart < n_particle; ipart++)
            {
                getline(inputfile_mixed_event, temp_string);
                stringstream temp2(temp_string);
                temp2 >> dummy >> dummy >> dummy >> dummy
                      >> dummy >> dummy >> dummy >> dummy
                      >> temp_mass >> urqmd_pid >> urqmd_iso3 >> urqmd_charge
                      >> dummy >> dummy >> parent_proc_type;
                pick_flag = decide_to_pick_UrQMD(
                        urqmd_pid, urqmd_iso3, urqmd_charge, parent_proc_type);
                if(pick_flag == 1)
                {
                     particle_info *temp_particle_info = new particle_info;
                     temp2 >> temp_particle_info->t
                           >> temp_particle_info->x 
                           >> temp_particle_info->y
                           >> temp_particle_info->z 
                           >> temp_particle_info->E
                           >> temp_particle_info->px 
                           >> temp_particle_info->py
                           >> temp_particle_info->pz ;
                     temp_particle_info->mass = temp_mass;
                     if(reject_decay_flag == 2 && parent_proc_type == 20)
                     {
                         double tau = sqrt(
                             temp_particle_info->t*temp_particle_info->t
                             - temp_particle_info->z*temp_particle_info->z);
                         if(tau > tau_reject)
                         {
                             pick_flag = 0;
                         }
                         else
                         {
                             pick_flag = 1;
                         }
                     }
                     if(pick_flag == 1)
                     {
                         (*particle_list_mixed_event)[idx]->push_back(
                                                        *temp_particle_info);
                     }
                     else
                     {
                         delete temp_particle_info;
                     }
                }
            }
        }
        else
            break;
    }
    return(0);
}

int particleSamples::read_in_particle_samples_UrQMD_3p3_mixed_event()
{
    // clean out the previous record
    for(int i = 0; i < particle_list_mixed_event->size(); i++)
        (*particle_list_mixed_event)[i]->clear();
    particle_list_mixed_event->clear();

    string temp_string;
    int n_particle;
    double dummy;
    int parent_proc_type;
    int ievent;
    int urqmd_pid, urqmd_iso3, urqmd_charge;
    double temp_mass;
    for(ievent = 0; ievent < event_buffer_size; ievent++)
    {
        getline(inputfile_mixed_event, temp_string);
        if(!inputfile_mixed_event.eof())
        {
            particle_list_mixed_event->push_back(new vector<particle_info> );

            // first skip the header
            for(int i = 0; i < 13; i++)
                getline(inputfile_mixed_event, temp_string);
            // then get number of particles within the event
            getline(inputfile_mixed_event, temp_string);
            stringstream temp1(temp_string);
            temp1 >> n_particle;
            // then get one useless line
            getline(inputfile_mixed_event, temp_string);  

            // clean out the previous record
            int idx = ievent;
            (*particle_list_mixed_event)[idx]->clear(); 

            int pick_flag = 0;
            for(int ipart = 0; ipart < n_particle; ipart++)
            {
                getline(inputfile_mixed_event, temp_string);
                stringstream temp2(temp_string);
                temp2 >> dummy >> dummy >> dummy >> dummy
                      >> dummy >> dummy >> dummy >> dummy
                      >> temp_mass >> urqmd_pid >> urqmd_iso3 >> urqmd_charge
                      >> dummy >> dummy >> parent_proc_type;
                pick_flag = decide_to_pick_UrQMD(
                        urqmd_pid, urqmd_iso3, urqmd_charge, parent_proc_type);
                if(pick_flag == 1)
                {
                     particle_info *temp_particle_info = new particle_info;
                     temp2 >> temp_particle_info->t
                           >> temp_particle_info->x 
                           >> temp_particle_info->y
                           >> temp_particle_info->z 
                           >> temp_particle_info->E
                           >> temp_particle_info->px 
                           >> temp_particle_info->py
                           >> temp_particle_info->pz ;
                     temp_particle_info->mass = temp_mass;
                     if(reject_decay_flag == 2 && parent_proc_type == 20)
                     {
                         double tau = sqrt(
                             temp_particle_info->t*temp_particle_info->t
                             - temp_particle_info->z*temp_particle_info->z);
                         if(tau > tau_reject)
                         {
                             pick_flag = 0;
                         }
                         else
                         {
                             pick_flag = 1;
                         }
                     }
                     if(pick_flag == 1)
                     {
                         (*particle_list_mixed_event)[idx]->push_back(
                                                        *temp_particle_info);
                     }
                     else
                     {
                         delete temp_particle_info;
                     }
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
    int parent_proc_type;
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
                      >> temp_mass >> urqmd_pid >> urqmd_iso3 >> urqmd_charge
                      >> dummy >> dummy >> parent_proc_type;

                pick_flag = decide_to_pick_UrQMD(
                        urqmd_pid, urqmd_iso3, urqmd_charge, parent_proc_type);
                if(pick_flag == 1)
                {
                    particle_info *temp_particle_info = new particle_info;
                    temp2 >> temp_particle_info->t
                          >> temp_particle_info->x 
                          >> temp_particle_info->y
                          >> temp_particle_info->z 
                          >> temp_particle_info->E
                          >> temp_particle_info->px 
                          >> temp_particle_info->py
                          >> temp_particle_info->pz ;
                    temp_particle_info->mass = temp_mass;
                    if(reject_decay_flag == 2 && parent_proc_type == 20)
                    {
                        double tau = sqrt(
                            temp_particle_info->t*temp_particle_info->t
                            - temp_particle_info->z*temp_particle_info->z);
                        if(tau > tau_reject)
                        {
                            pick_flag = 0;
                        }
                        else
                        {
                            pick_flag = 1;
                        }
                    }
                    if(pick_flag == 1)
                    {
                        (*particle_list)[idx]->push_back(*temp_particle_info);
                    }
                    else
                    {
                        delete temp_particle_info;
                    }
                }
            }
        }
        else
            break;
    }
    return(0);
}

int particleSamples::read_in_particle_samples_mixed_event_Sangwook()
{
    // clean out the previous record
    for(int i = 0; i < particle_list_mixed_event->size(); i++)
        (*particle_list_mixed_event)[i]->clear();
    particle_list_mixed_event->clear();

    string temp_string;
    int n_particle;
    double dummy;
    int parent_proc_type;
    int ievent;
    int urqmd_pid, urqmd_iso3, urqmd_charge;
    double temp_mass;
    for(ievent = 0; ievent < event_buffer_size; ievent++)
    {
        getline(inputfile_mixed_event, temp_string);
        if(!inputfile_mixed_event.eof())
        {
            particle_list_mixed_event->push_back(new vector<particle_info> );

            // get number of particles within the event
            stringstream temp1(temp_string);
            temp1 >> n_particle;
            // then get one useless line
            getline(inputfile_mixed_event, temp_string);

            // clean out the previous record
            int idx = ievent;
            (*particle_list_mixed_event)[idx]->clear(); 

            int pick_flag = 0;
            for(int ipart = 0; ipart < n_particle; ipart++)
            {
                getline(inputfile_mixed_event, temp_string);
                stringstream temp2(temp_string);
                temp2 >> dummy >> dummy >> dummy >> dummy
                      >> dummy >> dummy >> dummy >> dummy
                      >> temp_mass >> urqmd_pid >> urqmd_iso3 >> urqmd_charge
                      >> dummy >> dummy >> parent_proc_type;
                pick_flag = decide_to_pick_UrQMD(
                        urqmd_pid, urqmd_iso3, urqmd_charge, parent_proc_type);
                if(pick_flag == 1)
                {
                    particle_info *temp_particle_info = new particle_info;
                    temp2 >> temp_particle_info->t
                          >> temp_particle_info->x 
                          >> temp_particle_info->y
                          >> temp_particle_info->z 
                          >> temp_particle_info->E
                          >> temp_particle_info->px 
                          >> temp_particle_info->py
                          >> temp_particle_info->pz ;
                    temp_particle_info->mass = temp_mass;
                    if(reject_decay_flag == 2 && parent_proc_type == 20)
                    {
                        double tau = sqrt(
                            temp_particle_info->t*temp_particle_info->t
                            - temp_particle_info->z*temp_particle_info->z);
                        if(tau > tau_reject)
                        {
                            pick_flag = 0;
                        }
                        else
                        {
                            pick_flag = 1;
                        }
                    }
                    if(pick_flag == 1)
                    {
                        (*particle_list_mixed_event)[idx]->push_back(
                                                        *temp_particle_info);
                    }
                    else
                    {
                        delete temp_particle_info;
                    }
                }
            }
        }
        else
            break;
    }
    return(0);
}
