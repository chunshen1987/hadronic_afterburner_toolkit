// Copyright Chun Shen @ 2016
#ifndef SRC_particleSamples_h_
#define SRC_particleSamples_h_

#include <fstream>
#include <vector>

#include "./ParameterReader.h"
using namespace std;

struct particle_info {
    double mass;
    double px, py, pz, E;
    double x, y, z, t;
};

class particleSamples {
 private:
    ParameterReader *paraRdr;
    string path;                 // path for results folder
    ifstream inputfile;
    ifstream inputfile_mixed_event;
    int event_buffer_size;
    int read_in_mode;
    int run_mode;
    int particle_monval;
    int resonance_feed_down_flag;  // include Sigma0 feed down of Lambda
    int reconst_flag;  // reconst phi meson from (K^+, K^-) pairs
    int net_particle_flag;  // flag to collect net particle distribution
    int flag_isospin;
    int reject_decay_flag;
    double tau_reject;
    int particle_urqmd_id, particle_urqmd_isospin;

    int charged_hadron_pdg_list[6];
    int charged_hadron_urqmd_id_list[5];

    vector< vector<particle_info>* >* particle_list;
    vector< vector<particle_info>* >* anti_particle_list;
    vector< vector<particle_info>* >* particle_list_mixed_event;
    vector< vector<particle_info>* >* resonance_list;
    vector< vector<particle_info>* >* reconst_list_1;
    vector< vector<particle_info>* >* reconst_list_2;

 public:
    particleSamples(ParameterReader* paraRdr_in, string path_in);
    ~particleSamples();

    void initialize_charged_hadron_urqmd_id_list();
    void initialize_charged_hadron_pdg_list();
    void get_UrQMD_id(int monval);
    int decide_to_pick_UrQMD(int pid, int iso3, int charge,
                             int parent_proc_type);
    int decide_to_pick_UrQMD_resonance(int pid, int iso3, int charge);
    void decide_to_pick_UrQMD_reconst(
                int pid, int iso3, int charge, int parent_proc_type,
                int *flag1, int *flag2);
    int decide_to_pick_JAM(int pid);
    int decide_to_pick_UrQMD_anti_particles(int pid, int iso3);

    void perform_resonance_feed_down();
    void perform_two_body_decay(particle_info *mother,
                                particle_info* daughter1,
                                particle_info* daughter2);
    void perform_particle_reconstruction(); 
    int read_in_particle_samples();
    int read_in_particle_samples_mixed_event();
    int read_in_particle_samples_OSCAR();
    int read_in_particle_samples_OSCAR_mixed_event();
    int read_in_particle_samples_JAM();
    int read_in_particle_samples_JAM_mixed_event();
    int read_in_particle_samples_UrQMD();
    int read_in_particle_samples_UrQMD_mixed_event();
    int read_in_particle_samples_UrQMD_3p3();
    int read_in_particle_samples_UrQMD_3p3_mixed_event();
    int read_in_particle_samples_Sangwook();
    int read_in_particle_samples_mixed_event_Sangwook();

    bool end_of_file() {return(inputfile.eof());}
    bool end_of_file_mixed_event() {return(inputfile_mixed_event.eof());}

    int get_event_buffer_size() {return(event_buffer_size);}

    int get_number_of_events() {return(particle_list->size());}
    int get_number_of_events_anti_particle()
    {return(anti_particle_list->size());}
    int get_number_of_mixed_events()
    {return(particle_list_mixed_event->size());}

    int get_number_of_particles(int event_id)
    {return((*particle_list)[event_id]->size());}
    int get_number_of_particles_mixed_event(int event_id)
    {return((*particle_list_mixed_event)[event_id]->size());}
    int get_number_of_anti_particles(int event_id)
    {return((*anti_particle_list)[event_id]->size());}

    particle_info get_particle(int event_id, int part_id) 
    {return((*(*particle_list)[event_id])[part_id]);}
    particle_info get_anti_particle(int event_id, int part_id) 
    {return((*(*anti_particle_list)[event_id])[part_id]);}
    particle_info get_particle_from_mixed_event(int event_id, int part_id) 
    {return((*(*particle_list_mixed_event)[event_id])[part_id]);}
};

#endif  // SRC_particleSamples_h_
