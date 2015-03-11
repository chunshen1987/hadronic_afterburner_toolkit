#ifndef particleSamples_h
#define particleSamples_h

#include <fstream>
#include <vector>

#include "ParameterReader.h"
using namespace std;

struct particle_info
{
    double mass;
    double px, py, pz, E;
    double x, y, z, t;
};

class particleSamples
{
    private:
        ParameterReader *paraRdr;
        string path;                 // path for results folder
        ifstream inputfile;
        int event_buffer_size;
        int read_in_mode;
        int particle_monval;
        int particle_urqmd_id, particle_urqmd_isospin;

        vector< vector<particle_info>* >* particle_list;

    public:
        particleSamples(ParameterReader* paraRdr_in, string path_in);
        ~particleSamples();

        void get_UrQMD_id(int monval);
        int read_in_particle_samples();
        int read_in_particle_samples_OSCAR();
        int read_in_particle_samples_UrQMD();
        int read_in_particle_samples_Sangwook();
        bool end_of_file() {return(inputfile.eof());};
        int get_event_buffer_size() {return(event_buffer_size);};
        int get_number_of_events() {return(particle_list->size());};
        int get_number_of_particles(int event_id) {return((*particle_list)[event_id]->size());};
        particle_info get_particle(int event_id, int part_id) {return((*(*particle_list)[event_id])[part_id]);};

};

#endif
