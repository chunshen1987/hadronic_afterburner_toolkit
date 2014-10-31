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
        int end_event_idx;
        int read_in_mode;
        int *num_of_particles;
        int particle_monval;

        particle_info **particle_list;

    public:
        particleSamples(ParameterReader* paraRdr_in, string path_in);
        ~particleSamples();

        int read_in_particle_samples();
        int read_in_particle_samples_OSCAR();
        bool end_of_file() {return(inputfile.eof());};
        int get_number_of_particles(int event_id) {return(num_of_particles[event_id]);};
        particle_info* get_particle_list(int event_id) {return(particle_list[event_id]);};

};

#endif
