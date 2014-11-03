#ifndef single_particleSpectra_h
#define single_particleSpectra_h

#include "ParameterReader.h"
#include "particleSamples.h"

using namespace std;

class singleParticleSpectra
{
    private:
        ParameterReader *paraRdr;
        string path;
        particleSamples *particle_list;

        int order_max;

        int npT;
        double pT_min, pT_max, dpT;
        double *pT_array;

        int rap_type;
        double rap_min, rap_max;

        int total_number_of_events;
        double *Qn_vector_real, *Qn_vector_imag;
        double **Qn_diff_vector_real, **Qn_diff_vector_imag;


    public:
        singleParticleSpectra(ParameterReader *paraRdr_in, string path_in, particleSamples *particle_list_in);
        ~singleParticleSpectra();

        void calculate_Qn_vector_shell();
        void calculate_Qn_vector(int event_id);
        void output_Qn_vectors();

};

#endif
