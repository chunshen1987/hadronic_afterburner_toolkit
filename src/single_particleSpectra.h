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
        double *pT_array, *pT_mean_array, *pT_mean_array_err;

        int rap_type;
        double rap_min, rap_max;

        int total_number_of_events;
        double *Qn_vector_real, *Qn_vector_imag;
        double **Qn_diff_vector_real, **Qn_diff_vector_imag;
        double *Qn_vector_real_err, *Qn_vector_imag_err;
        double **Qn_diff_vector_real_err, **Qn_diff_vector_imag_err;

        int check_spatial_flag;
        int N_tau;
        double tau_min, tau_max, dtau;
        double *tau_array, *dNdtau_array;
        int N_xpt;
        double spatial_x_min, spatial_x_max, dspatial_x;
        double *xpt_array, *dNdx1_array;
        double *ypt_array, *dNdx2_array;
        int N_eta_s;
        double eta_s_min, eta_s_max, deta_s;
        double *eta_s_array, *dNdetas_array;

    public:
        singleParticleSpectra(ParameterReader *paraRdr_in, string path_in, particleSamples *particle_list_in);
        ~singleParticleSpectra();

        void calculate_Qn_vector_shell();
        void calculate_Qn_vector(int event_id);
        void output_Qn_vectors();
        
        void check_dNdSV(int event_id);
        void output_dNdSV();

};

#endif
