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

        int particle_monval;

        int order_max;

        int npT;
        double pT_min, pT_max, dpT;
        double *pT_array, *pT_mean_array, *pT_mean_array_err;

        int rap_type;
        double rap_min, rap_max;
        int rapidity_distribution_flag;
        int N_rap;
        double rapidity_dis_min, rapidity_dis_max, drap;
        double *rapidity_array, *dNdy_array;
        double vn_rapidity_dis_pT_min, vn_rapidity_dis_pT_max;
        double **vn_real_rapidity_dis_array, **vn_imag_rapidity_dis_array;
        double **vn_real_rapidity_dis_array_err;
        double **vn_imag_rapidity_dis_array_err;

        int total_number_of_events;
        double *Qn_vector_real, *Qn_vector_imag;
        double **Qn_diff_vector_real, **Qn_diff_vector_imag;
        double *Qn_vector_real_err, *Qn_vector_imag_err;
        double **Qn_diff_vector_real_err, **Qn_diff_vector_imag_err;

        int check_spatial_flag;
        int N_tau;
        double tau_min, tau_max, dtau, intrinsic_dtau;
        double *tau_array, *dNdtau_array;
        int N_xpt;
        double spatial_x_min, spatial_x_max, dspatial_x, intrinsic_dx;
        double *xpt_array, *dNdx1_array;
        double *ypt_array, *dNdx2_array;
        double **dNdtaudx1_array, **dNdtaudx2_array;
        int N_eta_s;
        double eta_s_min, eta_s_max, deta_s, intrinsic_detas;
        double *eta_s_array, *dNdetas_array;

    public:
        singleParticleSpectra(ParameterReader *paraRdr_in, string path_in, 
                              particleSamples *particle_list_in);
        ~singleParticleSpectra();

        void calculate_Qn_vector_shell();
        void calculate_Qn_vector(int event_id);
        void output_Qn_vectors();

        void calculate_rapidity_distribution(int event_id);
        void output_rapidity_distribution();
        
        void check_dNdSV(int event_id);
        void output_dNdSV();

};

#endif
