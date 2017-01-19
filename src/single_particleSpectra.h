#ifndef SRC_single_particleSpectra_h_
#define SRC_single_particleSpectra_h_

#include "./ParameterReader.h"
#include "./particleSamples.h"

using namespace std;

class singleParticleSpectra {
 private:
    ParameterReader *paraRdr;
    string path;
    particleSamples *particle_list;

    int particle_monval;

    double reconst_branching_ratio;

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

    int flag_correlation;
    double *Qn2_vector, *Qn2_vector_err;
    double **QnSP_diff_vector, **QnSP_diff_vector_err;
    int num_corr;
    double *C_nmk, *C_nmk_err;
    int flag_charge_dependence;
    double *C_nmk_ss, *C_nmk_ss_err, *C_nmk_os, *C_nmk_os_err;

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

    //! This is a driver function to compute the Qn flow vector
    void calculate_Qn_vector_shell();

    //! this function computes the pT-integrated and pT-differential Qn vector
    //! within a given rapidity region in one event
    void calculate_Qn_vector(int event_id,
            double *event_Qn_real, double *event_Qn_real_err,
            double *event_Qn_imag, double *event_Qn_imag_err,
            double **event_Qn_diff_real, double **event_Qn_diff_real_err,
            double **event_Qn_diff_imag, double **event_Qn_diff_imag_err);
    void calculate_Qn_vector_positive_charge(int event_id,
        double *event_Qn_real, double *event_Qn_real_err,
        double *event_Qn_imag, double *event_Qn_imag_err,
        double **event_Qn_diff_real, double **event_Qn_diff_real_err,
        double **event_Qn_diff_imag, double **event_Qn_diff_imag_err);
    void calculate_Qn_vector_negative_charge(int event_id,
        double *event_Qn_real, double *event_Qn_real_err,
        double *event_Qn_imag, double *event_Qn_imag_err,
        double **event_Qn_diff_real, double **event_Qn_diff_real_err,
        double **event_Qn_diff_imag, double **event_Qn_diff_imag_err);

    //! This function outputs the event averaged particle pT-spectra
    //! and flow coefficients
    void output_Qn_vectors();

    //! This function computes the 2-particle correlation for Qn vectors
    //! within one event
    //!     Real(Qn*conj(Qn)) for n = 0, 1, ... , order_max
    //!     Real(Qn(pT)*conj(Qn)) for n = 0, 1, ... , order_max
    //! self correlation is subtracted assuming full overlap
    void calculate_two_particle_correlation(
            double *event_Qn_real, double *event_Qn_imag,
            double **event_Qn_diff_real, double **event_Qn_diff_imag);

    //! This function outputs the event averaged two-particle flow correlation
    void output_two_particle_correlation();

    //! This function computes the 3-particle correlation for Qn vectors
    //! within one event
    //!     C_nmk = Real(Q_n*Q_m*conj(Q_k)) for (112), (123), (224), (235)
    //! self correlation is subtracted assuming Qk's sample >= Qn's and Qm's
    //! flag = 0: Qn = Qm <= Qk, flag = 1: Qn != Qm, Qn \in Qk, Qm \in Qk
    //! flag = 2: no overlap
    void calculate_three_particle_correlation(
            double *event_Q1_real, double *event_Q1_imag,
            double *event_Q2_real, double *event_Q2_imag,
            double *event_Q3_real, double *event_Q3_imag, int flag,
            double *corr, double *corr_err);

    //! This function outputs the event averaged three-particle correlation
    void output_three_particle_correlation();

    //! this function computes the pT-integrated Qn vector as a function of
    //! rapidity in one event
    void calculate_rapidity_distribution(int event_id,
                    double **event_Qn_real, double **event_Qn_real_err,
                    double **event_Qn_imag, double **event_Qn_imag_err);

    void output_rapidity_distribution();
    
    void check_dNdSV(int event_id);
    void output_dNdSV();

};

#endif  // SRC_single_particleSpectra_h_
