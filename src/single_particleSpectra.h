#ifndef SRC_single_particleSpectra_h_
#define SRC_single_particleSpectra_h_

#include <string>
#include <vector>
#include <memory>

#include "ParameterReader.h"
#include "particleSamples.h"
#include "Random.h"
#include "pretty_ostream.h"

class singleParticleSpectra {
 private:
    const std::string path_;
    std::shared_ptr<particleSamples> particle_list;

    std::shared_ptr<RandomUtil::Random> ran_gen_ptr;

    pretty_ostream messager;

    int particle_monval;

    double reconst_branching_ratio;

    int order_max;

    int npT;
    double pT_min, pT_max, dpT;
    std::vector<double> pT_array;
    std::vector<double> pT_mean_array;
    std::vector<double> pT_mean_array_err;

    int rap_type;
    double rap_min, rap_max;
    int rapidity_distribution_flag;
    int N_rap;
    double rapidity_dis_min, rapidity_dis_max, drap;
    std::vector<double> rapidity_array;
    std::vector<double> dNdy_array;

    double vn_rapidity_dis_pT_min, vn_rapidity_dis_pT_max;
    double **vn_real_rapidity_dis_array, **vn_imag_rapidity_dis_array;
    double **vn_real_rapidity_dis_array_err;
    double **vn_imag_rapidity_dis_array_err;

    int total_number_of_events;
    std::vector<double> Qn_vector_real;
    std::vector<double> Qn_vector_imag;
    std::vector<double> Qn_vector_real_err;
    std::vector<double> Qn_vector_imag_err;
    double **Qn_diff_vector_real, **Qn_diff_vector_imag;
    double **Qn_diff_vector_real_err, **Qn_diff_vector_imag_err;

    int flag_correlation;
    std::vector<double> Qn2_vector;
    std::vector<double> Qn2_vector_err;
    double **QnSP_diff_vector, **QnSP_diff_vector_err;
    double **QnSP_eta12, **QnSP_eta12_err;

    // charge dependent two-particle correlation
    std::vector<double> Cn2_ss;
    std::vector<double> Cn2_ss_err;
    std::vector<double> Cn2_os;
    std::vector<double> Cn2_os_err;
    // \Delta eta dependent two-particle correlation
    std::vector<std::vector<double>> Cn2_ss_eta12;
    std::vector<std::vector<double>> Cn2_ss_eta12_err;
    std::vector<std::vector<double>> Cn2_os_eta12;
    std::vector<std::vector<double>> Cn2_os_eta12_err;

    // 3-particle correlations
    int num_corr;
    std::vector<double> C_nmk;
    std::vector<double> C_nmk_err;
    std::vector<std::vector<double>>C_nmk_eta12;
    std::vector<std::vector<double>>C_nmk_eta13;
    std::vector<std::vector<double>>C_nmk_eta12_err;
    std::vector<std::vector<double>>C_nmk_eta13_err;
    int flag_charge_dependence;
    std::vector<double> C_nmk_ss;
    std::vector<double> C_nmk_ss_err;
    std::vector<double> C_nmk_os;
    std::vector<double> C_nmk_os_err;

    // 3-particle correlation with charge dep on 1 and 3 particles
    std::vector<double> C_nmk_ss_13;
    std::vector<double> C_nmk_ss_13_err;
    std::vector<double> C_nmk_os_13;
    std::vector<double> C_nmk_os_13_err;

    // \Delta \eta dependent 3-particle correlators
    std::vector<std::vector<double>> C_nmk_eta12_ss;
    std::vector<std::vector<double>> C_nmk_eta12_os;
    std::vector<std::vector<double>> C_nmk_eta13_ss;
    std::vector<std::vector<double>> C_nmk_eta13_os;
    std::vector<std::vector<double>> C_nmk_eta12_ss_err;
    std::vector<std::vector<double>> C_nmk_eta12_os_err;
    std::vector<std::vector<double>> C_nmk_eta13_ss_err;
    std::vector<std::vector<double>> C_nmk_eta13_os_err;

    // 4-particle symmetric cumulants
    int SC_num_corr;
    std::vector<double> SC_mn;
    std::vector<double> SC_mn_err;

    // 4-particle cumulant C_n{4}
    int num_Cn4;
    std::vector<double> Cn4;
    std::vector<double> Cn4_err;

    int check_spatial_flag;
    int N_tau;
    double tau_min, tau_max, dtau, intrinsic_dtau;
    double *tau_array, *dNdtau_array;
    int N_xpt;
    double spatial_x_min, spatial_x_max, dspatial_x, intrinsic_dx;
    double spatial_r_min, spatial_r_max, dspatial_r;
    double *xpt_array, *dNdx1_array;
    double *ypt_array, *dNdx2_array;
    double *rpt_array, *dNdr_array;
    double **dNdtaudx1_array, **dNdtaudx2_array;
    int N_eta_s;
    double eta_s_min, eta_s_max, deta_s, intrinsic_detas;
    double *eta_s_array, *dNdetas_array;

 public:
    singleParticleSpectra(const ParameterReader &paraRdr, std::string path,
                          std::shared_ptr<RandomUtil::Random> ran_gen);
    ~singleParticleSpectra();

    //! This is a driver function to compute the Qn flow vector
    void calculate_Qn_vector_shell(
                          std::shared_ptr<particleSamples> particle_list_in);

    int get_monval() const {return(particle_monval);}

    //! this function computes the pT-integrated and pT-differential Qn vector
    //! within a given rapidity region in one event
    void calculate_Qn_vector(int event_id,
        double pT_min_selected, double pT_max_selected,
        std::vector<double> &event_pT_mean, std::vector<double> &event_pT_mean_err,
        std::vector<double> &event_Qn_real, std::vector<double> &event_Qn_real_err,
        std::vector<double> &event_Qn_imag, std::vector<double> &event_Qn_imag_err,
        std::vector<std::vector<double>> &event_Qn_diff_real,
        std::vector<std::vector<double>> &event_Qn_diff_real_err,
        std::vector<std::vector<double>> &event_Qn_diff_imag,
        std::vector<std::vector<double>> &event_Qn_diff_imag_err);

    void calculate_Qn_vector_positive_charge(int event_id,
        std::vector<double> &event_Qn_real,
        std::vector<double> &event_Qn_real_err,
        std::vector<double> &event_Qn_imag,
        std::vector<double> &event_Qn_imag_err,
        std::vector<std::vector<double>> &event_Qn_diff_real,
        std::vector<std::vector<double>> &event_Qn_diff_real_err,
        std::vector<std::vector<double>> &event_Qn_diff_imag,
        std::vector<std::vector<double>> &event_Qn_diff_imag_err);

    void calculate_Qn_vector_negative_charge(int event_id,
        std::vector<double> &event_Qn_real,
        std::vector<double> &event_Qn_real_err,
        std::vector<double> &event_Qn_imag,
        std::vector<double> &event_Qn_imag_err,
        std::vector<std::vector<double>> &event_Qn_diff_real,
        std::vector<std::vector<double>> &event_Qn_diff_real_err,
        std::vector<std::vector<double>> &event_Qn_diff_imag,
        std::vector<std::vector<double>> &event_Qn_diff_imag_err);

    void output_spectra_and_Qn_results();

    //! This function outputs the event averaged particle pT-spectra
    //! and flow coefficients
    void output_Qn_vectors();

    //! This function computes the 2-particle correlation for Qn vectors
    //! within one event with same-sign and opposite-sign pairs
    void calculate_two_particle_correlation_charge_dep(
            std::vector<double> &event_Qn_p_real,
            std::vector<double> &event_Qn_p_imag,
            std::vector<double> &event_Qn_m_real,
            std::vector<double> &event_Qn_m_imag,
            std::vector<double> &corr_ss, std::vector<double> &corr_ss_err,
            std::vector<double> &corr_os, std::vector<double> &corr_os_err);

    //! This function computes the 2-particle correlation for Qn vectors
    //! within one event with particle dependence
    //!     Real(Qn*conj(Qn)) for n = 0, 1, ... , order_max
    //! self correlation is subtracted when flag == 0 (full overlap)
    void calculate_two_particle_correlation_charge_base(
        std::vector<double> &event_Q1_real,
        std::vector<double> &event_Q1_imag,
        std::vector<double> &event_Q2_real,
        std::vector<double> &event_Q2_imag, int flag,
        std::vector<double> &corr, std::vector<double> &corr_err);

    //! This function computes the 2-particle correlation for Qn vectors
    //! within one event
    //!     Real(Qn*conj(Qn)) for n = 0, 1, ... , order_max
    //!     Real(Qn(pT)*conj(Qn)) for n = 0, 1, ... , order_max
    //! self correlation is subtracted assuming full overlap
    void calculate_two_particle_correlation(
        std::vector<double> &event_Qn_real, std::vector<double> &event_Qn_imag,
        std::vector<std::vector<double>> &event_Qn_diff_real,
        std::vector<std::vector<double>> &event_Qn_diff_imag);

    //! This function computes the 2-particle correlation for Qn vectors
    //! as a function of \delta \eta within one event
    //!     Real(Qn(eta_1)*conj(Qn(eta_2))) for n = 0, 1, ... , order_max
    //! self correlation is subtracted when eta_1 = eta_2
    void calculate_two_particle_correlation_deltaeta(
        std::vector<std::vector<double>> &event_Qn_diff_real,
        std::vector<std::vector<double>> &event_Qn_diff_imag);

    //! This function computes the 2-particle correlation for Qn vectors
    //! as a function of \delta \eta within one event with charge depenedence
    //!     Real(Qn(eta_1)*conj(Qn(eta_2))) for n = 0, 1, ... , order_max
    //! self correlation is subtracted when eta_1 = eta_2 and flag == 0
    void calculate_two_particle_correlation_deltaeta_chdep(
            std::vector<std::vector<double>> &event_Qn_p_rap_real,
            std::vector<std::vector<double>> &event_Qn_p_rap_imag,
            std::vector<std::vector<double>> &event_Qn_m_rap_real,
            std::vector<std::vector<double>> &event_Qn_m_rap_imag,
            std::vector<std::vector<double>> &Cn2_ss_eta12,
            std::vector<std::vector<double>> &Cn2_ss_eta12_err,
            std::vector<std::vector<double>> &Cn2_os_eta12,
            std::vector<std::vector<double>> &Cn2_os_eta12_err);

    //! This function computes the 2-particle correlation for Qn vectors
    //! as a function of \delta \eta within one event with charge depenedence
    //!     Real(Qn(eta_1)*conj(Qn(eta_2))) for n = 0, 1, ... , order_max
    //! self correlation is subtracted when eta_1 = eta_2 and flag == 0
    void calculate_two_particle_correlation_deltaeta_chdep_base(
        std::vector<std::vector<double>> &event_Q1_diff_real,
        std::vector<std::vector<double>> &event_Q1_diff_imag,
        std::vector<std::vector<double>> &event_Q2_diff_real,
        std::vector<std::vector<double>> &event_Q2_diff_imag, int flag,
        std::vector<std::vector<double>> &corr,
        std::vector<std::vector<double>> &corr_err);

    //! This function outputs the event averaged two-particle flow correlation
    void output_two_particle_correlation();

    //! This function outputs the two-particle flow correlation as a function
    //! of delta eta between the two particles
    void output_two_particle_correlation_rap();

    //! This function computes the 3-particle correlation for Qn vectors
    //! within one event with same-sign and opposite-sign pairs
    void calculate_three_particle_correlation_charge_dep(
        std::vector<double> &event_Qn_p_real,
        std::vector<double> &event_Qn_p_imag,
        std::vector<double> &event_Qn_m_real,
        std::vector<double> &event_Qn_m_imag,
        std::vector<double> &event_Qn_real,
        std::vector<double> &event_Qn_imag,
        std::vector<double> &corr_ss,    std::vector<double> &corr_ss_err,
        std::vector<double> &corr_os,    std::vector<double> &corr_os_err,
        std::vector<double> &corr_ss_13, std::vector<double> &corr_ss_13_err,
        std::vector<double> &corr_os_13, std::vector<double> &corr_os_13_err);

    //! This function computes the 3-particle correlation for Qn vectors
    //! within one event
    //!     C_nmk = Real(Q_n*Q_m*conj(Q_k)) for (112), (123), (224), (235)
    //! self correlation is subtracted assuming Qk's sample >= Qn's and Qm's
    //! flag = 0: Qn = Qm <= Qk, flag = 1: Qn != Qm, Qn \in Qk, Qm \in Qk
    //! flag = 2: no overlap
    void calculate_three_particle_correlation(
        std::vector<double> &event_Q1_real, std::vector<double> &event_Q1_imag,
        std::vector<double> &event_Q2_real, std::vector<double> &event_Q2_imag,
        std::vector<double> &event_Q3_real, std::vector<double> &event_Q3_imag,
        int flag, std::vector<double> &corr, std::vector<double> &corr_err);

    //! This function computes the 3-particle correlation as a function of 
    //! relative rapidity between particle 1 and particle 2 within one event
    //!     C_nmk(eta_12) = Real(Q_n(eta1)*Q_m(eta2)*conj(Q_k))
    //! self correlation is subtracted assuming Qk's sample >= Qn's and Qm's
    //! Qn != Qm, Qn \in Qk, Qm \in Qk
    //! flag == 1 : eta_12, flag == 2: eta_13
    //! flag_ch == 1: same charge, flag_ch == 2: opposite charge
    void calculate_three_particle_correlation_deltaeta(
        std::vector<std::vector<double>> &event_Q1_real,
        std::vector<std::vector<double>> &event_Q1_imag,
        std::vector<std::vector<double>> &event_Q2_real,
        std::vector<std::vector<double>> &event_Q2_imag,
        std::vector<std::vector<double>> &event_Q3_real,
        std::vector<std::vector<double>> &event_Q3_imag, int flag, int flag_ch,
        std::vector<std::vector<double>> &corr_rap,
        std::vector<std::vector<double>> &corr_rap_err);


    //! This function computes the charge dependent 3-particle correlations
    void calculate_three_particle_correlation_deltaeta_chdep(
            std::vector<std::vector<double>> &event_Qn_p_rap_real,
            std::vector<std::vector<double>> &event_Qn_p_rap_imag,
            std::vector<std::vector<double>> &event_Qn_m_rap_real,
            std::vector<std::vector<double>> &event_Qn_m_rap_imag,
            std::vector<std::vector<double>> &event_Qn_rap_real,
            std::vector<std::vector<double>> &event_Qn_rap_imag,
            std::vector<std::vector<double>> &Cmnk_ss_eta12,
            std::vector<std::vector<double>> &Cmnk_ss_eta12_err,
            std::vector<std::vector<double>> &Cmnk_os_eta12,
            std::vector<std::vector<double>> &Cmnk_os_eta12_err,
            std::vector<std::vector<double>> &Cmnk_ss_eta13,
            std::vector<std::vector<double>> &Cmnk_ss_eta13_err,
            std::vector<std::vector<double>> &Cmnk_os_eta13,
            std::vector<std::vector<double>> &Cmnk_os_eta13_err);

    //! This function outputs the event averaged three-particle correlation
    void output_three_particle_correlation();

    //! This function outputs the rapidity dependent three-particle correlation
    void output_three_particle_correlation_rap();


    //! This function computes the 4-particle correlation for C_n{4}
    //! using Qn vectors within one event
    //!     C_n{4} = <Q1_n*conj(Q2_n)*Q3_n*conj(Q4_n)>
    //! for n = 0, 1, 2, 3, 4
    //! self correlation is subtracted
    void calculate_four_particle_correlation_Cn4(
        std::vector<double> &event_Qn_real, std::vector<double> &event_Qn_imag,
        std::vector<double> &corr, std::vector<double> &corr_err);

    //! This function computes the 4-particle correlation for
    //! symmetric cumulants using Qn vectors within one event
    //!     SC_mn = <Q1_m*conj(Q2_m)*Q3_n*conj(Q4_n)>
    //! for (32), (42), (52), (43), (53)
    //! self correlation is subtracted assuming Qk's sample >= Qn's and Qm's
    void calculate_four_particle_correlation_SC(
        std::vector<double> &event_Q1_real, std::vector<double> &event_Q1_imag,
        std::vector<double> &event_Q2_real, std::vector<double> &event_Q2_imag,
        std::vector<double> &event_Q3_real, std::vector<double> &event_Q3_imag,
        std::vector<double> &event_Q4_real, std::vector<double> &event_Q4_imag,
        int flag, std::vector<double> &corr, std::vector<double> &corr_err);


    //! This function outputs the event averaged four-particle Cn{4}
    void output_four_particle_Cn4_correlation();

    //! This function outputs the event averaged four-particle correlation
    void output_four_particle_SC_correlation();

    //! this function computes the pT-integrated Qn vector as a function of
    //! rapidity in one event
    void calculate_rapidity_distribution(int event_id,
        std::vector<std::vector<double>> &event_Qn_real,
        std::vector<std::vector<double>> &event_Qn_real_err,
        std::vector<std::vector<double>> &event_Qn_imag,
        std::vector<std::vector<double>> &event_Qn_imag_err, int flag);

    void output_rapidity_distribution();

    void check_dNdSV(int event_id);
    void output_dNdSV();

};

#endif  // SRC_single_particleSpectra_h_
