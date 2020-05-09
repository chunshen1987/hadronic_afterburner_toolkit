#ifndef HBT_correlation_h
#define HBT_correlation_h

#include <string>
#include <memory>
#include <vector>

#include "ParameterReader.h"
#include "particleSamples.h"
#include "Random.h"
#include "arsenal.h"
#include "pretty_ostream.h"

class HBT_correlation {
 private:
    const ParameterReader paraRdr_;
    const std::string path_;
    std::shared_ptr<particleSamples> particle_list;

    std::shared_ptr<RandomUtil::Random> ran_gen_ptr_;

    pretty_ostream messager;

    int qnpts;
    double q_min, q_max, delta_q;

    int azimuthal_flag_;
    int invariant_radius_flag_;
    double psi_ref;
    int n_KT, n_Kphi;
    double dKT, dKphi;
    double KT_min, KT_max;
    double Krap_min, Krap_max;
    double buffer_rapidity;
    std::vector<double> KT_array_, Kphi_array_;

    std::vector<double> q_out, q_side, q_long;
    int number_of_mixed_events_;
    int number_of_oversample_events_;
    unsigned long long int needed_number_of_pairs;
    unsigned long long int number_pairs_num, number_pairs_denorm;
    std::vector<unsigned long long int> number_of_pairs_numerator_KTdiff;
    std::vector<unsigned long long int> number_of_pairs_denormenator_KTdiff;
    unsigned long long int **number_of_pairs_numerator_KTKphidiff;
    unsigned long long int **number_of_pairs_denormenator_KTKphidiff;

    // arrays for the invariant radius
    double **q_inv_mean;
    double **correl_1d_inv_num, **correl_1d_inv_denorm;
    double **correl_1d_inv_num_count;

    double ****q_out_mean, ****q_side_mean, ****q_long_mean;
    double ****correl_3d_num, ****correl_3d_denorm;
    double ****correl_3d_num_count;
    double *****q_out_diff_mean, *****q_side_diff_mean;
    double *****q_long_diff_mean;
    double *****correl_3d_Kphi_diff_num, *****correl_3d_Kphi_diff_denorm;
    double *****correl_3d_Kphi_diff_num_count;

 public:
    HBT_correlation(ParameterReader &paraRdr, std::string path,
                    std::shared_ptr<RandomUtil::Random> ran_gen);
    ~HBT_correlation();

    double get_psi_ref() {return(psi_ref);};

    void set_particle_list(std::shared_ptr<particleSamples> particle_list_in) {
        particle_list = particle_list_in;
    }

    void calculate_flow_event_plane_angle(int n_order);
    void calculate_HBT_correlation_function(
                std::shared_ptr<particleSamples> particle_list_in);
    void combine_and_bin_particle_pairs(std::vector<int> event_list);
    void combine_and_bin_particle_pairs_mixed_events(
                            int event_id, std::vector<int> mixed_event_list);

    void output_HBTcorrelation();
    void output_correlation_function_inv();
    void output_correlation_function();
    void output_correlation_function_Kphi_differential();

    template <typename T>
    void create_a_2D_array(T **&arr2D, int nx, int ny);
    template <typename T>
    void create_a_3D_array(T ***&arr3D, int nx, int ny, int nz);
    template <typename T>
    void create_a_4D_array(T ****&arr4D, int n1, int n2, int n3, int n4);
    template <typename T>
    void create_a_5D_array(T *****&arr5D,
                           int n1, int n2, int n3, int n4, int n5);
    template <typename T>
    void delete_a_2D_array(T **&arr2D, int nx);
    template <typename T>
    void delete_a_3D_array(T ***&arr3D, int nx, int ny);
    template <typename T>
    void delete_a_4D_array(T ****&arr4D, int n1, int n2, int n3);
    template <typename T>
    void delete_a_5D_array(T *****&arr5D, int n1, int n2, int n3, int n4);
};

#endif
