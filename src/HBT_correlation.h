#ifndef HBT_correlation_h
#define HBT_correlation_h

#include "ParameterReader.h"
#include "particleSamples.h"

using namespace std;

struct particle_pair
{
    double K_perp, K_phi, K_rap;
    double q_out, q_side, q_long;
    double cos_qx;
};

class HBT_correlation
{
    private:
        ParameterReader *paraRdr;
        string path;
        particleSamples *particle_list;

        int qnpts;
        double q_min, q_max, delta_q;

        int azimuthal_flag;
        int n_KT, n_Kphi;
        double dKT, dKphi;
        double KT_min, KT_max;
        double Krap_min, Krap_max;
        double buffer_rapidity;
        double *KT_array, *Kphi_array;

        double *q_out, *q_side, *q_long;
        int number_of_mixed_events;
        int number_of_oversample_events;
        unsigned long long int needed_number_of_pairs;
        unsigned long long int number_pairs_num, number_pairs_denorm;
        unsigned long long int *number_of_pairs_numerator_KTdiff, *number_of_pairs_denormenator_KTdiff;
        unsigned long long int **number_of_pairs_numerator_KTKphidiff, **number_of_pairs_denormenator_KTKphidiff;
        double ****q_out_mean, ****q_side_mean, ****q_long_mean;
        double ****correl_3d_num, ****correl_3d_denorm;
        double ****correl_3d_num_count;
        double *****q_out_diff_mean, *****q_side_diff_mean, *****q_long_diff_mean;
        double *****correl_3d_Kphi_diff_num, *****correl_3d_Kphi_diff_denorm;
        double *****correl_3d_Kphi_diff_num_count;

    public:
        HBT_correlation(ParameterReader* paraRdr_in, string path_in, particleSamples *particle_list_in);
        ~HBT_correlation();

        void calculate_HBT_correlation_function();
        void combine_and_bin_particle_pairs(int* event_list);
        void combine_and_bin_particle_pairs_mixed_events(int event_id, int* mixed_event_list);
        void bin_into_correlation_function(int type, int num_pair, particle_pair* pairlist);

        void output_correlation_function();
        void output_correlation_function_Kphi_differential();

};

#endif
