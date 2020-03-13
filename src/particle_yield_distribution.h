#ifndef SRC_particle_yield_distribution_h_
#define SRC_particle_yield_distribution_h_

#include <string>

#include "ParameterReader.h"
#include "particleSamples.h"
#include "pretty_ostream.h"

class particle_yield_distribution {
 private:
    const ParameterReader paraRdr_;
    const std::string path_;
    std::shared_ptr<particleSamples> particle_list;

    pretty_ostream messager;

    int particle_monval;
    int net_particle_flag;

    double reconst_branching_ratio;

    double pT_min, pT_max;

    int rap_type;
    double rap_min, rap_max;

    int total_number_of_events;
    int n_max;
    int *number_of_events;

 public:
    particle_yield_distribution(ParameterReader &paraRdr, std::string path);
    ~particle_yield_distribution();

    void set_particle_list(std::shared_ptr<particleSamples> particle_list_in) {
        particle_list = particle_list_in;
    }
    void collect_particle_yield_distribution(
                        std::shared_ptr<particleSamples> particle_list_in);
    void collect_particle_yield(int event_id);
    void output_particle_yield_distribution();
};

#endif  // SRC_single_particleSpectra_h_
