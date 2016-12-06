// Copyright Chun Shen @ 2016
#ifndef SRC_particle_decay_h_
#define SRC_particle_decay_h_

#include <vector>
#include <cstring>
#include "./particle_info.h"

typedef struct {
    int decay_Npart;
    double branching_ratio;
    int decay_part[5];
} decay_channel_info;

typedef struct {
    int monval;     // Monte Carlo number according PDG
    std::string name;
    double mass;
    double width;
    int gspin;      // spin degeneracy
    int baryon;
    int strange;
    int charm;
    int bottom;
    int gisospin;   // isospin degeneracy
    int charge;
    int decays;     // amount of decays listed for this resonance
    int stable;     // defines whether this particle is considered as stable
    std::vector<decay_channel_info*> decay_channels;
    int sign;                   // Bose-Einstein or Dirac-Fermi statistics
} particle_decay_info;

class particle_decay {
 private:
     std::vector<particle_decay_info*> resonance_table;

 public:
    particle_decay();
    ~particle_decay();
    int read_resonances_list();
    void check_resonance_table();
    double get_particle_width(particle_info *part);
    void perform_two_body_decay(particle_info *mother,
                                particle_info *daughter1,
                                particle_info *daughter2);
    void perform_three_body_decay(particle_info *mother,
                                  particle_info *daughter1,
                                  particle_info *daughter2,
                                  particle_info *daughter3);
    double sample_breit_wigner(double mass, double width, double M_min);
};

#endif  // SRC_particle_decay_h_
