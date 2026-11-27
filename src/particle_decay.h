// Copyright Chun Shen @ 2016
#ifndef SRC_particle_decay_h_
#define SRC_particle_decay_h_

#include <atomic>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include "Random.h"
#include "particle_info.h"

typedef struct {
    int decay_Npart;
    double branching_ratio;
    int decay_part[5];
    double pole_threshold;
    double minimum_threshold;
    bool has_broad_daughter;
} decay_channel_info;

typedef struct {
    int monval;  // Monte Carlo number according PDG
    std::string name;
    double mass;
    double width;
    int gspin;  // spin degeneracy
    int baryon;
    int strange;
    int charm;
    int bottom;
    int gisospin;  // isospin degeneracy
    int charge;
    int decays;  // amount of decays listed for this resonance
    int stable;  // defines whether this particle is considered as stable
    std::vector<decay_channel_info> decay_channels;
    int sign;  // Bose-Einstein or Dirac-Fermi statistics
    double minimum_mass;
} particle_decay_info;

class particle_decay {
  private:
    std::map<int, particle_decay_info> resonance_table_;
    std::shared_ptr<RandomUtil::Random> ran_gen_ptr_;
    int weak_flag_;
    mutable std::atomic<long> unreachable_resonances_{0};
    std::atomic<long> mass_sampling_fallbacks_{0};

    const particle_decay_info &particle_data(int monval) const;
    void prepare_channel_kinematics();
    double compute_minimum_mass(int monval, std::map<int, int> &state);
    static bool is_broad(const particle_decay_info &particle);
    static bool channel_reachable(
        const decay_channel_info &channel, double mass);
    static bool has_reachable_channel(
        const particle_decay_info &particle, double mass);
    void sample_daughter_masses(
        double mass, std::vector<particle_info> &daughters,
        RandomUtil::Random &ran);
    void perform_n_body_decay(
        const particle_info &mother, std::vector<particle_info> &daughters,
        RandomUtil::Random &ran);

  public:
    particle_decay(std::shared_ptr<RandomUtil::Random> ran_gen, int weak_flag);
    ~particle_decay();

    //! This function reads in resonance decay table
    int read_resonances_list();

    //! This is a test function to check whether the resonance table is
    //! read in correctly
    void check_resonance_table() const;

    //! This function returns particle width in GeV
    double get_particle_width(const particle_info part) const;

    //! A particle with no channel reachable at its current mass is terminal.
    int check_particle_stable(const particle_info part) const;

    //! This function returns the electric charge of particle
    int get_particle_charge(const int monval) const;

    //! This function returns the baryon number of particle
    int get_particle_baryon_number(const int monval) const;

    //! This function returns the strange number of particle
    int get_particle_strange_number(const int monval) const;

    //! This function returns the particle mass for a given particle id
    double get_particle_mass(const int POI_monval) const;

    //! This is a shell function to perform resonance decays
    void perform_decays(
        particle_info &mother, std::vector<particle_info> &daughter_list);
    // Concurrent callers must provide their own random streams.
    void perform_decays(
        particle_info &mother, std::vector<particle_info> &daughter_list,
        RandomUtil::Random &ran);

    //! This function perform two body decay
    void perform_two_body_decay(
        particle_info &mother, particle_info &daughter1,
        particle_info &daughter2);
    void perform_two_body_decay(
        particle_info &mother, particle_info &daughter1,
        particle_info &daughter2, RandomUtil::Random &ran);

    //! This function perform 3 body decays
    void perform_three_body_decay(
        particle_info &mother, particle_info &daughter1,
        particle_info &daughter2, particle_info &daughter3);
    void perform_three_body_decay(
        particle_info &mother, particle_info &daughter1,
        particle_info &daughter2, particle_info &daughter3,
        RandomUtil::Random &ran);

    //! This function sample mother particle mass according to breit-wigner
    //! distribution
    double sample_breit_wigner(
        const double mass, const double width, const double M_min) const;
    double sample_breit_wigner(
        double mass, double width, double M_min, RandomUtil::Random &ran) const;
};

#endif  // SRC_particle_decay_h_
