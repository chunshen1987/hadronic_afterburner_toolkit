// Copyright 2018 @ Chun Shen

#ifndef BALANCEFUNCTION_H_
#define BALANCEFUNCTION_H_

#include <array>
#include <memory>
#include <string>
#include <vector>

#include "ParameterReader.h"
#include "Random.h"
#include "particleSamples.h"
#include "pretty_ostream.h"

class BalanceFunction {
  private:
    const ParameterReader paraRdr_;
    const std::string path_;
    std::shared_ptr<particleSamples> particle_list;

    std::weak_ptr<RandomUtil::Random> ran_gen_ptr;
    bool warned_single_event_batch_ = false;

    pretty_ostream messager;

    int particle_monval_a;
    int particle_monval_b;

    long int N_b, N_bbar;

    int rap_type_;
    int Bnpts;
    int Bnphi;
    double dphi;
    double Bphi_min;
    double Brap_min;
    double Brap_max;
    double drap;
    double BpT_min, BpT_max;
    std::vector<std::vector<double>> C_ab;
    std::vector<std::vector<double>> C_abarbbar;
    std::vector<std::vector<double>> C_abbar;
    std::vector<std::vector<double>> C_abarb;

    std::vector<std::vector<double>> C_mixed_ab;
    std::vector<std::vector<double>> C_mixed_abarbbar;
    std::vector<std::vector<double>> C_mixed_abbar;
    std::vector<std::vector<double>> C_mixed_abarb;

  public:
    BalanceFunction(
        const ParameterReader& paraRdr, const std::string path,
        std::shared_ptr<RandomUtil::Random> ran_gen);
    ~BalanceFunction() {};

    void set_particle_list(std::shared_ptr<particleSamples> particle_list_in) {
        particle_list = particle_list_in;
    }
    bool check_same_particle(
        const particle_info& lhs, const particle_info& rhs);
    void calculate_balance_function(
        std::shared_ptr<particleSamples> particle_list_in);
    void combine_and_bin_particle_pairs(
        std::vector<std::vector<double>>& hist,
        const std::vector<std::vector<particle_info>*>* plist_a,
        const std::vector<std::vector<particle_info>*>* plist_b);
    //! Bin pairs of each event in plist_a with one random event of plist_b.
    //! Set same_batch when plist_b is filtered from the same events as
    //! plist_a; event iev is then never paired with itself.
    void combine_and_bin_mixed_particle_pairs(
        std::vector<std::vector<double>>& hist,
        const std::vector<std::vector<particle_info>*>* plist_a,
        const std::vector<std::vector<particle_info>*>* plist_b,
        bool same_batch = false);
    int get_number_of_particles(
        const std::vector<std::vector<particle_info>*>* plist_b);
    void output_balance_function();
};

#endif  // BALANCEFUNCTION_H_
