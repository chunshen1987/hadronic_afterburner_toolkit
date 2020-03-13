// Copyright 2018 @ Chun Shen

#ifndef BALANCEFUNCTION_H_
#define BALANCEFUNCTION_H_

#include <string>
#include <vector>
#include <array>
#include <memory>

#include "ParameterReader.h"
#include "particleSamples.h"
#include "Random.h"
#include "pretty_ostream.h"

class BalanceFunction {
 private:
    const ParameterReader paraRdr_;
    const std::string path_;
    std::shared_ptr<particleSamples> particle_list;

    std::weak_ptr<RandomUtil::Random> ran_gen_ptr;

    pretty_ostream messager;

    int particle_monval_a;
    int particle_monval_b;
    bool same_species;

    long int N_b, N_bbar;

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
    BalanceFunction(const ParameterReader &paraRdr,
                    const std::string path,
                    std::shared_ptr<RandomUtil::Random> ran_gen);
    ~BalanceFunction() {};

    void set_particle_list(std::shared_ptr<particleSamples> particle_list_in) {
        particle_list = particle_list_in;
    }
    bool check_same_particle(const particle_info &lhs,
                             const particle_info &rhs);
    void calculate_balance_function(
                std::shared_ptr<particleSamples> particle_list_in);
    void combine_and_bin_particle_pairs(
                std::vector<std::vector<double>> &hist,
                const std::vector< std::vector<particle_info>* >* plist_a,
                const std::vector< std::vector<particle_info>* >* plist_b);
    void combine_and_bin_mixed_particle_pairs(
                std::vector<std::vector<double>> &hist,
                const std::vector< std::vector<particle_info>* >* plist_a,
                const std::vector< std::vector<particle_info>* >* plist_b);
    int get_number_of_particles(
                const std::vector< std::vector<particle_info>* >* plist_b);
    void output_balance_function();
};

#endif  // BALANCEFUNCTION_H_
