// Copyright 2018 @ Chun Shen

#ifndef BALANCEFUNCTION_H_
#define BALANCEFUNCTION_H_

#include <string>
#include <vector>
#include <array>

#include "ParameterReader.h"
#include "particleSamples.h"

class BalanceFunction {
 private:
    const ParameterReader *paraRdr;
    const std::string path;
    particleSamples *particle_list;

    int particle_monval_a;
    int particle_monval_b;
    bool same_species;

    int Bnpts;
    double Brap_min;
    double Brap_max;
    double drap;
    std::vector<double> Delta_y;
    std::vector<double> N_ab;
    std::vector<double> N_abarbbar;
    std::vector<double> N_abbar;
    std::vector<double> N_abarb;
    std::vector<double> B_delta_y;

 public:
    BalanceFunction(const ParameterReader *paraRdr_in,
                    const std::string path_in,
                    particleSamples *particle_list_in);
    ~BalanceFunction() {};

    void calculate_balance_function();
    void combine_and_bin_particle_pairs(
                std::vector<double> &hist,
                const std::vector< std::vector<particle_info>* >* plist_a,
                const std::vector< std::vector<particle_info>* >* plist_b);
    void combine_and_bin_particle_pairs1(
                std::vector<double> &hist,
                const std::vector< std::vector<particle_info>* >* plist_a,
                const std::vector< std::vector<particle_info>* >* plist_b);
    int get_number_of_particles(
                const std::vector< std::vector<particle_info>* >* plist_b);
    void output_balance_function();
};

#endif  // BALANCEFUNCTION_H_
