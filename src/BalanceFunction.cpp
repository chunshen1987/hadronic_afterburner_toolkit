// Copyright 2018 @ Chun Shen

#include "BalanceFunction.h"

#include <iostream>
#include <cmath>
#include <sstream>
#include <fstream>
#include <iomanip>

using std::cout;
using std::endl;

BalanceFunction::BalanceFunction(
    const ParameterReader *paraRdr_in, const std::string path_in,
    particleSamples *particle_list_in) :
    paraRdr(paraRdr_in), path(path_in) {

    particle_list = particle_list_in;
        
    particle_monval_a = paraRdr->getVal("particle_alpha");
    particle_monval_b = paraRdr->getVal("particle_beta");
    if (particle_monval_a == - particle_monval_b)
        same_species = true;
    else
        same_species = false;

    Bnpts    = paraRdr->getVal("Bnpts");
    Brap_min = 0.0;
    Brap_max = paraRdr->getVal("Brap_max");
    drap = (Brap_max - Brap_min)/(Bnpts - 1);

    Delta_y.assign(Bnpts, 0.);
    for (int i = 0; i < Bnpts; i++)
        Delta_y[i] = Brap_min + (i + 0.5)*drap;
    B_delta_y.assign(Bnpts, 0.);
    N_ab.assign(Bnpts, 0.);
    N_abarbbar.assign(Bnpts, 0.);
    N_abarb.assign(Bnpts, 0.);
    N_abbar.assign(Bnpts, 0.);
}


void BalanceFunction::calculate_balance_function() {
    int event_id = 0;
    const int buffer_size = particle_list->get_event_buffer_size();

    int N_b    = 0;
    int N_bbar = 0;
    while (!particle_list->end_of_file()) {
        cout << "Reading event: " << event_id + 1 << "-" 
             << event_id + buffer_size << " ... " << std::flush;
        particle_list->read_in_particle_samples();
        cout << " processing ..." << endl;
        auto plist_a    = particle_list->get_balance_function_particle_list_a();
        auto plist_b    = particle_list->get_balance_function_particle_list_b();
        auto plist_abar = particle_list->get_balance_function_particle_list_abar();
        auto plist_bbar = particle_list->get_balance_function_particle_list_bbar();
        
        N_b    += get_number_of_particles(plist_b);
        N_bbar += get_number_of_particles(plist_bbar);
        cout << "calculating N_ab ... " << endl;
        combine_and_bin_particle_pairs(N_ab, plist_a, plist_b);
        cout << "calculating N_abarbbar ... " << endl;
        combine_and_bin_particle_pairs(N_abarbbar, plist_abar, plist_bbar);
        cout << "calculating N_abbar ... " << endl;
        combine_and_bin_particle_pairs(N_abbar, plist_a, plist_bbar);
        cout << "calculating N_abarb ... " << endl;
        combine_and_bin_particle_pairs(N_abarb, plist_abar, plist_b);
        event_id += buffer_size;
    }
    
    for (int i = 0; i < Bnpts; i++) {
        B_delta_y[i] = ((N_ab[i] + N_abarbbar[i] - N_abbar[i] - N_abarb[i])
                        /(N_b + N_bbar));
    }
    output_balance_function();
}


void BalanceFunction::combine_and_bin_particle_pairs(
                std::vector<double> &hist,
                const std::vector< std::vector<particle_info>* >* plist_a,
                const std::vector< std::vector<particle_info>* >* plist_b) {
    int nev = plist_a->size();
    for (int iev = 0; iev < nev; iev++) {
        for (auto const& part_a: (*(*plist_a)[iev])) {
            for (auto const& part_b: (*(*plist_b)[iev])) {
                auto delta_y_local = std::abs(part_a.rap_y - part_b.rap_y);
                if (delta_y_local < 1e-15) continue;
                int y_bin_idx = static_cast<int>(
                                            (delta_y_local - Brap_min)/drap);
                if (y_bin_idx >= 0 && y_bin_idx < Bnpts) {
                    hist[y_bin_idx] += 1.;
                }
            }
        }
    }
}


void BalanceFunction::combine_and_bin_particle_pairs1(
                std::vector<double> &hist,
                const std::vector< std::vector<particle_info>* >* plist_a,
                const std::vector< std::vector<particle_info>* >* plist_b) {
    int nev = plist_a->size();
    for (int iev = 0; iev < nev; iev++) {
        int npart = (*plist_a)[iev]->size();
        for (int ipart = 0; ipart < npart; ipart++) {
            auto delta_y_local = std::abs(
                    (*(*plist_a)[iev])[ipart].rap_y - (*(*plist_b)[iev])[ipart].rap_y);
            int y_bin_idx = static_cast<int>((delta_y_local - Brap_min)/drap);
            if (y_bin_idx >= 0 && y_bin_idx < Bnpts) {
                hist[y_bin_idx] += 1.;
            }
        }
    }
}


int BalanceFunction::get_number_of_particles(
                const std::vector< std::vector<particle_info>* >* plist_b) {
    int N_b = 0;
    for (auto const& ev_i: (*plist_b))
        N_b += ev_i->size();
    return(N_b);
}

void BalanceFunction::output_balance_function() {
    std::ostringstream filename;
    filename << path << "/Balance_function_" << particle_monval_a << "_"
             << particle_monval_b << ".dat";
    std::ofstream output(filename.str().c_str(), std::ios::out);
    for (int i = 0; i < Bnpts; i++) {
        output << std::scientific << std::setw(18) << std::setprecision(8)
               << Delta_y[i] << "   " << B_delta_y[i] << endl;
    }
    output.close();
}
