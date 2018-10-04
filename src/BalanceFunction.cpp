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

    BpT_min  = paraRdr->getVal("BpT_min");
    BpT_max  = paraRdr->getVal("BpT_max");
    Bnpts    = paraRdr->getVal("Bnpts");
    Brap_min = 0.0;
    Brap_max = paraRdr->getVal("Brap_max");
    drap     = (Brap_max - Brap_min)/(Bnpts - 1);
    Bnphi    = 20;
    dphi     = 2*M_PI/Bnphi;

    C_ab.resize(Bnpts);
    C_abarbbar.resize(Bnpts);
    C_abarb.resize(Bnpts);
    C_abbar.resize(Bnpts);
    for (int i = 0; i < Bnpts; i++) {
        C_ab[i].assign(Bnphi, 0.);
        C_abarbbar[i].assign(Bnphi, 0.);
        C_abarb[i].assign(Bnphi, 0.);
        C_abbar[i].assign(Bnphi, 0.);
    }
}


void BalanceFunction::calculate_balance_function() {
    int event_id = 0;
    const int buffer_size = particle_list->get_event_buffer_size();

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
        combine_and_bin_particle_pairs(C_ab, plist_a, plist_b);
        cout << "calculating N_abarbbar ... " << endl;
        combine_and_bin_particle_pairs(C_abarbbar, plist_abar, plist_bbar);
        cout << "calculating N_abbar ... " << endl;
        combine_and_bin_particle_pairs(C_abbar, plist_a, plist_bbar);
        cout << "calculating N_abarb ... " << endl;
        combine_and_bin_particle_pairs(C_abarb, plist_abar, plist_b);
        event_id += buffer_size;
    }
    
    output_balance_function();
}


void BalanceFunction::combine_and_bin_particle_pairs(
                std::vector<std::vector<double>> &hist,
                const std::vector< std::vector<particle_info>* >* plist_a,
                const std::vector< std::vector<particle_info>* >* plist_b) {
    int nev = plist_a->size();
    for (int iev = 0; iev < nev; iev++) {
        for (auto const& part_a: (*(*plist_a)[iev])) {
            if (part_a.pT < BpT_min || part_a.pT > BpT_max) continue;
            for (auto const& part_b: (*(*plist_b)[iev])) {
                if (part_b.pT < BpT_min || part_b.pT > BpT_max) continue;
                auto delta_phi_local = part_a.phi_p - part_b.phi_p;
                int phi_idx = (static_cast<int>(delta_phi_local/dphi))%Bnphi;
                auto delta_y_local = std::abs(part_a.rap_y - part_b.rap_y);
                if (delta_y_local < 1e-15) continue;
                int y_bin_idx = static_cast<int>(
                                            (delta_y_local - Brap_min)/drap);
                if (y_bin_idx >= 0 && y_bin_idx < Bnpts) {
                    hist[y_bin_idx][phi_idx] += 1.;
                }
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
    // output the balance function as a function of \Delta y
    std::vector<double> Delta_y(Bnpts, 0.);
    for (int i = 0; i < Bnpts; i++)
        Delta_y[i] = Brap_min + (i + 0.5)*drap;
    std::vector<double> B_delta_y(Bnpts, 0.);
    for (int i = 0; i < Bnpts; i++) {
        for (int j = 0; j < Bnphi; j++) {
            B_delta_y[i] += (  C_ab[i][j]    + C_abarbbar[i][j]
                             - C_abbar[i][j] - C_abarb[i][j]   );
        }
        B_delta_y[i] /= static_cast<double>(N_b + N_bbar);
    }
    std::ostringstream filename;
    filename << path << "/Balance_function_" << particle_monval_a << "_"
             << particle_monval_b << "_delta_y.dat";
    std::ofstream output(filename.str().c_str(), std::ios::out);
    for (int i = 0; i < Bnpts; i++) {
        output << std::scientific << std::setw(18) << std::setprecision(8)
               << Delta_y[i] << "   " << B_delta_y[i] << endl;
    }
    output.close();
    
    // output the balance function as a function of \Delta phi
    std::vector<double> B_delta_phi(Bnphi, 0.);
    for (int j = 0; j < Bnphi; j++) {
        for (int i = 0; i < Bnpts; i++) {
            B_delta_phi[j] += (  C_ab[i][j]    + C_abarbbar[i][j]
                               - C_abbar[i][j] - C_abarb[i][j]   );
        }
        B_delta_phi[j] /= static_cast<double>(N_b + N_bbar);
    }
    std::ostringstream filename2;
    filename2 << path << "/Balance_function_" << particle_monval_a << "_"
              << particle_monval_b << "_delta_phi.dat";
    std::ofstream output2(filename.str().c_str(), std::ios::out);
    for (int j = 0; j < Bnphi; j++) {
        output2 << std::scientific << std::setw(18) << std::setprecision(8)
                << j*dphi << "   " << B_delta_phi[j] << endl;
    }
    output2.close();
}
