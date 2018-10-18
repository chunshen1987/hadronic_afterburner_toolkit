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
    Brap_max = paraRdr->getVal("Brap_max");
    Brap_min = -Brap_max;
    drap     = (Brap_max - Brap_min)/(Bnpts - 1);
    Bnphi    = 20;
    dphi     = 2.*M_PI/Bnphi;
    Bphi_min = -M_PI/2.;

    C_ab.resize(Bnpts);
    C_abarbbar.resize(Bnpts);
    C_abarb.resize(Bnpts);
    C_abbar.resize(Bnpts);
    C_mixed_ab.resize(Bnpts);
    C_mixed_abarbbar.resize(Bnpts);
    C_mixed_abarb.resize(Bnpts);
    C_mixed_abbar.resize(Bnpts);
    for (int i = 0; i < Bnpts; i++) {
        C_ab[i].assign(Bnphi, 0.);
        C_abarbbar[i].assign(Bnphi, 0.);
        C_abarb[i].assign(Bnphi, 0.);
        C_abbar[i].assign(Bnphi, 0.);
        C_mixed_ab[i].assign(Bnphi, 0.);
        C_mixed_abarbbar[i].assign(Bnphi, 0.);
        C_mixed_abarb[i].assign(Bnphi, 0.);
        C_mixed_abbar[i].assign(Bnphi, 0.);
    }

    N_b    = 0;
    N_bbar = 0;
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
        cout << "calculating C_ab ... " << endl;
        combine_and_bin_particle_pairs(C_ab, C_mixed_ab, plist_a, plist_b);
        cout << "calculating C_abarbbar ... " << endl;
        combine_and_bin_particle_pairs(C_abarbbar, C_mixed_abarbbar,
                                       plist_abar, plist_bbar);
        cout << "calculating C_abbar ... " << endl;
        combine_and_bin_particle_pairs(C_abbar, C_mixed_abbar,
                                       plist_a, plist_bbar);
        cout << "calculating C_abarb ... " << endl;
        combine_and_bin_particle_pairs(C_abarb, C_mixed_abarb,
                                       plist_abar, plist_b);
        event_id += buffer_size;
    }
    
    output_balance_function();
}


void BalanceFunction::combine_and_bin_particle_pairs(
                std::vector<std::vector<double>> &hist,
                std::vector<std::vector<double>> &hist_mixed,
                const std::vector< std::vector<particle_info>* >* plist_a,
                const std::vector< std::vector<particle_info>* >* plist_b) {
    int nev = plist_a->size();
    for (int iev = 0; iev < nev; iev++) {
        for (auto const& part_a: (*(*plist_a)[iev])) {
            if (part_a.pT < BpT_min || part_a.pT > BpT_max) continue;
            for (auto const& part_b: (*(*plist_b)[iev])) {
                if (part_b.pT < BpT_min || part_b.pT > BpT_max) continue;

                auto delta_phi_local = part_a.phi_p - part_b.phi_p;
                // randomized the angle to compute the combinatorial background
                auto delta_phi_mixed = delta_phi_local + drand48()*2.*M_PI;
                int phi_idx = ((static_cast<int>(
                            floor((delta_phi_local - Bphi_min)/dphi)))%Bnphi);
                if (phi_idx < 0) phi_idx += Bnphi;
                int phi_mixed_idx = ((static_cast<int>(
                            floor((delta_phi_mixed - Bphi_min)/dphi)))%Bnphi);
                if (phi_mixed_idx < 0) phi_mixed_idx += Bnphi;

                auto delta_y_local = part_a.rap_y - part_b.rap_y;
                if (delta_y_local < Brap_min) continue;
                if (std::abs(delta_y_local) < 1e-15) continue;

                int y_bin_idx = static_cast<int>(
                                            (delta_y_local - Brap_min)/drap);
                if (y_bin_idx >= 0 && y_bin_idx < Bnpts) {
                    hist[y_bin_idx][phi_idx]             += 1.;
                    hist_mixed[y_bin_idx][phi_mixed_idx] += 1.;
                }
            }
        }
    }
}


int BalanceFunction::get_number_of_particles(
                const std::vector< std::vector<particle_info>* >* plist_b) {
    int particle_number = 0;
    for (auto const& ev_i: (*plist_b)) {
        for (auto const& part_b: (*ev_i)) {
            if (part_b.pT < BpT_min || part_b.pT > BpT_max) continue;
            particle_number += 1;
        }
    }
    return(particle_number);
}

void BalanceFunction::output_balance_function() {
    // output the balance function as a function of \Delta y
    std::vector<double> Delta_y(Bnpts, 0.);
    for (int i = 0; i < Bnpts; i++)
        Delta_y[i] = Brap_min + (i + 0.5)*drap;
    std::vector<double> B_delta_y(Bnpts, 0.);
    std::vector<double> B_OS_delta_y(Bnpts, 0.);
    std::vector<double> B_SS_delta_y(Bnpts, 0.);
    for (int i = 0; i < Bnpts; i++) {
        for (int j = 0; j < Bnphi; j++) {
            B_OS_delta_y[i] += C_ab[i][j]    + C_abarbbar[i][j];
            B_SS_delta_y[i] += C_abbar[i][j] + C_abarb[i][j];
            B_delta_y[i]    += B_OS_delta_y[i] - B_SS_delta_y[i];
        }
        B_OS_delta_y[i] /= static_cast<double>(N_b + N_bbar);
        B_SS_delta_y[i] /= static_cast<double>(N_b + N_bbar);
        B_delta_y[i]    /= static_cast<double>(N_b + N_bbar);
    }
    std::ostringstream filename;
    filename << path << "/Balance_function_" << particle_monval_a << "_"
             << particle_monval_b << "_Delta_y.dat";
    std::ofstream output(filename.str().c_str(), std::ios::out);
    output << "# Delta y  B(Delta y)  C_OS(Delta y)  C_SS(Delta y)" << endl;
    for (int i = 0; i < Bnpts; i++) {
        output << std::scientific << std::setw(18) << std::setprecision(8)
               << Delta_y[i] << "   " << B_delta_y[i] << "  "
               << B_OS_delta_y[i] << "  " << B_SS_delta_y[i] << endl;
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
              << particle_monval_b << "_Delta_phi.dat";
    std::ofstream output2(filename2.str().c_str(), std::ios::out);
    for (int j = 0; j < Bnphi; j++) {
        output2 << std::scientific << std::setw(18) << std::setprecision(8)
                << Bphi_min + j*dphi << "   " << B_delta_phi[j] << endl;
    }
    output2.close();

    // output correlation functions as a function of \Delta Y and \Delta phi
    std::ostringstream filename3;
    filename3 << path << "/Balance_function_" << particle_monval_a << "_"
              << particle_monval_b << "_os_2D.dat";
    std::ofstream output3(filename3.str().c_str(), std::ios::out);
    for (int i = 0; i < Bnpts; i++) {
        for (int j = 0; j < Bnphi; j++) {
            output3 << std::scientific << std::setw(18) << std::setprecision(8)
                    << ((C_ab[i][j] + C_abarbbar[i][j])
                        /static_cast<double>(N_b + N_bbar)) << "  ";
        }
        output3 << endl;
    }
    output3.close();
    std::ostringstream filename4;
    filename4 << path << "/Balance_function_" << particle_monval_a << "_"
              << particle_monval_b << "_ss_2D.dat";
    std::ofstream output4(filename4.str().c_str(), std::ios::out);
    for (int i = 0; i < Bnpts; i++) {
        for (int j = 0; j < Bnphi; j++) {
            output4 << std::scientific << std::setw(18) << std::setprecision(8)
                    << ((C_abbar[i][j] + C_abarb[i][j])
                        /static_cast<double>(N_b + N_bbar)) << "  ";
        }
        output4 << endl;
    }
    output4.close();
    
    // output correlation functions as a function of \Delta Y and \Delta phi
    std::ostringstream filename5;
    filename5 << path << "/Correlation_function_" << particle_monval_a << "_"
              << particle_monval_b << "_os_2D.dat";
    std::ofstream output5(filename5.str().c_str(), std::ios::out);
    for (int i = 0; i < Bnpts; i++) {
        for (int j = 0; j < Bnphi; j++) {
            output5 << std::scientific << std::setw(18) << std::setprecision(8)
                    << ((C_ab[i][j] + C_abarbbar[i][j])
                        /(C_mixed_ab[i][j] + C_mixed_abarbbar[i][j])) << "  ";
        }
        output5 << endl;
    }
    output5.close();

    std::ostringstream filename6;
    filename6 << path << "/Correlation_function_" << particle_monval_a << "_"
              << particle_monval_b << "_ss_2D.dat";
    std::ofstream output6(filename6.str().c_str(), std::ios::out);
    for (int i = 0; i < Bnpts; i++) {
        for (int j = 0; j < Bnphi; j++) {
            output6 << std::scientific << std::setw(18) << std::setprecision(8)
                    << ((C_abbar[i][j] + C_abarb[i][j])
                        /(C_mixed_abbar[i][j] + C_mixed_abarb[i][j])) << "  ";
        }
        output6 << endl;
    }
    output6.close();
}
