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
        const ParameterReader &paraRdr, const std::string path,
        std::shared_ptr<RandomUtil::Random> ran_gen) :
    paraRdr_(paraRdr), path_(path) {

    ran_gen_ptr = ran_gen;

    particle_monval_a = paraRdr_.getVal("particle_alpha");
    particle_monval_b = paraRdr_.getVal("particle_beta");
    if (particle_monval_a == - particle_monval_b)
        same_species = true;
    else
        same_species = false;

    BpT_min  = paraRdr_.getVal("BpT_min");
    BpT_max  = paraRdr_.getVal("BpT_max");
    Bnpts    = paraRdr_.getVal("Bnpts");
    Brap_max = paraRdr_.getVal("Brap_max");
    Brap_min = -Brap_max;
    drap     = (Brap_max - Brap_min)/(Bnpts - 1);
    Bnphi    = 20;
    dphi     = 2.*M_PI/Bnphi;
    Bphi_min = -M_PI/2.;
    events   = paraRdr_.getVal("events");     // to normalize single-particle densities 
    rap_type_ = paraRdr_.getVal("rap_type");
    abs_y_max = paraRdr.getVal("abs_y_max"); // single particle density (Norm_b, Norm_bbar rapidity upper limit
    abs_y_min = paraRdr.getVal("abs_y_min"); // lower limit
  
    C_ab.resize(Bnpts);
    C_abarbbar.resize(Bnpts);
    C_abarb.resize(Bnpts);
    C_abbar.resize(Bnpts);
    C_mixed_ab.resize(Bnpts);
    C_mixed_abarbbar.resize(Bnpts);
    C_mixed_abarb.resize(Bnpts);
    C_mixed_abbar.resize(Bnpts);

    // delcaring R2 two-particle corr functions
    R_ab.resize(Bnpts);
    R_abarbbar.resize(Bnpts);
    R_abarb.resize(Bnpts);
    R_abbar.resize(Bnpts);
    
    for (int i = 0; i < Bnpts; i++) {
        C_ab[i].assign(Bnphi, 0.);
        C_abarbbar[i].assign(Bnphi, 0.);
        C_abarb[i].assign(Bnphi, 0.);
        C_abbar[i].assign(Bnphi, 0.);
        C_mixed_ab[i].assign(Bnphi, 0.);
        C_mixed_abarbbar[i].assign(Bnphi, 0.);
        C_mixed_abarb[i].assign(Bnphi, 0.);
        C_mixed_abbar[i].assign(Bnphi, 0.);


	// filling R2's with initial values
        R_ab[i].assign(Bnphi, 0.);
        R_abarbbar[i].assign(Bnphi, 0.);
        R_abarb[i].assign(Bnphi, 0.);
        R_abbar[i].assign(Bnphi, 0.);
    }

    N_b    = 0;
    N_bbar = 0;
    
    // declaring single-particle density counts 
    Normal_b = 0;
    Normal_bbar = 0;
}


void BalanceFunction::calculate_balance_function(
                std::shared_ptr<particleSamples> particle_list_in) {
    set_particle_list(particle_list_in);
    auto plist_a    = particle_list->get_balance_function_particle_list_a();
    auto plist_b    = particle_list->get_balance_function_particle_list_b();
    auto plist_abar = particle_list->get_balance_function_particle_list_abar();
    auto plist_bbar = particle_list->get_balance_function_particle_list_bbar();

    N_b    += get_number_of_particles(plist_b);
    N_bbar += get_number_of_particles(plist_bbar);

    // calling function defined below to get single-particle counts
    Normal_b += get_normalization_balance(plist_b);
    Normal_bbar += get_normalization_balance(plist_bbar);


    messager.info("calculating C_ab TEST CHECK ... ");
    combine_and_bin_particle_pairs(C_ab, plist_a, plist_b);
    messager.info("calculating C_abarbbar ... ");
    combine_and_bin_particle_pairs(C_abarbbar, plist_abar, plist_bbar);
    messager.info("calculating C_abbar ... ");
    combine_and_bin_particle_pairs(C_abbar, plist_a, plist_bbar);
    messager.info("calculating C_abarb ... ");
    combine_and_bin_particle_pairs(C_abarb, plist_abar, plist_b);

    messager.info(
        "calculating correlatoin function using mixed events ... ");

    auto plist_b_mixed_event    = particle_list->get_balance_function_particle_list_b_mixed_event();
    auto plist_bbar_mixed_event = particle_list->get_balance_function_particle_list_bbar_mixed_event();

    combine_and_bin_mixed_particle_pairs(
                    C_mixed_ab, plist_a, plist_b_mixed_event);
    combine_and_bin_mixed_particle_pairs(
                    C_mixed_abarbbar, plist_abar, plist_bbar_mixed_event);
    combine_and_bin_mixed_particle_pairs(
                    C_mixed_abbar, plist_a, plist_bbar_mixed_event);
    combine_and_bin_mixed_particle_pairs(
                    C_mixed_abarb, plist_abar, plist_b_mixed_event);

    // using nested loop to assign values to R's according to eqn 2.16 in thesis 
    for (int i=0; i < Bnpts; i++) {
      for (int j=0; j < Bnphi; j++) {
      
          R_ab[i][j] = (C_ab[i][j] / C_mixed_ab[i][j]) - 1;
          R_abarbbar[i][j] = (C_abarbbar[i][j] / C_mixed_abarbbar[i][j]) - 1;
          R_abarb[i][j] = (C_abarb[i][j] / C_mixed_abarb[i][j]) - 1;
          R_abbar[i][j] = (C_abbar[i][j] / C_mixed_abbar[i][j]) - 1;
        
      }
    
    }
}


bool BalanceFunction::check_same_particle(const particle_info &lhs,
                                          const particle_info &rhs) {
    bool flag = false;
    const double tol = 1e-15;
    if (lhs.monval == rhs.monval) {
        if (std::abs(lhs.E - rhs.E) < tol) {
            if (std::abs(lhs.px - rhs.px) < tol) {
                if (std::abs(lhs.py - rhs.py) < tol) {
                    if (std::abs(lhs.t - rhs.t) < tol) {
                        if (std::abs(lhs.x - rhs.x) < tol) {
                            if (std::abs(lhs.y - rhs.y) < tol) {
                                flag = true;
                            }
                        }
                    }
                }
            }
        }
    }
    return(flag);
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

                if (same_species) {
                    bool flag = check_same_particle(part_a, part_b);
                    if (flag) continue;
                }

                auto delta_phi_local = part_a.phi_p - part_b.phi_p;
                int phi_idx = ((static_cast<int>(
                            floor((delta_phi_local - Bphi_min)/dphi)))%Bnphi);
                if (phi_idx < 0) phi_idx += Bnphi;

                auto delta_y_local = part_a.rap_y - part_b.rap_y;
                if (rap_type_ == 0)
                    delta_y_local = part_a.rap_eta - part_b.rap_eta;
                if (delta_y_local < Brap_min) continue;

                int y_bin_idx = static_cast<int>(
                                            (delta_y_local - Brap_min)/drap);
                if (y_bin_idx >= 0 && y_bin_idx < Bnpts) {
                    hist[y_bin_idx][phi_idx] += 1.;
                }
            }
        }
    }
}


void BalanceFunction::combine_and_bin_mixed_particle_pairs(
                std::vector<std::vector<double>> &hist,
                const std::vector< std::vector<particle_info>* >* plist_a,
                const std::vector< std::vector<particle_info>* >* plist_b) {
    const int nev       = plist_a->size();
    const int nev_mixed = plist_b->size();
    for (int iev = 0; iev < nev; iev++) {
        const int iev_mixed = (
                        ran_gen_ptr.lock()->rand_int_uniform() % nev_mixed);
        const double global_random_rotation = (
                        ran_gen_ptr.lock()->rand_uniform()*2.*M_PI);
        for (auto const& part_a: (*(*plist_a)[iev])) {
            if (part_a.pT < BpT_min || part_a.pT > BpT_max) continue;
            for (auto const& part_b: (*(*plist_b)[iev_mixed])) {
                if (part_b.pT < BpT_min || part_b.pT > BpT_max) continue;

                auto delta_phi_local = (part_a.phi_p - part_b.phi_p
                                        + global_random_rotation);
                int phi_idx = ((static_cast<int>(
                            floor((delta_phi_local - Bphi_min)/dphi)))%Bnphi);
                if (phi_idx < 0) phi_idx += Bnphi;

                auto delta_y_local = part_a.rap_y - part_b.rap_y;
                if (rap_type_ == 0)
                    delta_y_local = part_a.rap_eta - part_b.rap_eta;
                if (delta_y_local < Brap_min) continue;

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
    int particle_number = 0;
    for (auto const& ev_i: (*plist_b)) {
        for (auto const& part_b: (*ev_i)) {
            if (part_b.pT < BpT_min || part_b.pT > BpT_max) continue;
            particle_number += 1;
        }
    }
    return(particle_number);
}


int BalanceFunction::get_normalization_balance(
	       	const std::vector< std::vector<particle_info>* >* plist_b) {
    int normal_number = 0;
    for (auto const& ev_i: (*plist_b)) {
        for (auto const& part_b: (*ev_i)) {
            if (part_b.pT < BpT_min || part_b.pT > BpT_max) continue;
            if (part_b.rap_y < abs_y_min || part_b.rap_y > abs_y_max) continue; // applying rapidity cut to single-particle density count
            normal_number += 1;
        }
    }
    return(normal_number);
}

void BalanceFunction::output_balance_function() {
    double N_OS       = 0.;
    double N_OS_mixed = 0.;
    double N_SS       = 0.;
    double N_SS_mixed = 0.;

    int Norm_b = Normal_b/events; // Norm_b and Norm_bbar are the final counts used in computing the balance function
    int Norm_bbar = Normal_bbar/events;
 
    std::vector<double> C2_OS_delta_y(Bnpts, 0.);
    std::vector<double> C2_SS_delta_y(Bnpts, 0.);
    std::vector<double> C2_OS_delta_y_mixed(Bnpts, 0.);
    std::vector<double> C2_SS_delta_y_mixed(Bnpts, 0.);
    
    // delcaring R2 arrays for y and phi (below) just as with C2's
    std::vector<double> R2_ab_delta_y(Bnpts, 0.);
    std::vector<double> R2_abarbbar_delta_y(Bnpts, 0.);
    std::vector<double> R2_abarb_delta_y(Bnpts, 0.);
    std::vector<double> R2_abbar_delta_y(Bnpts, 0.);
    
    for (int i = 0; i < Bnpts; i++) {
        for (int j = 0; j < Bnphi; j++) {
            C2_OS_delta_y[i]       += C_ab[i][j] + C_abarbbar[i][j];
            C2_SS_delta_y[i]       += C_abbar[i][j] + C_abarb[i][j];
            C2_OS_delta_y_mixed[i] += C_mixed_ab[i][j] + C_mixed_abarbbar[i][j];
            C2_SS_delta_y_mixed[i] += C_mixed_abbar[i][j] + C_mixed_abarb[i][j];
            
            // "integrating out" dependance on phi by applying sum
	    R2_ab_delta_y[i]       += R_ab[i][j];
            R2_abarbbar_delta_y[i] += R_abarbbar[i][j];
            R2_abarb_delta_y[i]    += R_abarb[i][j];
            R2_abbar_delta_y[i]    += R_abbar[i][j];
        }
      
        N_OS += C2_OS_delta_y[i];
        N_SS += C2_SS_delta_y[i];
        N_OS_mixed += C2_OS_delta_y_mixed[i];
        N_SS_mixed += C2_SS_delta_y_mixed[i];
    }
    std::cout << "R2_ab = " << R2_ab_delta_y[2]; // just a test print to screen when I wanted to make sure the program was edited 

    std::vector<double> C2_OS_delta_phi(Bnphi, 0.);
    std::vector<double> C2_SS_delta_phi(Bnphi, 0.);
    std::vector<double> C2_OS_delta_phi_mixed(Bnphi, 0.);
    std::vector<double> C2_SS_delta_phi_mixed(Bnphi, 0.);


    std::vector<double> R2_ab_delta_phi(Bnphi, 0.);
    std::vector<double> R2_abarbbar_delta_phi(Bnphi, 0.);
    std::vector<double> R2_abarb_delta_phi(Bnphi, 0.);
    std::vector<double> R2_abbar_delta_phi(Bnphi, 0.);
    
    for (int j = 0; j < Bnphi; j++) {
        for (int i = 0; i < Bnpts; i++) {
            C2_OS_delta_phi[j]       += C_ab[i][j] + C_abarbbar[i][j];
            C2_SS_delta_phi[j]       += C_abbar[i][j] + C_abarb[i][j];
            C2_OS_delta_phi_mixed[j] += C_mixed_ab[i][j] + C_mixed_abarbbar[i][j];
            C2_SS_delta_phi_mixed[j] += C_mixed_abbar[i][j] + C_mixed_abarb[i][j];
            
            // "integrating out" dependance on y by applying sum 
	    R2_ab_delta_phi[j]       += R_ab[i][j];
            R2_abarbbar_delta_phi[j] += R_abarbbar[i][j];
            R2_abarb_delta_phi[j]    += R_abarb[i][j];
            R2_abbar_delta_phi[j]    += R_abbar[i][j];
        }
    }

    // output the balance function as a function of \Delta y
    std::vector<double> Delta_y(Bnpts, 0.);
    for (int i = 0; i < Bnpts; i++)
        Delta_y[i] = Brap_min + (i + 0.5)*drap;
    std::ostringstream filename;
    filename << path_ << "/Balance_function_" << particle_monval_a << "_"
             << particle_monval_b << "_Delta_y.dat";
    std::ofstream output(filename.str().c_str(), std::ios::out);

    output << "# DeltaY  Delta_C2  C2(OS)  rho2(OS)  rho1^2(OS)  "
           << "C2(SS) rho2(SS)  rho1^2(SS)  R2_ab  R2_abarbbar  R2_abarb  " 
           << "R2_abbar  N_b  N_bbar  B  " << endl; // additional output histograms, note N_b N_bbar is the text on screen but the output below them is according to Norm_b Norm_bbar above
    for (int i = 0; i < Bnpts; i++) {
        double C2_OS = C2_OS_delta_y[i]/C2_OS_delta_y_mixed[i]*N_OS_mixed/N_OS;
        double C2_SS = C2_SS_delta_y[i]/C2_SS_delta_y_mixed[i]*N_SS_mixed/N_SS; 
        
        // computing balance function according to eqn 2.32 in thesis
        double B = 0.5 * ((Norm_b * R2_ab_delta_y[i]) - (Norm_bbar * R2_abbar_delta_y[i]) + (Norm_bbar * R2_abarbbar_delta_y[i]) - (Norm_b * R2_abarb_delta_y[i]));
        
        output << std::scientific << std::setw(38) << std::setprecision(8)
               << Delta_y[i] << "   " << C2_OS - C2_SS << "  "
               << C2_OS << "  " << C2_OS_delta_y[i] << "  " 
               << C2_OS_delta_y_mixed[i] << "  "
               << C2_SS << "  " << C2_SS_delta_y[i] << "  "
               << C2_SS_delta_y_mixed[i] << "  " << R2_ab_delta_y[i] << "  " 
               << R2_abarbbar_delta_y[i] << "  " << R2_abarb_delta_y[i] 
               << "  " << R2_abbar_delta_y[i] << "  " << Norm_b << "  "
	       << Norm_bbar << "  " << B << endl; // outputting additional histograms
    }
    output.close();

    // output the balance function as a function of \Delta phi
    std::vector<double> Delta_phi(Bnphi, 0.);
    for (int j = 0; j < Bnphi; j++)
        Delta_phi[j] = Bphi_min + (j + 0.5)*dphi;
    std::ostringstream filename2;
    filename2 << path_ << "/Balance_function_" << particle_monval_a << "_"
              << particle_monval_b << "_Delta_phi.dat";
    std::ofstream output2(filename2.str().c_str(), std::ios::out);

    output2 << "# Delta_phi  Delta_C2  C2(OS)  rho2(OS)  rho1^2(OS)  "
            << "C2(SS) rho2(SS)  rho1^2(SS)  R2_ab  R2_abarbbar  R2_abarb  "
            << "R_abbar  N_b  N_bbar  B  " << endl;
    for (int j = 0; j < Bnphi; j++) { // summing out dependence on y
        double C2_OS = (C2_OS_delta_phi[j]/C2_OS_delta_phi_mixed[j]
                        *N_OS_mixed/N_OS);
        double C2_SS = (C2_SS_delta_phi[j]/C2_SS_delta_phi_mixed[j]
                        *N_SS_mixed/N_SS);
                        
        double B = 0.5 * ((Norm_b * R2_ab_delta_phi[j]) - (Norm_bbar * R2_abbar_delta_phi[j]) + (Norm_bbar * R2_abarbbar_delta_phi[j]) - (Norm_b * R2_abarb_delta_phi[j]));
        
        output2 << std::scientific << std::setw(18) << std::setprecision(8)
                << Delta_phi[j] << "   " << C2_OS - C2_SS << "  "
                << C2_OS << "  " << C2_OS_delta_phi[j] << "  "
                << C2_OS_delta_phi_mixed[j] << "  "
                << C2_SS << "  " << C2_SS_delta_phi[j] << "  "
                << C2_SS_delta_phi_mixed[j] << "  " << R2_ab_delta_phi[j] << "  " 
                << R2_abarbbar_delta_phi[j] << "  " << R2_abarb_delta_phi[j] 
                << "  " << R2_abbar_delta_phi[j] << "  " << Norm_b << "  " 
                << Norm_bbar << "  " << B << endl;
    
    }
    output2.close();

    // output correlation functions as a function of \Delta Y and \Delta phi
    std::ostringstream filename5;
    filename5 << path_ << "/Correlation_function_" << particle_monval_a << "_"
              << particle_monval_b << "_2D.dat";
    std::ofstream output3(filename5.str().c_str(), std::ios::out);

    output3 << "# DY  Dphi  C2(OS)  rho2(OS)  rho1^2(OS)  "
            << "C2(SS)  rho2(SS)  rho1^2(SS)  R2_ab  R2_abarbbar  R2_abarb " 
            << "R2_abbar  N_b  N_bbar  B  " << endl;
    for (int i = 0; i < Bnpts; i++) { 
        for (int j = 0; j < Bnphi; j++) {
            double C2_OS = ( (C_ab[i][j] + C_abarbbar[i][j])
                            /(C_mixed_ab[i][j] + C_mixed_abarbbar[i][j] + 1e-15));
            double C2_SS = ( (C_abbar[i][j] + C_abarb[i][j])
                            /(C_mixed_abbar[i][j] + C_mixed_abarb[i][j] + 1e-15));
            
            double B = 0.5 * ((Norm_b * R_ab[i][j]) - (Norm_bbar * R_abbar[i][j]) + (Norm_bbar * R_abarbbar[i][j]) - (Norm_b * R_abarb[i][j]));
            
            output3 << std::scientific << std::setw(18) << std::setprecision(8)
                    << Delta_y[i] << "  " << Delta_phi[j] << "  "
                    << C2_OS*N_OS_mixed/N_OS << "  "
                    << C_ab[i][j] + C_abarbbar[i][j] << "  "
                    << C_mixed_ab[i][j] + C_mixed_abarbbar[i][j] << "  "
                    << C2_SS*N_SS_mixed/N_SS << "  "
                    << C_abbar[i][j] + C_abarb[i][j] << "  "
                    << C_mixed_abbar[i][j] + C_mixed_abarb[i][j] << R_ab[i][j] << "  " 
                    << R_abarbbar[i][j] << "  " << R_abarb[i][j] 
                    << "  " << R_abbar[i][j] << "  " << Norm_b << "  " 
                    << Norm_bbar << "  " << B << endl;
	}
    }
    output3.close();
}
