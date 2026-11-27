// Copyright 2018 @ Chun Shen

#include "BalanceFunction.h"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>

using std::cout;
using std::endl;

namespace {
double normalized_pair_ratio(
    double same, double mixed, double total_same, double total_mixed) {
    if (mixed <= 0.0 || total_same <= 0.0 || total_mixed <= 0.0) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    return (same / total_same) / (mixed / total_mixed);
}
}  // namespace

BalanceFunction::BalanceFunction(
    const ParameterReader& paraRdr, const std::string path,
    std::shared_ptr<RandomUtil::Random> ran_gen)
    : paraRdr_(paraRdr), path_(path) {
    ran_gen_ptr = ran_gen;

    particle_monval_a = paraRdr_.getVal("particle_alpha");
    particle_monval_b = paraRdr_.getVal("particle_beta");
    BpT_min = paraRdr_.getVal("BpT_min");
    BpT_max = paraRdr_.getVal("BpT_max");
    Bnpts = paraRdr_.getVal("Bnpts");
    Brap_max = paraRdr_.getVal("Brap_max");
    drap = 2. * std::abs(Brap_max) / (Bnpts - 1);
    Brap_min = -std::abs(Brap_max) - 0.5 * drap;
    Bnphi = 20;
    dphi = 2. * M_PI / Bnphi;
    Bphi_min = -M_PI / 2.;

    rap_type_ = paraRdr_.getVal("rap_type");

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

    N_b = 0;
    N_bbar = 0;
}

void BalanceFunction::calculate_balance_function(
    std::shared_ptr<particleSamples> particle_list_in) {
    set_particle_list(particle_list_in);
    auto plist_a = particle_list->get_balance_function_particle_list_a();
    auto plist_b = particle_list->get_balance_function_particle_list_b();
    auto plist_abar = particle_list->get_balance_function_particle_list_abar();
    auto plist_bbar = particle_list->get_balance_function_particle_list_bbar();
    auto plist_b_mixed_event =
        particle_list->get_balance_function_particle_list_b_mixed_event();
    auto plist_bbar_mixed_event =
        particle_list->get_balance_function_particle_list_bbar_mixed_event();
    if (!plist_a->empty() && plist_b_mixed_event->empty()) {
        throw std::runtime_error(
            "Balance-function analysis requires a mixed-event batch");
    }

    N_b += get_number_of_particles(plist_b);
    N_bbar += get_number_of_particles(plist_bbar);
    messager.info("calculating C_ab ... ");
    combine_and_bin_particle_pairs(C_ab, plist_a, plist_b);
    messager.info("calculating C_abarbbar ... ");
    combine_and_bin_particle_pairs(C_abarbbar, plist_abar, plist_bbar);
    messager.info("calculating C_abbar ... ");
    combine_and_bin_particle_pairs(C_abbar, plist_a, plist_bbar);
    messager.info("calculating C_abarb ... ");
    combine_and_bin_particle_pairs(C_abarb, plist_abar, plist_b);

    messager.info("calculating correlatoin function using mixed events ... ");

    // Without a separate mixed-event stream, the reference lists hold the
    // same events as the signal lists.
    const bool same_batch = !particle_list->uses_external_mixed_events();
    if (same_batch && plist_a->size() == 1 && !warned_single_event_batch_) {
        warned_single_event_batch_ = true;
        messager.warning(
            "Balance function: batches with one event have no independent "
            "mixed event; their mixed-event pairs are skipped. Increase "
            "event_buffer_size or set read_in_real_mixed_events = 1.");
    }
    combine_and_bin_mixed_particle_pairs(
        C_mixed_ab, plist_a, plist_b_mixed_event, same_batch);
    combine_and_bin_mixed_particle_pairs(
        C_mixed_abarbbar, plist_abar, plist_bbar_mixed_event, same_batch);
    combine_and_bin_mixed_particle_pairs(
        C_mixed_abbar, plist_a, plist_bbar_mixed_event, same_batch);
    combine_and_bin_mixed_particle_pairs(
        C_mixed_abarb, plist_abar, plist_b_mixed_event, same_batch);
}

bool BalanceFunction::check_same_particle(
    const particle_info& lhs, const particle_info& rhs) {
    return &lhs == &rhs
           || (lhs.event_track_id != 0
               && lhs.event_track_id == rhs.event_track_id);
}

void BalanceFunction::combine_and_bin_particle_pairs(
    std::vector<std::vector<double>>& hist,
    const std::vector<std::vector<particle_info>*>* plist_a,
    const std::vector<std::vector<particle_info>*>* plist_b) {
    int nev = plist_a->size();
    for (int iev = 0; iev < nev; iev++) {
        for (auto const& part_a : (*(*plist_a)[iev])) {
            if (part_a.pT < BpT_min || part_a.pT > BpT_max) continue;
            for (auto const& part_b : (*(*plist_b)[iev])) {
                if (part_b.pT < BpT_min || part_b.pT > BpT_max) continue;

                if (check_same_particle(part_a, part_b)) continue;

                auto delta_phi_local = part_a.phi_p - part_b.phi_p;
                int phi_idx =
                    ((static_cast<int>(
                         floor((delta_phi_local - Bphi_min) / dphi)))
                     % Bnphi);
                if (phi_idx < 0) phi_idx += Bnphi;

                auto delta_y_local = part_a.rap_y - part_b.rap_y;
                if (rap_type_ == 0)
                    delta_y_local = part_a.rap_eta - part_b.rap_eta;
                if (delta_y_local < Brap_min) continue;

                int y_bin_idx =
                    static_cast<int>((delta_y_local - Brap_min) / drap);
                if (y_bin_idx >= 0 && y_bin_idx < Bnpts) {
                    hist[y_bin_idx][phi_idx] += 1.;
                }
            }
        }
    }
}

void BalanceFunction::combine_and_bin_mixed_particle_pairs(
    std::vector<std::vector<double>>& hist,
    const std::vector<std::vector<particle_info>*>* plist_a,
    const std::vector<std::vector<particle_info>*>* plist_b, bool same_batch) {
    const int nev = plist_a->size();
    const int nev_mixed = plist_b->size();
    if (nev != 0 && nev_mixed == 0) {
        throw std::runtime_error(
            "Balance-function analysis requires a mixed-event batch");
    }
    if (same_batch && nev_mixed != nev) {
        throw std::runtime_error(
            "Balance function: same-batch mixing needs aligned event lists");
    }
    // Pairing an event with itself would add its own correlations, and for
    // like-sign lists every track paired with itself, to the reference.
    if (same_batch && nev_mixed < 2) return;
    for (int iev = 0; iev < nev; iev++) {
        int iev_mixed;
        if (same_batch) {
            // uniform over the nev_mixed - 1 other events
            iev_mixed =
                ran_gen_ptr.lock()->rand_int_uniform() % (nev_mixed - 1);
            if (iev_mixed >= iev) iev_mixed++;
        } else {
            iev_mixed = ran_gen_ptr.lock()->rand_int_uniform() % nev_mixed;
        }
        const double global_random_rotation =
            (ran_gen_ptr.lock()->rand_uniform() * 2. * M_PI);
        for (auto const& part_a : (*(*plist_a)[iev])) {
            if (part_a.pT < BpT_min || part_a.pT > BpT_max) continue;
            for (auto const& part_b : (*(*plist_b)[iev_mixed])) {
                if (part_b.pT < BpT_min || part_b.pT > BpT_max) continue;

                auto delta_phi_local =
                    (part_a.phi_p - part_b.phi_p + global_random_rotation);
                int phi_idx =
                    ((static_cast<int>(
                         floor((delta_phi_local - Bphi_min) / dphi)))
                     % Bnphi);
                if (phi_idx < 0) phi_idx += Bnphi;

                auto delta_y_local = part_a.rap_y - part_b.rap_y;
                if (rap_type_ == 0)
                    delta_y_local = part_a.rap_eta - part_b.rap_eta;
                if (delta_y_local < Brap_min) continue;

                int y_bin_idx =
                    static_cast<int>((delta_y_local - Brap_min) / drap);
                if (y_bin_idx >= 0 && y_bin_idx < Bnpts) {
                    hist[y_bin_idx][phi_idx] += 1.;
                }
            }
        }
    }
}

int BalanceFunction::get_number_of_particles(
    const std::vector<std::vector<particle_info>*>* plist_b) {
    int particle_number = 0;
    for (auto const& ev_i : (*plist_b)) {
        for (auto const& part_b : (*ev_i)) {
            if (part_b.pT < BpT_min || part_b.pT > BpT_max) continue;
            particle_number += 1;
        }
    }
    return (particle_number);
}

void BalanceFunction::output_balance_function() {
    double N_OS = 0.;
    double N_OS_mixed = 0.;
    double N_SS = 0.;
    double N_SS_mixed = 0.;
    std::vector<double> C2_OS_delta_y(Bnpts, 0.);
    std::vector<double> C2_SS_delta_y(Bnpts, 0.);
    std::vector<double> C2_OS_delta_y_mixed(Bnpts, 0.);
    std::vector<double> C2_SS_delta_y_mixed(Bnpts, 0.);
    for (int i = 0; i < Bnpts; i++) {
        for (int j = 0; j < Bnphi; j++) {
            C2_OS_delta_y[i] += C_ab[i][j] + C_abarbbar[i][j];
            C2_SS_delta_y[i] += C_abbar[i][j] + C_abarb[i][j];
            C2_OS_delta_y_mixed[i] += C_mixed_ab[i][j] + C_mixed_abarbbar[i][j];
            C2_SS_delta_y_mixed[i] += C_mixed_abbar[i][j] + C_mixed_abarb[i][j];
        }
        N_OS += C2_OS_delta_y[i];
        N_SS += C2_SS_delta_y[i];
        N_OS_mixed += C2_OS_delta_y_mixed[i];
        N_SS_mixed += C2_SS_delta_y_mixed[i];
    }

    std::vector<double> C2_OS_delta_phi(Bnphi, 0.);
    std::vector<double> C2_SS_delta_phi(Bnphi, 0.);
    std::vector<double> C2_OS_delta_phi_mixed(Bnphi, 0.);
    std::vector<double> C2_SS_delta_phi_mixed(Bnphi, 0.);
    for (int j = 0; j < Bnphi; j++) {
        for (int i = 0; i < Bnpts; i++) {
            C2_OS_delta_phi[j] += C_ab[i][j] + C_abarbbar[i][j];
            C2_SS_delta_phi[j] += C_abbar[i][j] + C_abarb[i][j];
            C2_OS_delta_phi_mixed[j] +=
                C_mixed_ab[i][j] + C_mixed_abarbbar[i][j];
            C2_SS_delta_phi_mixed[j] +=
                C_mixed_abbar[i][j] + C_mixed_abarb[i][j];
        }
    }

    // output the balance function as a function of \Delta y
    std::vector<double> Delta_y(Bnpts, 0.);
    for (int i = 0; i < Bnpts; i++) Delta_y[i] = Brap_min + (i + 0.5) * drap;
    std::ostringstream filename;
    filename << path_ << "/Balance_function_" << particle_monval_a << "_"
             << particle_monval_b << "_Delta_y.dat";
    std::ofstream output(filename.str().c_str(), std::ios::out);
    output << "# DeltaY  Delta_C2  C2(OS)  rho2(OS)  rho1^2(OS)  "
           << "C2(SS) rho2(SS)  rho1^2(SS)" << endl;
    for (int i = 0; i < Bnpts; i++) {
        double C2_OS = normalized_pair_ratio(
            C2_OS_delta_y[i], C2_OS_delta_y_mixed[i], N_OS, N_OS_mixed);
        double C2_SS = normalized_pair_ratio(
            C2_SS_delta_y[i], C2_SS_delta_y_mixed[i], N_SS, N_SS_mixed);
        output << std::scientific << std::setw(18) << std::setprecision(8)
               << Delta_y[i] << "   " << C2_OS - C2_SS << "  " << C2_OS << "  "
               << C2_OS_delta_y[i] << "  " << C2_OS_delta_y_mixed[i] << "  "
               << C2_SS << "  " << C2_SS_delta_y[i] << "  "
               << C2_SS_delta_y_mixed[i] << endl;
    }
    output.close();

    // output the balance function as a function of \Delta phi
    std::vector<double> Delta_phi(Bnphi, 0.);
    for (int j = 0; j < Bnphi; j++) Delta_phi[j] = Bphi_min + (j + 0.5) * dphi;
    std::ostringstream filename2;
    filename2 << path_ << "/Balance_function_" << particle_monval_a << "_"
              << particle_monval_b << "_Delta_phi.dat";
    std::ofstream output2(filename2.str().c_str(), std::ios::out);
    output2 << "# Delta_phi  Delta_C2  C2(OS)  rho2(OS)  rho1^2(OS)  "
            << "C2(SS) rho2(SS)  rho1^2(SS)" << endl;
    for (int j = 0; j < Bnphi; j++) {
        double C2_OS = normalized_pair_ratio(
            C2_OS_delta_phi[j], C2_OS_delta_phi_mixed[j], N_OS, N_OS_mixed);
        double C2_SS = normalized_pair_ratio(
            C2_SS_delta_phi[j], C2_SS_delta_phi_mixed[j], N_SS, N_SS_mixed);
        output2 << std::scientific << std::setw(18) << std::setprecision(8)
                << Delta_phi[j] << "   " << C2_OS - C2_SS << "  " << C2_OS
                << "  " << C2_OS_delta_phi[j] << "  "
                << C2_OS_delta_phi_mixed[j] << "  " << C2_SS << "  "
                << C2_SS_delta_phi[j] << "  " << C2_SS_delta_phi_mixed[j]
                << endl;
    }
    output2.close();

    // output correlation functions as a function of \Delta Y and \Delta phi
    std::ostringstream filename5;
    filename5 << path_ << "/Correlation_function_" << particle_monval_a << "_"
              << particle_monval_b << "_2D.dat";
    std::ofstream output3(filename5.str().c_str(), std::ios::out);
    output3 << "# DY  Dphi  C2(OS)  rho2(OS)  rho1^2(OS)  "
            << "C2(SS)  rho2(SS)  rho1^2(SS)" << endl;
    for (int i = 0; i < Bnpts; i++) {
        for (int j = 0; j < Bnphi; j++) {
            double C2_OS = normalized_pair_ratio(
                C_ab[i][j] + C_abarbbar[i][j],
                C_mixed_ab[i][j] + C_mixed_abarbbar[i][j], N_OS, N_OS_mixed);
            double C2_SS = normalized_pair_ratio(
                C_abbar[i][j] + C_abarb[i][j],
                C_mixed_abbar[i][j] + C_mixed_abarb[i][j], N_SS, N_SS_mixed);
            output3 << std::scientific << std::setw(18) << std::setprecision(8)
                    << Delta_y[i] << "  " << Delta_phi[j] << "  " << C2_OS
                    << "  " << C_ab[i][j] + C_abarbbar[i][j] << "  "
                    << C_mixed_ab[i][j] + C_mixed_abarbbar[i][j] << "  "
                    << C2_SS << "  " << C_abbar[i][j] + C_abarb[i][j] << "  "
                    << C_mixed_abbar[i][j] + C_mixed_abarb[i][j] << endl;
        }
    }
    output3.close();
}
