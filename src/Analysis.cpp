// Copyright Chun Shen @ 2020

#include "Analysis.h"

#include <map>
#include <memory>
#include <string>
#include <vector>

#include "BalanceFunction.h"
#include "HBT_correlation.h"
#include "particle_yield_distribution.h"
#include "single_particleSpectra.h"

using std::vector;

Analysis::Analysis(std::string path) : path_(path) {}

void Analysis::UpdateParameterDict(std::string param_filename) {
    paraRdr_.readFromFile(param_filename);
    paraRdr_.echo();
}

void Analysis::UpdateParameterDict(
    std::string param_filename, int argc, char *argv[]) {
    paraRdr_.readFromFile(param_filename);
    paraRdr_.readFromArguments(argc, argv);

    // rapidity shift from file
    if (paraRdr_.getVal("readRapidityShiftFromFile", 0) == 1) {
        setRapiditySfhitFromFile();
    }
    paraRdr_.echo();
}

void Analysis::setRapiditySfhitFromFile() {
    std::ifstream infile("rapidity_shift.dat");
    if (!infile) {
        messager.warning("Failed to open the file 'rapidity_shift.dat'.");
        return;
    }

    double rapidity_shift, dummy;
    infile >> dummy >> dummy >> rapidity_shift;
    infile.close();
    paraRdr_.setVal("rapidity_shift", rapidity_shift);
}

void Analysis::InitializeAnalysis() {
    int randomSeed = paraRdr_.getVal("randomSeed", -1);
    ran_gen_ptr_ = std::make_shared<RandomUtil::Random>(randomSeed);
}

void Analysis::PerformAnalysis() {
    InitializeAnalysis();
    if (paraRdr_.getVal("analyze_flow") != 0) {
        messager.info("Analyze flow observables ...");
        particle_list_ =
            std::make_shared<particleSamples>(paraRdr_, path_, ran_gen_ptr_);
        if (paraRdr_.getVal("analyze_flow") == 1) {
            FlowAnalysis();
        } else if (paraRdr_.getVal("analyze_flow") == 2) {
            FlowAnalysis_3D();
        } else if (paraRdr_.getVal("analyze_flow") == 3) {
            FlowAnalysis_RHIC();
        } else if (paraRdr_.getVal("analyze_flow") == 4) {
            FlowAnalysis_LHC();
        }
        particle_list_.reset();

        if (paraRdr_.getVal("resonance_weak_feed_down_flag") == 1) {
            // analysis multi-strange particles
            paraRdr_.setVal("resonance_weak_feed_down_flag", 0);
            particle_list_ = std::make_shared<particleSamples>(
                paraRdr_, path_, ran_gen_ptr_);
            if (paraRdr_.getVal("analyze_flow") == 1) {
                FlowAnalysis();
            } else if (paraRdr_.getVal("analyze_flow") == 2) {
                FlowAnalysis_3D();
            } else if (paraRdr_.getVal("analyze_flow") == 3) {
                FlowAnalysis_RHIC();
            } else if (paraRdr_.getVal("analyze_flow") == 4) {
                FlowAnalysis_LHC();
            }
            particle_list_.reset();
        }
    }

    if (paraRdr_.getVal("analyze_HBT") == 1) {
        messager.info("Analyze HBT ...");
        particle_list_ =
            std::make_shared<particleSamples>(paraRdr_, path_, ran_gen_ptr_);
        HBTAnalysis();
        particle_list_.reset();
    }

    if (paraRdr_.getVal("analyze_balance_function") == 1) {
        messager.info("Analyze balance functions ...");
        particle_list_ =
            std::make_shared<particleSamples>(paraRdr_, path_, ran_gen_ptr_);
        BalanceFunctionAnalysis();
        particle_list_.reset();
    }

    if (paraRdr_.getVal("analyze_ebe_yield") == 1) {
        messager.info("Analyze event-by-event particle yield distriution ...");
        particle_list_ =
            std::make_shared<particleSamples>(paraRdr_, path_, ran_gen_ptr_);
        ParticleYieldDistributionAnalysis();
        particle_list_.reset();
    }
}

void Analysis::FlowAnalysis() {
    ParameterReader paraRdr = paraRdr_;
    int compute_correlation = paraRdr.getVal("compute_correlation");
    if (compute_correlation == 1) paraRdr.setVal("compute_correlation", 0);
    int flag_charge_dependence = paraRdr.getVal("flag_charge_dependence");
    if (flag_charge_dependence == 1)
        paraRdr.setVal("flag_charge_dependence", 0);

    // first define all the analysis sets
    std::vector<singleParticleSpectra *> spvn;

    // all hadrons (99999) first (collect E_T)
    paraRdr.setVal("particle_monval", 99999);
    paraRdr.setVal("rap_type", 0);
    paraRdr.setVal("rapidity_distribution", 1);
    paraRdr.setVal("vn_rapidity_dis_pT_min", 0.);
    paraRdr.setVal("vn_rapidity_dis_pT_max", 4.);
    paraRdr.setVal("rapidityPTDistributionFlag", 1);
    spvn.push_back(new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_));

    // charged hadron first
    paraRdr.setVal("particle_monval", 9999);
    paraRdr.setVal("rap_type", 0);
    paraRdr.setVal("rapidity_distribution", 1);
    paraRdr.setVal("rapidityPTDistributionFlag", 0);
    paraRdr.setVal("rap_min", -0.5);
    paraRdr.setVal("rap_max", 0.5);
    paraRdr.setVal("vn_rapidity_dis_pT_min", 0.15);
    paraRdr.setVal("vn_rapidity_dis_pT_max", 2.0);
    spvn.push_back(new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_));
    paraRdr.setVal("vn_rapidity_dis_pT_min", 0.2);
    paraRdr.setVal("vn_rapidity_dis_pT_max", 3.0);
    spvn.push_back(new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_));
    // STAR rapidity cut
    paraRdr.setVal("rap_min", -1.0);
    paraRdr.setVal("rap_max", -0.5);
    spvn.push_back(new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_));
    paraRdr.setVal("rap_min", 0.5);
    paraRdr.setVal("rap_max", 1.0);
    spvn.push_back(new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_));
    paraRdr.setVal("rap_min", -1.0);
    paraRdr.setVal("rap_max", 1.0);
    spvn.push_back(new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_));
    // ALICE rapidity cut
    paraRdr.setVal("rap_min", -0.8);
    paraRdr.setVal("rap_max", 0.8);
    spvn.push_back(new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_));
    // STAR flow rapidity decorrelation cut
    paraRdr.setVal("vn_rapidity_dis_pT_min", 0.4);
    paraRdr.setVal("vn_rapidity_dis_pT_max", 4.0);
    paraRdr.setVal("rap_min", -1.0);
    paraRdr.setVal("rap_max", 1.0);
    spvn.push_back(new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_));
    // CMS default cut
    paraRdr.setVal("vn_rapidity_dis_pT_min", 0.3);
    paraRdr.setVal("vn_rapidity_dis_pT_max", 3.0);
    paraRdr.setVal("rap_min", -2.4);
    paraRdr.setVal("rap_max", 2.4);
    spvn.push_back(new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_));
    // ATLAS UPC cut
    paraRdr.setVal("vn_rapidity_dis_pT_min", 0.4);
    paraRdr.setVal("vn_rapidity_dis_pT_max", 2.0);
    paraRdr.setVal("rap_min", -2.5);
    paraRdr.setVal("rap_max", 2.5);
    spvn.push_back(new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_));
    // ATLAS default cut
    paraRdr.setVal("vn_rapidity_dis_pT_min", 0.5);
    paraRdr.setVal("vn_rapidity_dis_pT_max", 5.0);
    paraRdr.setVal("rap_min", -2.5);
    paraRdr.setVal("rap_max", 2.5);
    spvn.push_back(new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_));
    paraRdr.setVal("rap_min", 0.5);
    paraRdr.setVal("rap_max", 2.5);
    spvn.push_back(new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_));
    paraRdr.setVal("rap_min", -2.5);
    paraRdr.setVal("rap_max", -0.5);
    spvn.push_back(new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_));
    // ATLAS forward
    paraRdr.setVal("rap_min", -4.9);
    paraRdr.setVal("rap_max", -3.1);
    spvn.push_back(new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_));
    paraRdr.setVal("rap_min", 3.1);
    paraRdr.setVal("rap_max", 4.9);
    spvn.push_back(new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_));
    // ALICE V0A
    paraRdr.setVal("vn_rapidity_dis_pT_min", 0.2);
    paraRdr.setVal("vn_rapidity_dis_pT_max", 3.0);
    paraRdr.setVal("rap_min", -5.1);
    paraRdr.setVal("rap_max", -2.8);
    spvn.push_back(new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_));
    paraRdr.setVal("rap_min", 2.8);
    paraRdr.setVal("rap_max", 5.1);
    spvn.push_back(new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_));
    // ALICE V0C
    paraRdr.setVal("rap_min", 1.7);
    paraRdr.setVal("rap_max", 3.7);
    spvn.push_back(new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_));
    paraRdr.setVal("rap_min", -3.7);
    paraRdr.setVal("rap_max", -1.7);
    spvn.push_back(new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_));
    // PHENIX BBC
    paraRdr.setVal("vn_rapidity_dis_pT_min", 0.2);
    paraRdr.setVal("vn_rapidity_dis_pT_max", 2.0);
    paraRdr.setVal("rap_min", -3.9);
    paraRdr.setVal("rap_max", -3.1);
    spvn.push_back(new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_));
    paraRdr.setVal("rap_min", 3.1);
    paraRdr.setVal("rap_max", 3.9);
    spvn.push_back(new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_));

    // now identified particle
    paraRdr.setVal("rap_type", 1);
    paraRdr.setVal("rapidity_distribution", 1);
    paraRdr.setVal("vn_rapidity_dis_pT_min", 0.2);
    paraRdr.setVal("vn_rapidity_dis_pT_max", 3.0);
    paraRdr.setVal("rap_min", -0.5);
    paraRdr.setVal("rap_max", 0.5);
    paraRdr.setVal("particle_monval", 211);
    spvn.push_back(new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_));
    paraRdr.setVal("particle_monval", -211);
    spvn.push_back(new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_));
    paraRdr.setVal("particle_monval", 321);
    spvn.push_back(new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_));
    paraRdr.setVal("particle_monval", -321);
    spvn.push_back(new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_));
    paraRdr.setVal("particle_monval", 2212);
    spvn.push_back(new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_));
    paraRdr.setVal("particle_monval", -2212);
    spvn.push_back(new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_));
    if (paraRdr.getVal("resonance_weak_feed_down_flag") == 0) {
        paraRdr.setVal("particle_monval", 3122);
        spvn.push_back(new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_));
        paraRdr.setVal("particle_monval", -3122);
        spvn.push_back(new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_));
        paraRdr.setVal("particle_monval", 3312);
        spvn.push_back(new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_));
        paraRdr.setVal("particle_monval", -3312);
        spvn.push_back(new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_));
        paraRdr.setVal("particle_monval", 3334);
        spvn.push_back(new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_));
        paraRdr.setVal("particle_monval", -3334);
        spvn.push_back(new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_));
        paraRdr.setVal("particle_monval", 333);
        spvn.push_back(new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_));
    }
    if (paraRdr.getVal("collect_neutral_particles", 0) == 1) {
        paraRdr.setVal("particle_monval", 111);
        spvn.push_back(new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_));
        paraRdr.setVal("particle_monval", 311);
        spvn.push_back(new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_));
        paraRdr.setVal("particle_monval", -311);
        spvn.push_back(new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_));
        paraRdr.setVal("particle_monval", 2112);
        spvn.push_back(new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_));
        paraRdr.setVal("particle_monval", -2112);
        spvn.push_back(new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_));
    }
    paraRdr.setVal("rap_type", 0);
    paraRdr.setVal("particle_monval", 2212);
    spvn.push_back(new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_));
    paraRdr.setVal("particle_monval", -2212);
    spvn.push_back(new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_));

    // lastly, if we want to compute multi-particle correlations within
    // the same UrQMD events
    if (compute_correlation == 1) {
        paraRdr.setVal("compute_correlation", 1);
        paraRdr.setVal("flag_charge_dependence", flag_charge_dependence);
        paraRdr.setVal("particle_monval", 9999);
        paraRdr.setVal("rap_type", 0);
        if (paraRdr.getVal("compute_corr_rap_dep") == 1) {
            // turn on flag to compute the rapidity dependent muti-particle
            // correlations
            paraRdr.setVal("rapidity_distribution", 1);
            paraRdr.setVal("rapidity_dis_min", -2.0);
            paraRdr.setVal("rapidity_dis_max", 2.0);
            paraRdr.setVal("n_rap", 41);
        } else {
            paraRdr.setVal("rapidity_distribution", 0);
        }
        paraRdr.setVal("vn_rapidity_dis_pT_min", 0.2);
        paraRdr.setVal("vn_rapidity_dis_pT_max", 2.0);
        paraRdr.setVal("rap_min", -1.0);
        paraRdr.setVal("rap_max", 1.0);
        spvn.push_back(new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_));
        paraRdr.setVal("rap_min", -2.0);
        paraRdr.setVal("rap_max", 2.0);
        spvn.push_back(new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_));
    }

    // start the loop
    int event_id = 0;
    while (!particle_list_->end_of_file()) {
        messager << "Reading event: " << event_id + 1 << " ... ";
        messager.flush("info");
        int nev = particle_list_->read_in_particle_samples();
        messager << "nev = " << nev;
        messager.flush("info");
        messager.info(" processing ...");
        for (auto &ipart : spvn) {
            int particle_monval = ipart->get_monval();
            particle_list_->filter_particles_from_events(particle_monval);
            ipart->calculate_Qn_vector_shell(particle_list_);
        }
        messager.info("done!");
        event_id += nev;
    }
    for (auto &ipart : spvn) ipart->output_spectra_and_Qn_results();
    spvn.clear();
}

void Analysis::FlowAnalysis_multistrange_particles() {
    ParameterReader paraRdr = paraRdr_;
    paraRdr.setVal("compute_correlation", 0);
    paraRdr.setVal("flag_charge_dependence", 0);
    paraRdr.setVal("rap_type", 1);
    paraRdr.setVal("rapidity_distribution", 1);
    paraRdr.setVal("rapidityPTDistributionFlag", 0);
    paraRdr.setVal("rap_min", -0.5);
    paraRdr.setVal("rap_max", 0.5);
    paraRdr.setVal("vn_rapidity_dis_pT_min", 0.2);
    paraRdr.setVal("vn_rapidity_dis_pT_max", 3.0);

    // first define all the analysis sets
    std::vector<singleParticleSpectra *> spvn;
    paraRdr.setVal("particle_monval", 3122);
    spvn.push_back(new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_));
    paraRdr.setVal("particle_monval", -3122);
    spvn.push_back(new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_));
    paraRdr.setVal("particle_monval", 3312);
    spvn.push_back(new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_));
    paraRdr.setVal("particle_monval", -3312);
    spvn.push_back(new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_));
    paraRdr.setVal("particle_monval", 3334);
    spvn.push_back(new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_));
    paraRdr.setVal("particle_monval", -3334);
    spvn.push_back(new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_));
    paraRdr.setVal("particle_monval", 333);
    spvn.push_back(new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_));

    // start the loop
    int event_id = 0;
    while (!particle_list_->end_of_file()) {
        messager << "Reading event: " << event_id + 1 << " ... ";
        messager.flush("info");
        int nev = particle_list_->read_in_particle_samples();
        messager << "nev = " << nev;
        messager.flush("info");
        messager.info(" processing ...");
        for (auto &ipart : spvn) {
            int particle_monval = ipart->get_monval();
            particle_list_->filter_particles_from_events(particle_monval);
            ipart->calculate_Qn_vector_shell(particle_list_);
        }
        messager.info("done!");
        event_id += nev;
    }
    for (auto &ipart : spvn) ipart->output_spectra_and_Qn_results();
    spvn.clear();
}

void Analysis::FlowAnalysis_3D() {
    ParameterReader paraRdr = paraRdr_;
    paraRdr.setVal("compute_correlation", 0);
    paraRdr.setVal("flag_charge_dependence", 0);

    // first define all the analysis sets
    std::vector<singleParticleSpectra *> spvn;

    // all hadrons (99999) first (collect E_T)
    paraRdr.setVal("particle_monval", 99999);
    paraRdr.setVal("rap_type", 0);
    paraRdr.setVal("rapidity_distribution", 1);
    paraRdr.setVal("vn_rapidity_dis_pT_min", 0.);
    paraRdr.setVal("vn_rapidity_dis_pT_max", 4.);
    paraRdr.setVal("rapidityPTDistributionFlag", 1);
    spvn.push_back(new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_));

    // charged hadron first
    paraRdr.setVal("particle_monval", 9999);
    paraRdr.setVal("rap_type", 0);
    paraRdr.setVal("rapidity_distribution", 1);
    paraRdr.setVal("vn_rapidity_dis_pT_min", 0.2);
    paraRdr.setVal("vn_rapidity_dis_pT_max", 3.0);
    paraRdr.setVal("rapidityPTDistributionFlag", 0);
    paraRdr.setVal("rap_min", -2.5);
    paraRdr.setVal("rap_max", 2.5);
    spvn.push_back(new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_));
    paraRdr.setVal("rap_min", 0.5);
    paraRdr.setVal("rap_max", 2.5);
    spvn.push_back(new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_));
    paraRdr.setVal("rap_min", -2.5);
    paraRdr.setVal("rap_max", -0.5);
    spvn.push_back(new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_));
    paraRdr.setVal("rap_min", -0.5);
    paraRdr.setVal("rap_max", 0.5);
    paraRdr.setVal("rapidityPTDistributionFlag", 1);
    spvn.push_back(new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_));

    // now identified particle
    paraRdr.setVal("rap_type", 1);
    paraRdr.setVal("rapidityPTDistributionFlag", 0);
    paraRdr.setVal("particle_monval", 211);
    spvn.push_back(new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_));
    paraRdr.setVal("particle_monval", -211);
    spvn.push_back(new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_));
    paraRdr.setVal("particle_monval", 321);
    spvn.push_back(new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_));
    paraRdr.setVal("particle_monval", -321);
    spvn.push_back(new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_));
    paraRdr.setVal("particle_monval", 2212);
    spvn.push_back(new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_));
    paraRdr.setVal("particle_monval", -2212);
    spvn.push_back(new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_));
    if (paraRdr.getVal("resonance_weak_feed_down_flag") == 0) {
        paraRdr.setVal("particle_monval", 3122);
        spvn.push_back(new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_));
        paraRdr.setVal("particle_monval", -3122);
        spvn.push_back(new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_));
        paraRdr.setVal("particle_monval", 3312);
        spvn.push_back(new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_));
        paraRdr.setVal("particle_monval", -3312);
        spvn.push_back(new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_));
        paraRdr.setVal("particle_monval", 3334);
        spvn.push_back(new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_));
        paraRdr.setVal("particle_monval", -3334);
        spvn.push_back(new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_));
        paraRdr.setVal("particle_monval", 333);
        spvn.push_back(new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_));
    }

    // start the loop
    int event_id = 0;
    while (!particle_list_->end_of_file()) {
        messager << "Reading event: " << event_id + 1 << " ... ";
        messager.flush("info");
        int nev = particle_list_->read_in_particle_samples();
        messager << "nev = " << nev;
        messager.flush("info");
        messager.info(" processing ...");
        for (auto &ipart : spvn) {
            int particle_monval = ipart->get_monval();
            particle_list_->filter_particles_from_events(particle_monval);
            ipart->calculate_Qn_vector_shell(particle_list_);
        }
        messager.info("done!");
        event_id += nev;
    }
    for (auto &ipart : spvn) ipart->output_spectra_and_Qn_results();
    spvn.clear();
}

void Analysis::FlowAnalysis_RHIC() {
    ParameterReader paraRdr = paraRdr_;
    int compute_correlation = paraRdr.getVal("compute_correlation");
    if (compute_correlation == 1) paraRdr.setVal("compute_correlation", 0);
    int flag_charge_dependence = paraRdr.getVal("flag_charge_dependence");
    if (flag_charge_dependence == 1)
        paraRdr.setVal("flag_charge_dependence", 0);

    // first define all the analysis sets
    std::vector<singleParticleSpectra *> spvn;

    // all hadrons (99999) first (collect E_T)
    paraRdr.setVal("particle_monval", 99999);
    paraRdr.setVal("rap_type", 0);
    paraRdr.setVal("rapidity_distribution", 1);
    paraRdr.setVal("vn_rapidity_dis_pT_min", 0.);
    paraRdr.setVal("vn_rapidity_dis_pT_max", 4.);
    paraRdr.setVal("rapidityPTDistributionFlag", 1);
    spvn.push_back(new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_));

    // charged hadron first
    paraRdr.setVal("particle_monval", 9999);
    paraRdr.setVal("rap_type", 0);
    paraRdr.setVal("rapidity_distribution", 1);
    paraRdr.setVal("vn_rapidity_dis_pT_min", 0.2);
    paraRdr.setVal("vn_rapidity_dis_pT_max", 3.0);
    paraRdr.setVal("rapidityPTDistributionFlag", 1);
    spvn.push_back(new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_));

    // speical triggers follows
    paraRdr.setVal("rapidityPTDistributionFlag", 0);
    paraRdr.setVal("rap_min", -0.5);
    paraRdr.setVal("rap_max", 0.5);
    spvn.push_back(new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_));
    // STAR rapidity cut
    paraRdr.setVal("rap_min", -1.0);
    paraRdr.setVal("rap_max", -0.5);
    spvn.push_back(new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_));
    paraRdr.setVal("rap_min", 0.5);
    paraRdr.setVal("rap_max", 1.0);
    spvn.push_back(new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_));
    // PHENIX BBC
    paraRdr.setVal("rap_min", -3.9);
    paraRdr.setVal("rap_max", -3.1);
    spvn.push_back(new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_));
    paraRdr.setVal("rap_min", 3.1);
    paraRdr.setVal("rap_max", 3.9);
    spvn.push_back(new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_));

    // now identified particle
    paraRdr.setVal("vn_rapidity_dis_pT_min", 0.2);
    paraRdr.setVal("vn_rapidity_dis_pT_max", 3.0);
    if (paraRdr.getVal("pidwithRapidityPTDistribution", 0) == 1) {
        paraRdr.setVal("rapidityPTDistributionFlag", 1);
        paraRdr.setVal("rap_type", 1);
        paraRdr.setVal("particle_monval", 211);
        spvn.push_back(new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_));
        paraRdr.setVal("particle_monval", -211);
        spvn.push_back(new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_));
        paraRdr.setVal("particle_monval", 321);
        spvn.push_back(new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_));
        paraRdr.setVal("particle_monval", -321);
        spvn.push_back(new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_));
        paraRdr.setVal("particle_monval", 2212);
        spvn.push_back(new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_));
        paraRdr.setVal("particle_monval", -2212);
        spvn.push_back(new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_));
    }
    paraRdr.setVal("rap_type", 1);
    paraRdr.setVal("rapidity_distribution", 1);
    paraRdr.setVal("rap_min", -0.5);
    paraRdr.setVal("rap_max", 0.5);
    paraRdr.setVal("particle_monval", 211);
    spvn.push_back(new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_));
    paraRdr.setVal("particle_monval", -211);
    spvn.push_back(new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_));
    paraRdr.setVal("particle_monval", 321);
    spvn.push_back(new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_));
    paraRdr.setVal("particle_monval", -321);
    spvn.push_back(new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_));
    paraRdr.setVal("particle_monval", 2212);
    spvn.push_back(new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_));
    paraRdr.setVal("particle_monval", -2212);
    spvn.push_back(new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_));
    if (paraRdr.getVal("resonance_weak_feed_down_flag") == 0) {
        paraRdr.setVal("particle_monval", 3122);
        spvn.push_back(new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_));
        paraRdr.setVal("particle_monval", -3122);
        spvn.push_back(new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_));
        paraRdr.setVal("particle_monval", 3312);
        spvn.push_back(new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_));
        paraRdr.setVal("particle_monval", -3312);
        spvn.push_back(new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_));
        paraRdr.setVal("particle_monval", 3334);
        spvn.push_back(new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_));
        paraRdr.setVal("particle_monval", -3334);
        spvn.push_back(new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_));
        paraRdr.setVal("particle_monval", 333);
        spvn.push_back(new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_));
    }
    if (paraRdr.getVal("collect_neutral_particles", 0) == 1) {
        paraRdr.setVal("particle_monval", 111);
        spvn.push_back(new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_));
        paraRdr.setVal("particle_monval", 311);
        spvn.push_back(new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_));
        paraRdr.setVal("particle_monval", -311);
        spvn.push_back(new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_));
        paraRdr.setVal("particle_monval", 2112);
        spvn.push_back(new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_));
        paraRdr.setVal("particle_monval", -2112);
        spvn.push_back(new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_));
    }

    // start the loop
    int event_id = 0;
    while (!particle_list_->end_of_file()) {
        messager << "Reading event: " << event_id + 1 << " ... ";
        messager.flush("info");
        int nev = particle_list_->read_in_particle_samples();
        messager << "nev = " << nev;
        messager.flush("info");
        messager.info(" processing ...");
        for (auto &ipart : spvn) {
            int particle_monval = ipart->get_monval();
            particle_list_->filter_particles_from_events(particle_monval);
            ipart->calculate_Qn_vector_shell(particle_list_);
        }
        messager.info("done!");
        event_id += nev;
    }
    for (auto &ipart : spvn) ipart->output_spectra_and_Qn_results();
    spvn.clear();
}

void Analysis::FlowAnalysis_LHC() {
    ParameterReader paraRdr = paraRdr_;
    int compute_correlation = paraRdr.getVal("compute_correlation");
    if (compute_correlation == 1) paraRdr.setVal("compute_correlation", 0);
    int flag_charge_dependence = paraRdr.getVal("flag_charge_dependence");
    if (flag_charge_dependence == 1)
        paraRdr.setVal("flag_charge_dependence", 0);

    // first define all the analysis sets
    std::vector<singleParticleSpectra *> spvn;

    // all hadrons (99999) first (collect E_T)
    paraRdr.setVal("particle_monval", 99999);
    paraRdr.setVal("rap_type", 0);
    paraRdr.setVal("rapidity_distribution", 1);
    paraRdr.setVal("vn_rapidity_dis_pT_min", 0.);
    paraRdr.setVal("vn_rapidity_dis_pT_max", 4.);
    paraRdr.setVal("rapidityPTDistributionFlag", 1);
    spvn.push_back(new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_));

    // charged hadron first
    paraRdr.setVal("particle_monval", 9999);
    paraRdr.setVal("rap_type", 0);
    paraRdr.setVal("rapidity_distribution", 1);
    paraRdr.setVal("rapidityPTDistributionFlag", 1);
    paraRdr.setVal("vn_rapidity_dis_pT_min", 0.2);
    paraRdr.setVal("vn_rapidity_dis_pT_max", 3.0);
    spvn.push_back(new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_));

    // special trigger follows
    paraRdr.setVal("rapidityPTDistributionFlag", 0);
    paraRdr.setVal("rap_min", -0.5);
    paraRdr.setVal("rap_max", 0.5);
    paraRdr.setVal("vn_rapidity_dis_pT_min", 0.2);
    paraRdr.setVal("vn_rapidity_dis_pT_max", 3.0);
    spvn.push_back(new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_));
    // ALICE rapidity cut
    paraRdr.setVal("rap_min", -0.8);
    paraRdr.setVal("rap_max", -0.4);
    spvn.push_back(new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_));
    paraRdr.setVal("rap_min", 0.4);
    paraRdr.setVal("rap_max", 0.8);
    spvn.push_back(new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_));
    paraRdr.setVal("rap_min", -0.4);
    paraRdr.setVal("rap_max", 0.4);
    spvn.push_back(new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_));
    // ALICE V0A
    paraRdr.setVal("rap_min", -5.1);
    paraRdr.setVal("rap_max", -2.8);
    spvn.push_back(new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_));
    paraRdr.setVal("rap_min", 2.8);
    paraRdr.setVal("rap_max", 5.1);
    spvn.push_back(new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_));
    // ALICE V0C
    paraRdr.setVal("rap_min", 1.7);
    paraRdr.setVal("rap_max", 3.7);
    spvn.push_back(new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_));
    paraRdr.setVal("rap_min", -3.7);
    paraRdr.setVal("rap_max", -1.7);
    spvn.push_back(new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_));
    // CMS default cut
    paraRdr.setVal("vn_rapidity_dis_pT_min", 0.3);
    paraRdr.setVal("vn_rapidity_dis_pT_max", 3.0);
    paraRdr.setVal("rap_min", -2.4);
    paraRdr.setVal("rap_max", -0.5);
    spvn.push_back(new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_));
    paraRdr.setVal("rap_min", 0.5);
    paraRdr.setVal("rap_max", 2.4);
    spvn.push_back(new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_));
    // ATLAS default cut
    paraRdr.setVal("vn_rapidity_dis_pT_min", 0.5);
    paraRdr.setVal("vn_rapidity_dis_pT_max", 5.0);
    paraRdr.setVal("rap_min", 0.5);
    paraRdr.setVal("rap_max", 2.5);
    spvn.push_back(new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_));
    paraRdr.setVal("rap_min", -2.5);
    paraRdr.setVal("rap_max", -0.5);
    spvn.push_back(new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_));
    // ATLAS forward
    paraRdr.setVal("rap_min", -4.9);
    paraRdr.setVal("rap_max", -3.1);
    spvn.push_back(new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_));
    paraRdr.setVal("rap_min", 3.1);
    paraRdr.setVal("rap_max", 4.9);
    spvn.push_back(new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_));

    // now identified particle
    paraRdr.setVal("vn_rapidity_dis_pT_min", 0.2);
    paraRdr.setVal("vn_rapidity_dis_pT_max", 3.0);
    if (paraRdr.getVal("pidwithRapidityPTDistribution", 0) == 1) {
        paraRdr.setVal("rapidityPTDistributionFlag", 1);
        paraRdr.setVal("rap_type", 1);
        paraRdr.setVal("particle_monval", 211);
        spvn.push_back(new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_));
        paraRdr.setVal("particle_monval", -211);
        spvn.push_back(new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_));
        paraRdr.setVal("particle_monval", 321);
        spvn.push_back(new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_));
        paraRdr.setVal("particle_monval", -321);
        spvn.push_back(new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_));
        paraRdr.setVal("particle_monval", 2212);
        spvn.push_back(new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_));
        paraRdr.setVal("particle_monval", -2212);
        spvn.push_back(new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_));
    }
    paraRdr.setVal("rapidityPTDistributionFlag", 0);
    if (paraRdr.getVal("pidwithPseudoRapCuts", 0) == 1) {
        paraRdr.setVal("rap_type", 0);
        paraRdr.setVal("rapidity_distribution", 1);
        paraRdr.setVal("rap_min", -0.8);
        paraRdr.setVal("rap_max", 0.8);
        paraRdr.setVal("particle_monval", 211);
        spvn.push_back(new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_));
        paraRdr.setVal("particle_monval", -211);
        spvn.push_back(new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_));
        paraRdr.setVal("particle_monval", 321);
        spvn.push_back(new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_));
        paraRdr.setVal("particle_monval", -321);
        spvn.push_back(new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_));
        paraRdr.setVal("particle_monval", 2212);
        spvn.push_back(new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_));
        paraRdr.setVal("particle_monval", -2212);
        spvn.push_back(new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_));
    }
    paraRdr.setVal("rap_type", 1);
    paraRdr.setVal("rapidity_distribution", 1);
    paraRdr.setVal("rap_min", -0.5);
    paraRdr.setVal("rap_max", 0.5);
    paraRdr.setVal("particle_monval", 211);
    spvn.push_back(new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_));
    paraRdr.setVal("particle_monval", -211);
    spvn.push_back(new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_));
    paraRdr.setVal("particle_monval", 321);
    spvn.push_back(new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_));
    paraRdr.setVal("particle_monval", -321);
    spvn.push_back(new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_));
    paraRdr.setVal("particle_monval", 2212);
    spvn.push_back(new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_));
    paraRdr.setVal("particle_monval", -2212);
    spvn.push_back(new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_));
    if (paraRdr.getVal("resonance_weak_feed_down_flag") == 0) {
        paraRdr.setVal("particle_monval", 3122);
        spvn.push_back(new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_));
        paraRdr.setVal("particle_monval", -3122);
        spvn.push_back(new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_));
        paraRdr.setVal("particle_monval", 3312);
        spvn.push_back(new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_));
        paraRdr.setVal("particle_monval", -3312);
        spvn.push_back(new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_));
        paraRdr.setVal("particle_monval", 3334);
        spvn.push_back(new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_));
        paraRdr.setVal("particle_monval", -3334);
        spvn.push_back(new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_));
        paraRdr.setVal("particle_monval", 333);
        spvn.push_back(new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_));
    }
    if (paraRdr.getVal("collect_neutral_particles", 0) == 1) {
        paraRdr.setVal("particle_monval", 111);
        spvn.push_back(new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_));
        paraRdr.setVal("particle_monval", 311);
        spvn.push_back(new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_));
        paraRdr.setVal("particle_monval", -311);
        spvn.push_back(new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_));
        paraRdr.setVal("particle_monval", 2112);
        spvn.push_back(new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_));
        paraRdr.setVal("particle_monval", -2112);
        spvn.push_back(new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_));
    }

    // lastly, if we want to compute multi-particle correlations within
    // the same UrQMD events
    if (compute_correlation == 1) {
        paraRdr.setVal("compute_correlation", 1);
        paraRdr.setVal("flag_charge_dependence", flag_charge_dependence);
        paraRdr.setVal("particle_monval", 9999);
        paraRdr.setVal("rap_type", 0);
        if (paraRdr.getVal("compute_corr_rap_dep") == 1) {
            // turn on flag to compute the rapidity dependent muti-particle
            // correlations
            paraRdr.setVal("rapidity_distribution", 1);
            paraRdr.setVal("rapidity_dis_min", -2.0);
            paraRdr.setVal("rapidity_dis_max", 2.0);
            paraRdr.setVal("n_rap", 41);
        } else {
            paraRdr.setVal("rapidity_distribution", 0);
        }
        paraRdr.setVal("vn_rapidity_dis_pT_min", 0.2);
        paraRdr.setVal("vn_rapidity_dis_pT_max", 2.0);
        paraRdr.setVal("rap_min", -1.0);
        paraRdr.setVal("rap_max", 1.0);
        spvn.push_back(new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_));
        paraRdr.setVal("rap_min", -2.0);
        paraRdr.setVal("rap_max", 2.0);
        spvn.push_back(new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_));
    }

    // start the loop
    int event_id = 0;
    while (!particle_list_->end_of_file()) {
        messager << "Reading event: " << event_id + 1 << " ... ";
        messager.flush("info");
        int nev = particle_list_->read_in_particle_samples();
        messager << "nev = " << nev;
        messager.flush("info");
        messager.info(" processing ...");
        for (auto &ipart : spvn) {
            int particle_monval = ipart->get_monval();
            particle_list_->filter_particles_from_events(particle_monval);
            ipart->calculate_Qn_vector_shell(particle_list_);
        }
        messager.info("done!");
        event_id += nev;
    }
    for (auto &ipart : spvn) ipart->output_spectra_and_Qn_results();
    spvn.clear();
}

void Analysis::HBTAnalysis() {
    HBT_correlation HBT_analysis(paraRdr_, path_, ran_gen_ptr_);
    // start the loop
    int event_id = 0;
    while (!particle_list_->end_of_file()) {
        messager << "Reading event: " << event_id + 1 << " ... ";
        messager.flush("info");
        particle_list_->read_in_particle_samples_and_filter();
        particle_list_->read_in_particle_samples_mixed_event_and_filter();
        int nev = particle_list_->get_number_of_events();
        messager << "nev = " << nev;
        messager.flush("info");
        messager.info(" processing ...");
        HBT_analysis.calculate_HBT_correlation_function(particle_list_);
        messager.info("done!");
        event_id += nev;
    }
    HBT_analysis.output_HBTcorrelation();
}

void Analysis::ParticleYieldDistributionAnalysis() {
    particle_yield_distribution partN_dis(paraRdr_, path_);
    // start the loop
    int event_id = 0;
    while (!particle_list_->end_of_file()) {
        messager << "Reading event: " << event_id + 1 << " ... ";
        messager.flush("info");
        particle_list_->read_in_particle_samples_and_filter();
        int nev = particle_list_->get_number_of_events();
        messager << "nev = " << nev;
        messager.flush("info");
        messager.info(" processing ...");
        partN_dis.collect_particle_yield_distribution(particle_list_);
        messager.info("done!");
        event_id += nev;
    }
    partN_dis.output_particle_yield_distribution();
}

void Analysis::BalanceFunctionAnalysis() {
    BalanceFunction BF_analysis(paraRdr_, path_, ran_gen_ptr_);
    // start the loop
    int event_id = 0;
    while (!particle_list_->end_of_file()) {
        messager << "Reading event: " << event_id + 1 << " ... ";
        messager.flush("info");
        particle_list_->read_in_particle_samples_and_filter();
        particle_list_->read_in_particle_samples_mixed_event_and_filter();
        int nev = particle_list_->get_number_of_events();
        messager << "nev = " << nev;
        messager.flush("info");
        messager.info(" processing ...");
        BF_analysis.calculate_balance_function(particle_list_);
        messager.info("done!");
        event_id += nev;
    }
    BF_analysis.output_balance_function();
}
