// Copyright Chun Shen @ 2020

#include "Analysis.h"

#include <map>
#include <memory>
#include <set>
#include <stdexcept>
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
    const int flow_mode = paraRdr_.getVal("analyze_flow");
    if (flow_mode != 0) {
        messager.info("Analyze flow observables ...");
        perform_flow(paraRdr_, flow_mode);
        if (paraRdr_.getVal("resonance_weak_feed_down_flag") == 1) {
            // Collect multistrange hadrons before their weak decays.
            auto config = paraRdr_;
            config.setVal("resonance_weak_feed_down_flag", 0);
            perform_flow(config, flow_mode);
        }
    }

    if (paraRdr_.getVal("analyze_HBT") == 1) {
        InitializeAnalysis();
        messager.info("Analyze HBT ...");
        particle_list_ =
            std::make_shared<particleSamples>(paraRdr_, path_, ran_gen_ptr_);
        HBTAnalysis();
        particle_list_.reset();
    }

    if (paraRdr_.getVal("analyze_balance_function") == 1) {
        InitializeAnalysis();
        messager.info("Analyze balance functions ...");
        particle_list_ =
            std::make_shared<particleSamples>(paraRdr_, path_, ran_gen_ptr_);
        BalanceFunctionAnalysis();
        particle_list_.reset();
    }

    if (paraRdr_.getVal("analyze_ebe_yield") == 1) {
        InitializeAnalysis();
        messager.info("Analyze event-by-event particle yield distriution ...");
        particle_list_ =
            std::make_shared<particleSamples>(paraRdr_, path_, ran_gen_ptr_);
        ParticleYieldDistributionAnalysis();
        particle_list_.reset();
    }
}

void Analysis::perform_flow(const ParameterReader &config, int mode) {
    InitializeAnalysis();
    switch (mode) {
        case 1:
            FlowAnalysis(config);
            break;
        case 2:
            FlowAnalysis_3D(config);
            break;
        case 3:
            FlowAnalysis_RHIC(config);
            break;
        case 4:
            FlowAnalysis_LHC(config);
            break;
        default:
            throw std::invalid_argument("Invalid analyze_flow mode");
    }
}

void Analysis::run_flow_spectra(
    const ParameterReader &config,
    const vector<std::unique_ptr<singleParticleSpectra>> &spectra) {
    std::map<std::string, std::set<std::string>> claimed_outputs;
    // Preserve the filename of the last matching definition.
    for (auto it = spectra.rbegin(); it != spectra.rend(); ++it) {
        (*it)->plan_outputs(claimed_outputs);
    }

    auto reader_config = config;
    particle_list_ =
        std::make_shared<particleSamples>(reader_config, path_, ran_gen_ptr_);
    int event_id = 0;
    while (!particle_list_->end_of_file()) {
        messager << "Reading event: " << event_id + 1 << " ... ";
        messager.flush("info");
        const int nev = particle_list_->read_in_particle_samples();
        messager << "nev = " << nev;
        messager.flush("info");
        if (nev == 0) continue;

        messager.info(" processing ...");
        bool have_selection = false;
        int previous_pid = 0;
        for (const auto &spectrum : spectra) {
            const int pid = spectrum->get_monval();
            if (!have_selection || previous_pid != pid) {
                particle_list_->filter_particles_from_events(pid);
                previous_pid = pid;
                have_selection = true;
            }
            spectrum->calculate_Qn_vector_shell(particle_list_);
        }
        messager.info("done!");
        event_id += nev;
    }
    for (const auto &spectrum : spectra) {
        spectrum->output_spectra_and_Qn_results();
    }
    particle_list_.reset();
}

void Analysis::FlowAnalysis() {
    InitializeAnalysis();
    FlowAnalysis(paraRdr_);
}

void Analysis::FlowAnalysis(const ParameterReader &config) {
    ParameterReader paraRdr = config;
    int compute_correlation = paraRdr.getVal("compute_correlation");
    if (compute_correlation == 1) paraRdr.setVal("compute_correlation", 0);
    int flag_charge_dependence = paraRdr.getVal("flag_charge_dependence");
    if (flag_charge_dependence == 1)
        paraRdr.setVal("flag_charge_dependence", 0);

    // first define all the analysis sets
    std::vector<std::unique_ptr<singleParticleSpectra>> spvn;
    const auto add_spectra = [&] {
        spvn.push_back(
            std::unique_ptr<singleParticleSpectra>(
                new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_)));
    };

    // all hadrons (99999) first (collect E_T)
    paraRdr.setVal("particle_monval", 99999);
    paraRdr.setVal("rap_type", 0);
    paraRdr.setVal("rapidity_distribution", 1);
    paraRdr.setVal("vn_rapidity_dis_pT_min", 0.);
    paraRdr.setVal("vn_rapidity_dis_pT_max", 4.);
    paraRdr.setVal("rapidityPTDistributionFlag", 1);
    add_spectra();

    // charged hadron first
    paraRdr.setVal("particle_monval", 9999);
    paraRdr.setVal("rap_type", 0);
    paraRdr.setVal("rapidity_distribution", 1);
    paraRdr.setVal("rapidityPTDistributionFlag", 0);
    paraRdr.setVal("rap_min", -0.5);
    paraRdr.setVal("rap_max", 0.5);
    paraRdr.setVal("vn_rapidity_dis_pT_min", 0.15);
    paraRdr.setVal("vn_rapidity_dis_pT_max", 2.0);
    add_spectra();
    paraRdr.setVal("vn_rapidity_dis_pT_min", 0.2);
    paraRdr.setVal("vn_rapidity_dis_pT_max", 3.0);
    add_spectra();
    // STAR rapidity cut
    paraRdr.setVal("rap_min", -1.0);
    paraRdr.setVal("rap_max", -0.5);
    add_spectra();
    paraRdr.setVal("rap_min", 0.5);
    paraRdr.setVal("rap_max", 1.0);
    add_spectra();
    paraRdr.setVal("rap_min", -1.0);
    paraRdr.setVal("rap_max", 1.0);
    add_spectra();
    // ALICE rapidity cut
    paraRdr.setVal("rap_min", -0.8);
    paraRdr.setVal("rap_max", 0.8);
    add_spectra();
    // STAR flow rapidity decorrelation cut
    paraRdr.setVal("vn_rapidity_dis_pT_min", 0.4);
    paraRdr.setVal("vn_rapidity_dis_pT_max", 4.0);
    paraRdr.setVal("rap_min", -1.0);
    paraRdr.setVal("rap_max", 1.0);
    add_spectra();
    // CMS default cut
    paraRdr.setVal("vn_rapidity_dis_pT_min", 0.3);
    paraRdr.setVal("vn_rapidity_dis_pT_max", 3.0);
    paraRdr.setVal("rap_min", -2.4);
    paraRdr.setVal("rap_max", 2.4);
    add_spectra();
    // ATLAS UPC cut
    paraRdr.setVal("vn_rapidity_dis_pT_min", 0.4);
    paraRdr.setVal("vn_rapidity_dis_pT_max", 2.0);
    paraRdr.setVal("rap_min", -2.5);
    paraRdr.setVal("rap_max", 2.5);
    add_spectra();
    // ATLAS default cut
    paraRdr.setVal("vn_rapidity_dis_pT_min", 0.5);
    paraRdr.setVal("vn_rapidity_dis_pT_max", 5.0);
    paraRdr.setVal("rap_min", -2.5);
    paraRdr.setVal("rap_max", 2.5);
    add_spectra();
    paraRdr.setVal("rap_min", 0.5);
    paraRdr.setVal("rap_max", 2.5);
    add_spectra();
    paraRdr.setVal("rap_min", -2.5);
    paraRdr.setVal("rap_max", -0.5);
    add_spectra();
    // ATLAS forward
    paraRdr.setVal("rap_min", -4.9);
    paraRdr.setVal("rap_max", -3.1);
    add_spectra();
    paraRdr.setVal("rap_min", 3.1);
    paraRdr.setVal("rap_max", 4.9);
    add_spectra();
    // ALICE V0A
    paraRdr.setVal("vn_rapidity_dis_pT_min", 0.2);
    paraRdr.setVal("vn_rapidity_dis_pT_max", 3.0);
    paraRdr.setVal("rap_min", -5.1);
    paraRdr.setVal("rap_max", -2.8);
    add_spectra();
    paraRdr.setVal("rap_min", 2.8);
    paraRdr.setVal("rap_max", 5.1);
    add_spectra();
    // ALICE V0C
    paraRdr.setVal("rap_min", 1.7);
    paraRdr.setVal("rap_max", 3.7);
    add_spectra();
    paraRdr.setVal("rap_min", -3.7);
    paraRdr.setVal("rap_max", -1.7);
    add_spectra();
    // PHENIX BBC
    paraRdr.setVal("vn_rapidity_dis_pT_min", 0.2);
    paraRdr.setVal("vn_rapidity_dis_pT_max", 2.0);
    paraRdr.setVal("rap_min", -3.9);
    paraRdr.setVal("rap_max", -3.1);
    add_spectra();
    paraRdr.setVal("rap_min", 3.1);
    paraRdr.setVal("rap_max", 3.9);
    add_spectra();

    // now identified particle
    paraRdr.setVal("rap_type", 1);
    paraRdr.setVal("rapidity_distribution", 1);
    paraRdr.setVal("vn_rapidity_dis_pT_min", 0.2);
    paraRdr.setVal("vn_rapidity_dis_pT_max", 3.0);
    paraRdr.setVal("rap_min", -0.5);
    paraRdr.setVal("rap_max", 0.5);
    paraRdr.setVal("particle_monval", 211);
    add_spectra();
    paraRdr.setVal("particle_monval", -211);
    add_spectra();
    paraRdr.setVal("particle_monval", 321);
    add_spectra();
    paraRdr.setVal("particle_monval", -321);
    add_spectra();
    paraRdr.setVal("particle_monval", 2212);
    add_spectra();
    paraRdr.setVal("particle_monval", -2212);
    add_spectra();
    if (paraRdr.getVal("resonance_weak_feed_down_flag") == 0) {
        paraRdr.setVal("particle_monval", 3122);
        add_spectra();
        paraRdr.setVal("particle_monval", -3122);
        add_spectra();
        paraRdr.setVal("particle_monval", 3312);
        add_spectra();
        paraRdr.setVal("particle_monval", -3312);
        add_spectra();
        paraRdr.setVal("particle_monval", 3334);
        add_spectra();
        paraRdr.setVal("particle_monval", -3334);
        add_spectra();
        paraRdr.setVal("particle_monval", 333);
        add_spectra();
    }
    if (paraRdr.getVal("collect_neutral_particles", 0) == 1) {
        paraRdr.setVal("particle_monval", 111);
        add_spectra();
        paraRdr.setVal("particle_monval", 311);
        add_spectra();
        paraRdr.setVal("particle_monval", -311);
        add_spectra();
        paraRdr.setVal("particle_monval", 2112);
        add_spectra();
        paraRdr.setVal("particle_monval", -2112);
        add_spectra();
    }
    paraRdr.setVal("rap_type", 0);
    paraRdr.setVal("particle_monval", 2212);
    add_spectra();
    paraRdr.setVal("particle_monval", -2212);
    add_spectra();

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
        add_spectra();
        paraRdr.setVal("rap_min", -2.0);
        paraRdr.setVal("rap_max", 2.0);
        add_spectra();
    }

    run_flow_spectra(config, spvn);
}

void Analysis::FlowAnalysis_multistrange_particles() {
    InitializeAnalysis();
    FlowAnalysis_multistrange_particles(paraRdr_);
}

void Analysis::FlowAnalysis_multistrange_particles(
    const ParameterReader &config) {
    ParameterReader paraRdr = config;
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
    std::vector<std::unique_ptr<singleParticleSpectra>> spvn;
    const auto add_spectra = [&] {
        spvn.push_back(
            std::unique_ptr<singleParticleSpectra>(
                new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_)));
    };
    paraRdr.setVal("particle_monval", 3122);
    add_spectra();
    paraRdr.setVal("particle_monval", -3122);
    add_spectra();
    paraRdr.setVal("particle_monval", 3312);
    add_spectra();
    paraRdr.setVal("particle_monval", -3312);
    add_spectra();
    paraRdr.setVal("particle_monval", 3334);
    add_spectra();
    paraRdr.setVal("particle_monval", -3334);
    add_spectra();
    paraRdr.setVal("particle_monval", 333);
    add_spectra();

    run_flow_spectra(config, spvn);
}

void Analysis::FlowAnalysis_3D() {
    InitializeAnalysis();
    FlowAnalysis_3D(paraRdr_);
}

void Analysis::FlowAnalysis_3D(const ParameterReader &config) {
    ParameterReader paraRdr = config;
    paraRdr.setVal("compute_correlation", 0);
    paraRdr.setVal("flag_charge_dependence", 0);

    // first define all the analysis sets
    std::vector<std::unique_ptr<singleParticleSpectra>> spvn;
    const auto add_spectra = [&] {
        spvn.push_back(
            std::unique_ptr<singleParticleSpectra>(
                new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_)));
    };

    // all hadrons (99999) first (collect E_T)
    paraRdr.setVal("particle_monval", 99999);
    paraRdr.setVal("rap_type", 0);
    paraRdr.setVal("rapidity_distribution", 1);
    paraRdr.setVal("vn_rapidity_dis_pT_min", 0.);
    paraRdr.setVal("vn_rapidity_dis_pT_max", 4.);
    paraRdr.setVal("rapidityPTDistributionFlag", 1);
    add_spectra();

    // charged hadron first
    paraRdr.setVal("particle_monval", 9999);
    paraRdr.setVal("rap_type", 0);
    paraRdr.setVal("rapidity_distribution", 1);
    paraRdr.setVal("vn_rapidity_dis_pT_min", 0.2);
    paraRdr.setVal("vn_rapidity_dis_pT_max", 3.0);
    paraRdr.setVal("rapidityPTDistributionFlag", 0);
    paraRdr.setVal("rap_min", -2.5);
    paraRdr.setVal("rap_max", 2.5);
    add_spectra();
    paraRdr.setVal("rap_min", 0.5);
    paraRdr.setVal("rap_max", 2.5);
    add_spectra();
    paraRdr.setVal("rap_min", -2.5);
    paraRdr.setVal("rap_max", -0.5);
    add_spectra();
    paraRdr.setVal("rap_min", -0.5);
    paraRdr.setVal("rap_max", 0.5);
    paraRdr.setVal("rapidityPTDistributionFlag", 1);
    add_spectra();

    // now identified particle
    paraRdr.setVal("rap_type", 1);
    paraRdr.setVal("rapidityPTDistributionFlag", 0);
    paraRdr.setVal("particle_monval", 211);
    add_spectra();
    paraRdr.setVal("particle_monval", -211);
    add_spectra();
    paraRdr.setVal("particle_monval", 321);
    add_spectra();
    paraRdr.setVal("particle_monval", -321);
    add_spectra();
    paraRdr.setVal("particle_monval", 2212);
    add_spectra();
    paraRdr.setVal("particle_monval", -2212);
    add_spectra();
    if (paraRdr.getVal("resonance_weak_feed_down_flag") == 0) {
        paraRdr.setVal("particle_monval", 3122);
        add_spectra();
        paraRdr.setVal("particle_monval", -3122);
        add_spectra();
        paraRdr.setVal("particle_monval", 3312);
        add_spectra();
        paraRdr.setVal("particle_monval", -3312);
        add_spectra();
        paraRdr.setVal("particle_monval", 3334);
        add_spectra();
        paraRdr.setVal("particle_monval", -3334);
        add_spectra();
        paraRdr.setVal("particle_monval", 333);
        add_spectra();
    }

    run_flow_spectra(config, spvn);
}

void Analysis::FlowAnalysis_RHIC() {
    InitializeAnalysis();
    FlowAnalysis_RHIC(paraRdr_);
}

void Analysis::FlowAnalysis_RHIC(const ParameterReader &config) {
    ParameterReader paraRdr = config;
    int compute_correlation = paraRdr.getVal("compute_correlation");
    if (compute_correlation == 1) paraRdr.setVal("compute_correlation", 0);
    int flag_charge_dependence = paraRdr.getVal("flag_charge_dependence");
    if (flag_charge_dependence == 1)
        paraRdr.setVal("flag_charge_dependence", 0);

    // first define all the analysis sets
    std::vector<std::unique_ptr<singleParticleSpectra>> spvn;
    const auto add_spectra = [&] {
        spvn.push_back(
            std::unique_ptr<singleParticleSpectra>(
                new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_)));
    };

    // all hadrons (99999) first (collect E_T)
    paraRdr.setVal("particle_monval", 99999);
    paraRdr.setVal("rap_type", 0);
    paraRdr.setVal("rapidity_distribution", 1);
    paraRdr.setVal("vn_rapidity_dis_pT_min", 0.);
    paraRdr.setVal("vn_rapidity_dis_pT_max", 4.);
    paraRdr.setVal("rapidityPTDistributionFlag", 1);
    add_spectra();

    // charged hadron first
    paraRdr.setVal("particle_monval", 9999);
    paraRdr.setVal("rap_type", 0);
    paraRdr.setVal("rapidity_distribution", 1);
    paraRdr.setVal("vn_rapidity_dis_pT_min", 0.2);
    paraRdr.setVal("vn_rapidity_dis_pT_max", 3.0);
    paraRdr.setVal("rapidityPTDistributionFlag", 1);
    add_spectra();

    // speical triggers follows
    paraRdr.setVal("rapidityPTDistributionFlag", 0);
    paraRdr.setVal("rap_min", -0.5);
    paraRdr.setVal("rap_max", 0.5);
    add_spectra();
    // STAR rapidity cut
    paraRdr.setVal("rap_min", -1.0);
    paraRdr.setVal("rap_max", -0.5);
    add_spectra();
    paraRdr.setVal("rap_min", 0.5);
    paraRdr.setVal("rap_max", 1.0);
    add_spectra();
    // PHENIX BBC
    paraRdr.setVal("rap_min", -3.9);
    paraRdr.setVal("rap_max", -3.1);
    add_spectra();
    paraRdr.setVal("rap_min", 3.1);
    paraRdr.setVal("rap_max", 3.9);
    add_spectra();

    // now identified particle
    paraRdr.setVal("vn_rapidity_dis_pT_min", 0.2);
    paraRdr.setVal("vn_rapidity_dis_pT_max", 3.0);
    if (paraRdr.getVal("pidwithRapidityPTDistribution", 0) == 1) {
        paraRdr.setVal("rapidityPTDistributionFlag", 1);
        paraRdr.setVal("rap_type", 1);
        paraRdr.setVal("particle_monval", 211);
        add_spectra();
        paraRdr.setVal("particle_monval", -211);
        add_spectra();
        paraRdr.setVal("particle_monval", 321);
        add_spectra();
        paraRdr.setVal("particle_monval", -321);
        add_spectra();
        paraRdr.setVal("particle_monval", 2212);
        add_spectra();
        paraRdr.setVal("particle_monval", -2212);
        add_spectra();
    }
    paraRdr.setVal("rap_type", 1);
    paraRdr.setVal("rapidity_distribution", 1);
    paraRdr.setVal("rap_min", -0.5);
    paraRdr.setVal("rap_max", 0.5);
    paraRdr.setVal("particle_monval", 211);
    add_spectra();
    paraRdr.setVal("particle_monval", -211);
    add_spectra();
    paraRdr.setVal("particle_monval", 321);
    add_spectra();
    paraRdr.setVal("particle_monval", -321);
    add_spectra();
    paraRdr.setVal("particle_monval", 2212);
    add_spectra();
    paraRdr.setVal("particle_monval", -2212);
    add_spectra();
    if (paraRdr.getVal("resonance_weak_feed_down_flag") == 0) {
        paraRdr.setVal("particle_monval", 3122);
        add_spectra();
        paraRdr.setVal("particle_monval", -3122);
        add_spectra();
        paraRdr.setVal("particle_monval", 3312);
        add_spectra();
        paraRdr.setVal("particle_monval", -3312);
        add_spectra();
        paraRdr.setVal("particle_monval", 3334);
        add_spectra();
        paraRdr.setVal("particle_monval", -3334);
        add_spectra();
        paraRdr.setVal("particle_monval", 333);
        add_spectra();
    }
    if (paraRdr.getVal("collect_neutral_particles", 0) == 1) {
        paraRdr.setVal("particle_monval", 111);
        add_spectra();
        paraRdr.setVal("particle_monval", 311);
        add_spectra();
        paraRdr.setVal("particle_monval", -311);
        add_spectra();
        paraRdr.setVal("particle_monval", 2112);
        add_spectra();
        paraRdr.setVal("particle_monval", -2112);
        add_spectra();
    }

    run_flow_spectra(config, spvn);
}

void Analysis::FlowAnalysis_LHC() {
    InitializeAnalysis();
    FlowAnalysis_LHC(paraRdr_);
}

void Analysis::FlowAnalysis_LHC(const ParameterReader &config) {
    ParameterReader paraRdr = config;
    int compute_correlation = paraRdr.getVal("compute_correlation");
    if (compute_correlation == 1) paraRdr.setVal("compute_correlation", 0);
    int flag_charge_dependence = paraRdr.getVal("flag_charge_dependence");
    if (flag_charge_dependence == 1)
        paraRdr.setVal("flag_charge_dependence", 0);

    // first define all the analysis sets
    std::vector<std::unique_ptr<singleParticleSpectra>> spvn;
    const auto add_spectra = [&] {
        spvn.push_back(
            std::unique_ptr<singleParticleSpectra>(
                new singleParticleSpectra(paraRdr, path_, ran_gen_ptr_)));
    };

    // all hadrons (99999) first (collect E_T)
    paraRdr.setVal("particle_monval", 99999);
    paraRdr.setVal("rap_type", 0);
    paraRdr.setVal("rapidity_distribution", 1);
    paraRdr.setVal("vn_rapidity_dis_pT_min", 0.);
    paraRdr.setVal("vn_rapidity_dis_pT_max", 4.);
    paraRdr.setVal("rapidityPTDistributionFlag", 1);
    add_spectra();

    // charged hadron first
    paraRdr.setVal("particle_monval", 9999);
    paraRdr.setVal("rap_type", 0);
    paraRdr.setVal("rapidity_distribution", 1);
    paraRdr.setVal("rapidityPTDistributionFlag", 1);
    paraRdr.setVal("vn_rapidity_dis_pT_min", 0.2);
    paraRdr.setVal("vn_rapidity_dis_pT_max", 3.0);
    add_spectra();

    // special trigger follows
    paraRdr.setVal("rapidityPTDistributionFlag", 0);
    paraRdr.setVal("rap_min", -0.5);
    paraRdr.setVal("rap_max", 0.5);
    paraRdr.setVal("vn_rapidity_dis_pT_min", 0.2);
    paraRdr.setVal("vn_rapidity_dis_pT_max", 3.0);
    add_spectra();
    // ALICE rapidity cut
    paraRdr.setVal("rap_min", -0.8);
    paraRdr.setVal("rap_max", -0.4);
    add_spectra();
    paraRdr.setVal("rap_min", 0.4);
    paraRdr.setVal("rap_max", 0.8);
    add_spectra();
    paraRdr.setVal("rap_min", -0.4);
    paraRdr.setVal("rap_max", 0.4);
    add_spectra();
    // ALICE V0A
    paraRdr.setVal("rap_min", -5.1);
    paraRdr.setVal("rap_max", -2.8);
    add_spectra();
    paraRdr.setVal("rap_min", 2.8);
    paraRdr.setVal("rap_max", 5.1);
    add_spectra();
    // ALICE V0C
    paraRdr.setVal("rap_min", 1.7);
    paraRdr.setVal("rap_max", 3.7);
    add_spectra();
    paraRdr.setVal("rap_min", -3.7);
    paraRdr.setVal("rap_max", -1.7);
    add_spectra();
    // CMS default cut
    paraRdr.setVal("vn_rapidity_dis_pT_min", 0.3);
    paraRdr.setVal("vn_rapidity_dis_pT_max", 3.0);
    paraRdr.setVal("rap_min", -2.4);
    paraRdr.setVal("rap_max", -0.5);
    add_spectra();
    paraRdr.setVal("rap_min", 0.5);
    paraRdr.setVal("rap_max", 2.4);
    add_spectra();
    // ATLAS default cut
    paraRdr.setVal("vn_rapidity_dis_pT_min", 0.5);
    paraRdr.setVal("vn_rapidity_dis_pT_max", 5.0);
    paraRdr.setVal("rap_min", 0.5);
    paraRdr.setVal("rap_max", 2.5);
    add_spectra();
    paraRdr.setVal("rap_min", -2.5);
    paraRdr.setVal("rap_max", -0.5);
    add_spectra();
    // ATLAS forward
    paraRdr.setVal("rap_min", -4.9);
    paraRdr.setVal("rap_max", -3.1);
    add_spectra();
    paraRdr.setVal("rap_min", 3.1);
    paraRdr.setVal("rap_max", 4.9);
    add_spectra();

    // now identified particle
    paraRdr.setVal("rap_type", 1);
    paraRdr.setVal("rapidity_distribution", 1);
    paraRdr.setVal("rapidityPTDistributionFlag", 0);
    paraRdr.setVal("rap_min", -0.5);
    paraRdr.setVal("rap_max", 0.5);
    paraRdr.setVal("vn_rapidity_dis_pT_min", 0.0);
    paraRdr.setVal("vn_rapidity_dis_pT_max", 4.0);
    paraRdr.setVal("particle_monval", 211);
    add_spectra();
    paraRdr.setVal("particle_monval", -211);
    add_spectra();
    paraRdr.setVal("particle_monval", 321);
    add_spectra();
    paraRdr.setVal("particle_monval", -321);
    add_spectra();
    paraRdr.setVal("particle_monval", 2212);
    add_spectra();
    paraRdr.setVal("particle_monval", -2212);
    add_spectra();
    if (paraRdr.getVal("resonance_weak_feed_down_flag") == 0) {
        paraRdr.setVal("particle_monval", 3122);
        add_spectra();
        paraRdr.setVal("particle_monval", -3122);
        add_spectra();
        paraRdr.setVal("particle_monval", 3312);
        add_spectra();
        paraRdr.setVal("particle_monval", -3312);
        add_spectra();
        paraRdr.setVal("particle_monval", 3334);
        add_spectra();
        paraRdr.setVal("particle_monval", -3334);
        add_spectra();
        paraRdr.setVal("particle_monval", 333);
        add_spectra();
    }
    if (paraRdr.getVal("collect_neutral_particles", 0) == 1) {
        paraRdr.setVal("particle_monval", 111);
        add_spectra();
        paraRdr.setVal("particle_monval", 311);
        add_spectra();
        paraRdr.setVal("particle_monval", -311);
        add_spectra();
        paraRdr.setVal("particle_monval", 2112);
        add_spectra();
        paraRdr.setVal("particle_monval", -2112);
        add_spectra();
    }
    if (paraRdr.getVal("pidwithPseudoRapCuts", 0) == 1) {
        paraRdr.setVal("rap_type", 0);
        paraRdr.setVal("rapidity_distribution", 1);
        paraRdr.setVal("rap_min", -0.8);
        paraRdr.setVal("rap_max", 0.8);
        paraRdr.setVal("particle_monval", 211);
        add_spectra();
        paraRdr.setVal("particle_monval", -211);
        add_spectra();
        paraRdr.setVal("particle_monval", 321);
        add_spectra();
        paraRdr.setVal("particle_monval", -321);
        add_spectra();
        paraRdr.setVal("particle_monval", 2212);
        add_spectra();
        paraRdr.setVal("particle_monval", -2212);
        add_spectra();
    }

    paraRdr.setVal("rap_type", 1);
    paraRdr.setVal("rap_min", -0.5);
    paraRdr.setVal("rap_max", 0.5);
    paraRdr.setVal("rapidity_distribution", 1);
    paraRdr.setVal("vn_rapidity_dis_pT_min", 0.2);
    paraRdr.setVal("vn_rapidity_dis_pT_max", 3.0);
    paraRdr.setVal("particle_monval", 211);
    add_spectra();
    paraRdr.setVal("particle_monval", -211);
    add_spectra();
    paraRdr.setVal("particle_monval", 321);
    add_spectra();
    paraRdr.setVal("particle_monval", -321);
    add_spectra();
    paraRdr.setVal("particle_monval", 2212);
    add_spectra();
    paraRdr.setVal("particle_monval", -2212);
    add_spectra();
    if (paraRdr.getVal("pidwithRapidityPTDistribution", 0) == 1) {
        paraRdr.setVal("rapidityPTDistributionFlag", 1);
        paraRdr.setVal("rap_type", 1);
        paraRdr.setVal("particle_monval", 211);
        add_spectra();
        paraRdr.setVal("particle_monval", -211);
        add_spectra();
        paraRdr.setVal("particle_monval", 321);
        add_spectra();
        paraRdr.setVal("particle_monval", -321);
        add_spectra();
        paraRdr.setVal("particle_monval", 2212);
        add_spectra();
        paraRdr.setVal("particle_monval", -2212);
        add_spectra();
    }

    // lastly, if we want to compute multi-particle correlations within
    // the same UrQMD events
    if (compute_correlation == 1) {
        paraRdr.setVal("vn_rapidity_dis_pT_min", 0.2);
        paraRdr.setVal("vn_rapidity_dis_pT_max", 3.0);
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
        add_spectra();
        paraRdr.setVal("rap_min", -2.0);
        paraRdr.setVal("rap_max", 2.0);
        add_spectra();
    }

    run_flow_spectra(config, spvn);
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
