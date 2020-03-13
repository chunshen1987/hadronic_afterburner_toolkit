// Copyright Chun Shen @ 2020

#include <string>
#include <vector>
#include <memory>
#include <map>
#include "Analysis.h"
#include "single_particleSpectra.h"

using std::vector;

Analysis::Analysis(std::string path) : path_(path) {}


void Analysis::UpdateParameterDict(std::string param_filename) {
    paraRdr_.readFromFile(param_filename);
    paraRdr_.echo();
}


void Analysis::UpdateParameterDict(std::string param_filename,
                                   int argc, char *argv[]) {
    paraRdr_.readFromFile(param_filename);
    paraRdr_.readFromArguments(argc, argv);
    paraRdr_.echo();
}


void Analysis::InitializeAnalysis() {
    run_mode_      = paraRdr_.getVal("run_mode");
    int randomSeed = paraRdr_.getVal("randomSeed");
    ran_gen_ptr_   = std::make_shared<RandomUtil::Random> (randomSeed);
    particle_list_ = std::make_shared<particleSamples> (paraRdr_, path_,
                                                        ran_gen_ptr_);
}


void Analysis::PerformAnalysis() {
    InitializeAnalysis();
    if (run_mode_ == 0) {
        // collect single particle spectra and vn
        FlowAnalysis();
    }
}


void Analysis::FlowAnalysis() {
    // first define all the analysis sets
    std::vector<singleParticleSpectra*> spvn;
    // charged hadron first
    paraRdr_.setVal("particle_monval", 9999);
    paraRdr_.setVal("rap_type", 0);
    paraRdr_.setVal("rap_min", -0.5); paraRdr_.setVal("rap_max", 0.5);
    paraRdr_.setVal("vn_rapidity_dis_pT_min", 0.15);
    paraRdr_.setVal("vn_rapidity_dis_pT_max", 2.0);
    spvn.push_back(new singleParticleSpectra(paraRdr_, path_, ran_gen_ptr_));
    paraRdr_.setVal("vn_rapidity_dis_pT_min", 0.2);
    paraRdr_.setVal("vn_rapidity_dis_pT_max", 3.0);
    spvn.push_back(new singleParticleSpectra(paraRdr_, path_, ran_gen_ptr_));
    paraRdr_.setVal("rap_min", -1.0); paraRdr_.setVal("rap_max", -0.1);
    spvn.push_back(new singleParticleSpectra(paraRdr_, path_, ran_gen_ptr_));
    paraRdr_.setVal("rap_min", 0.1); paraRdr_.setVal("rap_max", 1.0);
    spvn.push_back(new singleParticleSpectra(paraRdr_, path_, ran_gen_ptr_));
    paraRdr_.setVal("rap_min", -1.0); paraRdr_.setVal("rap_max", 1.0);
    spvn.push_back(new singleParticleSpectra(paraRdr_, path_, ran_gen_ptr_));
    paraRdr_.setVal("rap_min", -0.8); paraRdr_.setVal("rap_max", 0.8);
    spvn.push_back(new singleParticleSpectra(paraRdr_, path_, ran_gen_ptr_));
    paraRdr_.setVal("rap_min", -2.4); paraRdr_.setVal("rap_max", 2.4);
    spvn.push_back(new singleParticleSpectra(paraRdr_, path_, ran_gen_ptr_));
    paraRdr_.setVal("rap_min", 0.5); paraRdr_.setVal("rap_max", 2.0);
    spvn.push_back(new singleParticleSpectra(paraRdr_, path_, ran_gen_ptr_));
    paraRdr_.setVal("rap_min", -2.0); paraRdr_.setVal("rap_max", -0.5);
    spvn.push_back(new singleParticleSpectra(paraRdr_, path_, ran_gen_ptr_));

    // now identified particle
    paraRdr_.setVal("rap_type", 1);
    paraRdr_.setVal("rap_min", -0.5); paraRdr_.setVal("rap_max", 0.5);
    paraRdr_.setVal("particle_monval", 211);
    spvn.push_back(new singleParticleSpectra(paraRdr_, path_, ran_gen_ptr_));
    paraRdr_.setVal("particle_monval", -211);
    spvn.push_back(new singleParticleSpectra(paraRdr_, path_, ran_gen_ptr_));
    paraRdr_.setVal("particle_monval", 321);
    spvn.push_back(new singleParticleSpectra(paraRdr_, path_, ran_gen_ptr_));
    paraRdr_.setVal("particle_monval", -321);
    spvn.push_back(new singleParticleSpectra(paraRdr_, path_, ran_gen_ptr_));
    paraRdr_.setVal("particle_monval", 2212);
    spvn.push_back(new singleParticleSpectra(paraRdr_, path_, ran_gen_ptr_));
    paraRdr_.setVal("particle_monval", -2212);
    spvn.push_back(new singleParticleSpectra(paraRdr_, path_, ran_gen_ptr_));
    paraRdr_.setVal("resonance_weak_feed_down_flag", 0);
    paraRdr_.setVal("particle_monval", 3122);
    spvn.push_back(new singleParticleSpectra(paraRdr_, path_, ran_gen_ptr_));
    paraRdr_.setVal("particle_monval", -3122);
    spvn.push_back(new singleParticleSpectra(paraRdr_, path_, ran_gen_ptr_));
    paraRdr_.setVal("particle_monval", 3312);
    spvn.push_back(new singleParticleSpectra(paraRdr_, path_, ran_gen_ptr_));
    paraRdr_.setVal("particle_monval", -3312);
    spvn.push_back(new singleParticleSpectra(paraRdr_, path_, ran_gen_ptr_));
    paraRdr_.setVal("particle_monval", 3334);
    spvn.push_back(new singleParticleSpectra(paraRdr_, path_, ran_gen_ptr_));
    paraRdr_.setVal("particle_monval", -3334);
    spvn.push_back(new singleParticleSpectra(paraRdr_, path_, ran_gen_ptr_));
    paraRdr_.setVal("particle_monval", 333);
    spvn.push_back(new singleParticleSpectra(paraRdr_, path_, ran_gen_ptr_));

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
    }
    for (auto &ipart : spvn)
        ipart->output_spectra_and_Qn_results();
}
