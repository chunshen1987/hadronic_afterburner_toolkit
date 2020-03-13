// Copyright Chun Shen @ 2020

#include <string>
#include <vector>
#include <memory>
#include "Analysis.h"
#include "single_particleSpectra.h"

using std::vector;

Analysis::Analysis(std::string path) : path_(path) {

}


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


void Analysis::FlowAnalysis() {
    // start the loop
    int event_id = 0;
    singleParticleSpectra pion_p(paraRdr_, path_, ran_gen_ptr_);
    paraRdr_.setVal("particle_monval", 321);
    singleParticleSpectra kaon_p(paraRdr_, path_, ran_gen_ptr_);
    while (!particle_list_->end_of_file()) {
        messager << "Reading event: " << event_id + 1 << " ... ";
        messager.flush("info");
        int nev = particle_list_->read_in_particle_samples();
        messager << "nev = " << nev;
        messager.flush("info");
        messager.info(" processing ...");
        int particle_monval = 211;
        particle_list_->filter_particles_from_events(particle_monval);
        pion_p.calculate_Qn_vector_shell(particle_list_);
        particle_monval = 321;
        particle_list_->filter_particles_from_events(particle_monval);
        kaon_p.calculate_Qn_vector_shell(particle_list_);
    }
    pion_p.output_spectra_and_Qn_results();
    kaon_p.output_spectra_and_Qn_results();
}
