// Copyright Chun Shen @ 2020

#include <string>
#include <memory>
#include "Analysis.h"


Analysis::Analysis(std::string path) : path_(path) {

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
