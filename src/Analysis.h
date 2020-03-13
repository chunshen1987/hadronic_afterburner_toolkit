#ifndef SRC_analysis_h_
#define SRC_analysis_h_

#include <string>
#include <memory>
#include "ParameterReader.h"
#include "particleSamples.h"
#include "Random.h"
#include "pretty_ostream.h"

class Analysis {
 private:
    const std::string path_;
    ParameterReader paraRdr_;
    std::shared_ptr<RandomUtil::Random> ran_gen_ptr_;
    std::shared_ptr<particleSamples> particle_list_;
    pretty_ostream messager;

    int run_mode_;

 public:
    Analysis(std::string path);
    ~Analysis() {};

    void UpdateParameterDict(std::string param_filename);
    void UpdateParameterDict(std::string param_filename,
                             int argc, char *argv[]);

    void InitializeAnalysis();
    void FlowAnalysis();
};

#endif  // SRC_analysis_h_
