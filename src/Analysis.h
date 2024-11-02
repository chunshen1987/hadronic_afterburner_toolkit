#ifndef SRC_analysis_h_
#define SRC_analysis_h_

#include <memory>
#include <string>

#include "ParameterReader.h"
#include "Random.h"
#include "particleSamples.h"
#include "pretty_ostream.h"

class Analysis {
  private:
    const std::string path_;
    ParameterReader paraRdr_;
    std::shared_ptr<RandomUtil::Random> ran_gen_ptr_;
    std::shared_ptr<particleSamples> particle_list_;
    pretty_ostream messager;

  public:
    Analysis(std::string path);
    ~Analysis() {};

    void UpdateParameterDict(std::string param_filename);
    void UpdateParameterDict(
        std::string param_filename, int argc, char *argv[]);
    void setRapiditySfhitFromFile();

    void InitializeAnalysis();
    void PerformAnalysis();
    void FlowAnalysis();
    void FlowAnalysis_3D();
    void FlowAnalysis_RHIC();
    void FlowAnalysis_LHC();
    void FlowAnalysis_multistrange_particles();
    void HBTAnalysis();
    void ParticleYieldDistributionAnalysis();
    void BalanceFunctionAnalysis();
};

#endif  // SRC_analysis_h_
