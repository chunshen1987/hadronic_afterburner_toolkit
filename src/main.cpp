//=============================================================================
//  calculate the HBT correlation functions from MC samples of particles
//
//  Programmer: Chun Shen
//       Email: chunshen@physics.mcgill.ca
//
//        Date: 10/31/14
//
//  To do:
//
//=============================================================================

#include<iostream>
#include<string>

#include "Stopwatch.h"
#include "Analysis.h"

using namespace std;

int main(int argc, char *argv[]) {
    // program title
    cout << "*********************************************************";
    cout << endl
         << "               iHBT_afterbuner         " << endl
         << endl
         << "  Ver 1.0   ----- Chun Shen, 10/2014   " << endl;
    cout << endl << "*********************************************************"
         << endl;

    const std::string path="results";

    Stopwatch sw;
    Stopwatch sw_total;
    sw_total.tic();
    sw.tic();

    Analysis flow_analysis(path);
    flow_analysis.UpdateParameterDict("parameters.dat", argc, argv);
    flow_analysis.PerformAnalysis();

    //// Read-in parameters
    //ParameterReader paraRdr;
    //paraRdr.readFromFile("parameters.dat");
    //paraRdr.readFromArguments(argc, argv);
    //paraRdr.echo();
    //int randomSeed = paraRdr.getVal("randomSeed");
    //std::shared_ptr<RandomUtil::Random> ran_gen_ptr(
    //                                new RandomUtil::Random(randomSeed));

    //auto particle_list = std::make_shared<particleSamples> (paraRdr, path,
    //                                                        ran_gen_ptr);
    //if (run_mode == 0) {
    //    // collect single particle spectra and vn
    //    singleParticleSpectra testSP(paraRdr, path, ran_gen_ptr);
    //    int event_id = 0;
    //    while (!particle_list->end_of_file()) {
    //        std::cout << "Reading event: " << event_id + 1 << " ... "
    //                  << std::endl;
    //        int nev = particle_list->read_in_particle_samples();
    //        std::cout << "nev = " << nev << std::endl;
    //        std::cout << " processing ..." << std::endl;
    //        int monval = paraRdr.getVal("particle_monval");
    //        particle_list->filter_particles_from_events(monval);
    //        testSP.calculate_Qn_vector_shell(particle_list);
    //    }
    //    testSP.output_spectra_and_Qn_results();
    //} else if (run_mode == 1) {
    //    // compute HBT correlation function and HBT radii
    //    HBT_correlation test(paraRdr, path, ran_gen_ptr, particle_list);
    //    test.calculate_HBT_correlation_function();
    //} else if (run_mode == 2) {
    //    // collect event-by-event particle yield distribution
    //    particle_yield_distribution test_dis(paraRdr, path, particle_list);
    //    test_dis.collect_particle_yield_distribution();
    //} else if (run_mode == 3) {
    //    // compute balance function
    //    BalanceFunction test(paraRdr, path, ran_gen_ptr, particle_list);
    //    test.calculate_balance_function();
    //} else {
    //    cout << "Error: unrecognized run_mode: " << run_mode << endl;
    //    exit(1);
    //}

    sw_total.toc();
    cout << "Program totally finished in " << sw_total.takeTime() << " sec."
         << endl;
    return 0;
}
