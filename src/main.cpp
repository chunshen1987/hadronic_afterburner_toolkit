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

    sw_total.toc();
    cout << "Program totally finished in " << sw_total.takeTime() << " sec."
         << endl;
    return 0;
}
