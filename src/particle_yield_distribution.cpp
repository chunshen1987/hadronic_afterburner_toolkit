// Copyright Chun Shen @ 2016
#include "./particle_yield_distribution.h"

#include <stdlib.h>

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

#include "./parameters.h"

using namespace std;

particle_yield_distribution::particle_yield_distribution(
    ParameterReader &paraRdr, std::string path)
    : paraRdr_(paraRdr), path_(path) {
    particle_monval = paraRdr_.getVal("particle_monval");
    net_particle_flag = paraRdr_.getVal("net_particle_flag");
    n_max = 6000;
    number_of_events = new int[n_max];
    for (int i = 0; i < n_max; i++) number_of_events[i] = 0;

    pT_min = paraRdr_.getVal("pT_min");
    pT_max = paraRdr_.getVal("pT_max");

    total_number_of_events = 0;

    rap_type = paraRdr_.getVal("rap_type");
    rap_min = paraRdr_.getVal("rap_min");
    rap_max = paraRdr_.getVal("rap_max");

    if (particle_monval == 9999)  // use pseudo-rapidity for all charged hadrons
        rap_type = 0;
}

particle_yield_distribution::~particle_yield_distribution() {
    delete[] number_of_events;
}

void particle_yield_distribution::collect_particle_yield_distribution(
    std::shared_ptr<particleSamples> particle_list_in) {
    set_particle_list(particle_list_in);
    int nev = particle_list->get_number_of_events();
    for (int iev = 0; iev < nev; iev++) {
        collect_particle_yield(iev);
    }
    total_number_of_events += nev;
}

void particle_yield_distribution::collect_particle_yield(int event_id) {
    const auto passes_cuts = [this](const particle_info &particle) {
        const double rapidity =
            rap_type == 0 ? particle.rap_eta : particle.rap_y;
        const double pT = std::hypot(particle.px, particle.py);
        return rapidity > rap_min && rapidity < rap_max && pT > pT_min
               && pT < pT_max;
    };
    int count = 0;
    for (int i = 0; i < particle_list->get_number_of_particles(event_id); i++) {
        count += passes_cuts(particle_list->get_particle(event_id, i));
    }
    int count_anti_particle = 0;
    if (net_particle_flag == 1) {
        const int anti_count =
            particle_list->get_number_of_anti_particles(event_id);
        for (int i = 0; i < anti_count; i++) {
            count_anti_particle +=
                passes_cuts(particle_list->get_anti_particle(event_id, i));
        }
    }
    int idx = count - count_anti_particle + n_max / 2;
    if (idx < 0 || idx > n_max - 1) {
        cout << "particle number is out of bound!" << endl;
        cout << "bound: [" << -n_max / 2 << " , " << n_max / 2 << "]" << endl;
        cout << "count = " << count
             << ", count_anti_particle = " << count_anti_particle << endl;
        exit(1);
    }
    number_of_events[idx]++;
}

void particle_yield_distribution::output_particle_yield_distribution() {
    // this function outputs the particle yield distribution

    ostringstream filename;
    if (net_particle_flag == 0) {
        if (rap_type == 0)
            filename << path_ << "/particle_" << particle_monval
                     << "_yield_distribution_eta.dat";
        else
            filename << path_ << "/particle_" << particle_monval
                     << "_yield_distribution_y.dat";
    } else {
        if (rap_type == 0)
            filename << path_ << "/particle_0" << particle_monval
                     << "_yield_distribution_eta.dat";
        else
            filename << path_ << "/particle_0" << particle_monval
                     << "_yield_distribution_y.dat";
    }
    ofstream output(filename.str().c_str());

    // output header
    output << "# N  P(N)" << endl;

    for (int i = 0; i < n_max; i++) {
        const double p_N = total_number_of_events > 0
                               ? static_cast<double>(number_of_events[i])
                                     / total_number_of_events
                               : 0.0;
        output << scientific << setw(18) << setprecision(8) << i - n_max / 2
               << "   " << p_N << endl;
    }
    output.close();
}
