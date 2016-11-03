// Copyright Chun Shen @ 2016
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <iomanip>
#include "./parameters.h"
#include "./particle_yield_distribution.h"

using namespace std;

particle_yield_distribution::particle_yield_distribution(
                                ParameterReader *paraRdr_in, string path_in, 
                                particleSamples *particle_list_in) {
    paraRdr = paraRdr_in;
    path = path_in;
    particle_list = particle_list_in;

    particle_monval = paraRdr->getVal("particle_monval");

    if (particle_monval == 333) {
        // phi(1020) is reconstructed from (K^+, K^-) pairs
        reconst_branching_ratio = 0.489;
    } else {
        reconst_branching_ratio = 1.0;
    }

    n_max = 2000;
    number_of_events = new int[n_max];
    for (int i = 0; i < n_max; i++)
        number_of_events[i] = 0;

    pT_min = paraRdr->getVal("pT_min");
    pT_max = paraRdr->getVal("pT_max");
    
    total_number_of_events = 0;

    rap_type = paraRdr->getVal("rap_type");
    rap_min = paraRdr->getVal("rap_min");
    rap_max = paraRdr->getVal("rap_max");

    if (particle_monval == 9999)  // use pseudo-rapidity for all charged hadrons
        rap_type = 0;
}

particle_yield_distribution::~particle_yield_distribution() {
    delete [] number_of_events;
}

void particle_yield_distribution::collect_particle_yield_distribution() {
    int event_id = 0;
    int buffer_size = particle_list->get_event_buffer_size();
    while (!particle_list->end_of_file()) {
        cout << "Reading event: " << event_id+1 
             << "-" << event_id + buffer_size << " ... " << flush;
        particle_list->read_in_particle_samples();
        cout << " processing ..." << flush;
        int nev = particle_list->get_number_of_events();
        for (int iev = 0; iev < nev; iev++) {
            event_id++;
            collect_particle_yield(iev);
        }
        cout << " done!" << endl;
    }
    total_number_of_events = event_id;
    output_particle_yield_distribution();
}

void particle_yield_distribution::collect_particle_yield(int event_id) {
    int number_of_particles = particle_list->get_number_of_particles(event_id);

    int count = 0;
    for (int i = 0; i < number_of_particles; i++) {
        double pz_local = particle_list->get_particle(event_id, i).pz;
        double E_local = particle_list->get_particle(event_id, i).E;

        double rap_local;
        if (rap_type == 0) {
            double mass = particle_list->get_particle(event_id, i).mass;
            double pmag = sqrt(E_local*E_local - mass*mass);
            rap_local = 0.5*log((pmag + pz_local)/(pmag - pz_local));
        } else {
            rap_local = 0.5*log((E_local + pz_local)/(E_local - pz_local));
        }

        if (rap_local > rap_min && rap_local < rap_max) {
            double px_local = particle_list->get_particle(event_id, i).px;
            double py_local = particle_list->get_particle(event_id, i).py;
            double p_perp = sqrt(px_local*px_local + py_local*py_local);
            if (p_perp > pT_min && p_perp < pT_max) {
                count++;
            }
        }
    }
    number_of_events[count-1]++;
}

void particle_yield_distribution::output_particle_yield_distribution() {
    // this function outputs the particle yield distribution

    ostringstream filename;
    filename << path << "/particle_" << particle_monval
             << "_yield_distribution.dat";
    ofstream output(filename.str().c_str());

    // output header
    output << "# N  P(N)" << endl;

    for (int i = 0; i < n_max; i++) {
        double p_N = (number_of_events[i]/total_number_of_events
                      /reconst_branching_ratio);
        output << scientific << setw(18) << setprecision(8)
               << i << "   " << p_N << endl;
    }
    output.close();
}

