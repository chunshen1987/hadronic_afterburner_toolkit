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
                ParameterReader &paraRdr, std::string path) :
        paraRdr_(paraRdr), path_(path) {
    particle_monval = paraRdr_.getVal("particle_monval");
    if (particle_monval == 333) {
        // phi(1020) is reconstructed from (K^+, K^-) pairs
        reconst_branching_ratio = 0.489;
    } else {
        reconst_branching_ratio = 1.0;
    }

    net_particle_flag = paraRdr_.getVal("net_particle_flag");
    n_max = 6000;
    number_of_events = new int[n_max];
    for (int i = 0; i < n_max; i++)
        number_of_events[i] = 0;

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
    delete [] number_of_events;
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
    int count_anti_particle = 0;
    if (net_particle_flag == 1) {
        int number_of_anti_particles =
            particle_list->get_number_of_anti_particles(event_id);
        for (int i = 0; i < number_of_anti_particles; i++) {
            double pz_local = particle_list->get_anti_particle(event_id, i).pz;
            double E_local = particle_list->get_anti_particle(event_id, i).E;

            double rap_local;
            if (rap_type == 0) {
                double mass = particle_list->get_anti_particle(event_id, i).mass;
                double pmag = sqrt(E_local*E_local - mass*mass);
                rap_local = 0.5*log((pmag + pz_local)/(pmag - pz_local));
            } else {
                rap_local = 0.5*log((E_local + pz_local)/(E_local - pz_local));
            }

            if (rap_local > rap_min && rap_local < rap_max) {
                double px_local = particle_list->get_anti_particle(event_id, i).px;
                double py_local = particle_list->get_anti_particle(event_id, i).py;
                double p_perp = sqrt(px_local*px_local + py_local*py_local);
                if (p_perp > pT_min && p_perp < pT_max) {
                    count_anti_particle++;
                }
            }
        }
    }
    int idx = count - count_anti_particle + n_max/2;
    if (idx < 0 || idx > n_max - 1) {
        cout << "particle number is out of bound!" << endl;
        cout << "bound: [" << -n_max/2 << " , " << n_max/2 << "]" << endl;
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
        double p_N = (static_cast<double>(number_of_events[i])
                      /static_cast<double>(total_number_of_events)
                      /reconst_branching_ratio);
        output << scientific << setw(18) << setprecision(8)
               << i - n_max/2 << "   " << p_N << endl;
    }
    output.close();
}

