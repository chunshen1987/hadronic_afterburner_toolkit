// Copyright Chun Shen @ 2016

#include <sys/time.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include "./particle_decay.h"

using namespace std;

particle_decay::particle_decay() {
    resonance_table.clear();
}

particle_decay::~particle_decay() {
    for (unsigned int i = 0; i < resonance_table.size(); i++) {
        for (int j = 0; j < resonance_table[i]->decays; j++) {
            delete resonance_table[i]->decay_channels[j];
        }
        delete resonance_table[i];
    }
    resonance_table.clear();
    timeval a;
    gettimeofday(&a, 0);
    int random_seed = a.tv_usec;
    srand48(random_seed);
}

int particle_decay::read_resonances_list() {
    double eps = 1e-15;
    cout << " -- Read in particle resonance decay table..." << endl;
    ifstream resofile("EOS/pdg.dat");
    int dummy_int;
    int particle_monval;
    resofile >> particle_monval;
    while (!resofile.eof()) {
        // add one resonance
        particle_decay_info *temp_resonance = new particle_decay_info;
        temp_resonance->monval = particle_monval;
        resofile >> temp_resonance->name;
        resofile >> temp_resonance->mass;
        resofile >> temp_resonance->width;
        resofile >> temp_resonance->gspin;        // spin degeneracy
        resofile >> temp_resonance->baryon;
        resofile >> temp_resonance->strange;
        resofile >> temp_resonance->charm;
        resofile >> temp_resonance->bottom;
        resofile >> temp_resonance->gisospin;     // isospin degeneracy
        resofile >> temp_resonance->charge;
        resofile >> temp_resonance->decays;
        for (int j = 0; j < temp_resonance->decays; j++) {
            // read in its decay channels
            decay_channel_info *temp_decay_channel = new decay_channel_info;
            resofile >> dummy_int;
            resofile >> temp_decay_channel->decay_Npart;
            resofile >> temp_decay_channel->branching_ratio;
            resofile >> temp_decay_channel->decay_part[0];
            resofile >> temp_decay_channel->decay_part[1];
            resofile >> temp_decay_channel->decay_part[2];
            resofile >> temp_decay_channel->decay_part[3];
            resofile >> temp_decay_channel->decay_part[4];
            temp_resonance->decay_channels.push_back(temp_decay_channel);
        }

        // decide whether particle is stable under strong interactions
        if (temp_resonance->decay_channels[0]->decay_Npart == 1) {
            temp_resonance->stable = 1;
        } else {
            temp_resonance->stable = 0;
        }
        resonance_table.push_back(temp_resonance);

        if (temp_resonance->baryon == 1) {
            // add anti-particle entry for baryon
            particle_decay_info *temp_anti_resonance = new particle_decay_info;
            temp_anti_resonance->monval = -temp_resonance->monval;
            ostringstream antiname;
            antiname << "Anti-" << temp_resonance->name;
            temp_anti_resonance->name = antiname.str();
            temp_anti_resonance->mass = temp_resonance->mass;
            temp_anti_resonance->width = temp_resonance->width;
            temp_anti_resonance->gspin = temp_resonance->gspin;
            temp_anti_resonance->baryon = -temp_resonance->baryon;
            temp_anti_resonance->strange = -temp_resonance->strange;
            temp_anti_resonance->charm = -temp_resonance->charm;
            temp_anti_resonance->bottom = -temp_resonance->bottom;
            temp_anti_resonance->gisospin = temp_resonance->gisospin;
            temp_anti_resonance->charge = -temp_resonance->charge;
            temp_anti_resonance->decays = temp_resonance->decays;
            temp_anti_resonance->stable = temp_resonance->stable;
            for (int j = 0; j < temp_resonance->decays; j++) {
                // add anti-particle decay channels
                decay_channel_info *temp_anti_decay_channel =
                                                    new decay_channel_info;
                temp_anti_decay_channel->decay_Npart = 
                        temp_resonance->decay_channels[j]->decay_Npart;
                temp_anti_decay_channel->branching_ratio = 
                        temp_resonance->decay_channels[j]->branching_ratio;
                for (int k = 0; k < 5; k++) {
                    int decay_part_monval =
                            temp_resonance->decay_channels[j]->decay_part[k];
                    if (decay_part_monval == 0) {
                        // a null entry
                        temp_anti_decay_channel->decay_part[k] = 0;
                    } else {
                        // find the index for decay particle in the
                        // current resonance table
                        int idx;
                        int reso_table_length = resonance_table.size();
                        for (idx = 0; idx < reso_table_length; idx++) {
                            if (resonance_table[idx]->monval
                                == decay_part_monval) {
                                break;
                            }
                        }
                        double temp_br =
                            temp_resonance->decay_channels[j]->branching_ratio;
                        if (idx == reso_table_length
                            && temp_resonance->stable == 0 && temp_br > eps) {
                            cout << "Error: can not find decay particle index "
                                 << "for anti-baryon!" << endl;
                            cout << "particle monval : " 
                                 << decay_part_monval << endl;
                            exit(1);
                        }
                        if (resonance_table[idx]->baryon == 0
                            && resonance_table[idx]->charge == 0 
                            && resonance_table[idx]->strange == 0) {
                            temp_anti_decay_channel->decay_part[k] =
                                                        decay_part_monval;
                        } else {
                            temp_anti_decay_channel->decay_part[k] =
                                                        - decay_part_monval;
                        }
                    }
                }
                temp_anti_resonance->decay_channels.push_back(
                                                    temp_anti_decay_channel);
            }
            resonance_table.push_back(temp_anti_resonance);
        }
        resofile >> particle_monval;
    }
    resofile.close();  // close the pdg file
    int Nparticle = resonance_table.size();
    for (int i = 0; i < Nparticle; i++) {
        // determine the quantum statistics for bosons and fermions
        if (resonance_table[i]->baryon == 0) {
            resonance_table[i]->sign = -1;
        } else {
            resonance_table[i]->sign = 1;
        }
    }
    return(Nparticle);
}

void particle_decay::check_resonance_table() {
    for (unsigned int i = 0; i < resonance_table.size(); i++) {
        cout << "name: " << resonance_table[i]->name << endl;
        cout << "monval: " << resonance_table[i]->monval << endl;
        cout << "mass: " << resonance_table[i]->mass << endl;
        cout << "stable: " << resonance_table[i]->stable << endl;
        cout << "# of decay channels: " << resonance_table[i]->decays << endl;
        for (int j = 0; j < resonance_table[i]->decays; j++) {
            cout << j << ": # of daughters: "
                 << resonance_table[i]->decay_channels[j]->decay_Npart << endl;
            cout << j << ": branching_ratio: "
                 << resonance_table[i]->decay_channels[j]->branching_ratio
                 << endl;
            cout << j << ": decay part: ";
            for (int k = 0;
                 k < resonance_table[i]->decay_channels[j]->decay_Npart;
                 k++) {
                 cout << resonance_table[i]->decay_channels[j]->decay_part[k]
                      << "  ";
            }
            cout << endl;
        }
    }
}

double particle_decay::get_particle_width(particle_info *part) {
    double width = 0.0;
    for (unsigned int i = 0; i < resonance_table.size(); i++) {
        if (part->monval == resonance_table[i]->monval) {
            width = resonance_table[i]->width;
            break;
        }
    }
    return(width);
}

void particle_decay::perform_two_body_decay(particle_info *mother,
                                            particle_info *daughter1,
                                            particle_info *daughter2) {
    // this function perform two body decay
    double M_pole = mother->mass;
    double M_width = get_particle_width(mother);
    double m1 = daughter1->mass;
    double m2 = daughter2->mass;
    double M_min = m1 + m2;
    if (M_pole < (m1 + m2)) {
        cout << "Error:particleSamples::perform_two_body_decay:"
             << "can not found decays!" << endl;
        cout << "M = " << M_pole << ", m1 = " << m1 << ", m2 = " << m2 << endl;
        exit(1);
    }
    double M_sampled = sample_breit_wigner(M_pole, M_width, M_min);
    double temp = M_sampled*M_sampled - m1*m1 - m2*m2;
    double p_lrf = sqrt(temp*temp - 4*m1*m1*m2*m2)/(2*M_sampled);

    // randomly pick emission angle
    double phi = drand48()*2*M_PI;
    double cos_theta = 2.*(drand48() - 0.5);
    double sin_theta = sqrt(1. - cos_theta*cos_theta);

    // compute daughter particles' energy and momentum in the rest frame
    double E1_lrf = sqrt(p_lrf*p_lrf + m1*m1);
    double p1_lrf_x = p_lrf*sin_theta*cos(phi);
    double p1_lrf_y = p_lrf*sin_theta*sin(phi);
    double p1_lrf_z = p_lrf*cos_theta;
    double E2_lrf = sqrt(p_lrf*p_lrf + m2*m2);
    double p2_lrf_x = -p1_lrf_x;
    double p2_lrf_y = -p1_lrf_y;
    double p2_lrf_z = -p1_lrf_z;
    
    // compute mother velocity
    double vx = mother->px/mother->E;
    double vy = mother->py/mother->E;
    double vz = mother->pz/mother->E;

    // perform the boost
    double v2 = vx*vx + vy*vy + vz*vz;
    double gamma = 1./sqrt(1. - v2);
    double gamma_m_1 = gamma - 1.;
    double vp1 = vx*p1_lrf_x + vy*p1_lrf_y + vz*p1_lrf_z;
    double vp2 = vx*p2_lrf_x + vy*p2_lrf_y + vz*p2_lrf_z;
    daughter1->E = gamma*(E1_lrf - vp1);
    daughter1->px = p1_lrf_x + (gamma_m_1*vp1/v2 - gamma*E1_lrf)*vx;
    daughter1->py = p1_lrf_y + (gamma_m_1*vp1/v2 - gamma*E1_lrf)*vy;
    daughter1->pz = p1_lrf_z + (gamma_m_1*vp1/v2 - gamma*E1_lrf)*vz;
    daughter2->E = gamma*(E2_lrf - vp2);
    daughter2->px = p2_lrf_x + (gamma_m_1*vp2/v2 - gamma*E2_lrf)*vx;
    daughter2->py = p2_lrf_y + (gamma_m_1*vp2/v2 - gamma*E2_lrf)*vy;
    daughter2->pz = p2_lrf_z + (gamma_m_1*vp2/v2 - gamma*E2_lrf)*vz;

    double life_time = 1e10;
    if (M_width > 1e-10) {
        // compute life-time = gamma*1/\Gamma
        double tau0 = mother->E/(M_sampled)*1./(M_width);
        life_time = -tau0*log(drand48());
    }

    daughter1->t = mother->t + life_time;
    daughter1->x = mother->x + mother->x/mother->E*life_time;
    daughter1->y = mother->y + mother->py/mother->E*life_time;
    daughter1->z = mother->z + mother->pz/mother->E*life_time;
    daughter2->t = mother->t + life_time;
    daughter2->x = mother->x + mother->px/mother->E*life_time;
    daughter2->y = mother->y + mother->py/mother->E*life_time;
    daughter2->z = mother->z + mother->pz/mother->E*life_time;

    return;
}

void particle_decay::perform_three_body_decay(particle_info *mother,
                                              particle_info *daughter1,
                                              particle_info *daughter2,
                                              particle_info *daughter3) {
    // this function perform 3 body decays
    double M_pole = mother->mass;
    double M_width = get_particle_width(mother);
    double m1 = daughter1->mass;
    double m2 = daughter2->mass;
    double m3 = daughter3->mass;
    double M_min = m1 + m2 + m3;
    if (M_pole < M_min) {
        cout << "Error:particleSamples::perform_two_body_decay:"
             << "can not found decays!" << endl;
        cout << "M = " << M_pole << ", m1 = " << m1 << ", m2 = " << m2
             << ", m3 = " << m3 << endl;
        exit(1);
    }
    double M_sampled = sample_breit_wigner(M_pole, M_width, M_min);
    // generate lrf E1, E2, and theta12 using accept and reject method
    double E1_lrf, E2_lrf, E3_lrf, p1_lrf, p2_lrf, cos12_lrf;
    do {
        do {
            E1_lrf = drand48()*(M_sampled - m1 - m2 - m3) + m1;
            E2_lrf = drand48()*(M_sampled - m1 - m2 - m3) + m2;
        } while (E1_lrf + E2_lrf > M_sampled);
        p1_lrf = sqrt(E1_lrf*E1_lrf - m1*m1);
        p2_lrf = sqrt(E2_lrf*E2_lrf - m2*m2);
        double E3_lrf = M_sampled - E1_lrf - E2_lrf;
        cos12_lrf = (E3_lrf*E3_lrf - p1_lrf*p1_lrf - p2_lrf*p2_lrf - m3*m3)
                     /(2.*p1_lrf*p2_lrf);
    } while (cos12_lrf < - 1.0 || cos12_lrf > 1.0);
    // now we get the a good sample

    // sample the lifetime
    double life_time = 1e10;
    if (M_width > 1e-10) {
        double tau = mother->E/(M_sampled)*1./M_width;
        life_time = -tau*log(drand48());
    }
    // compute the decay position
    double decay_time = mother->t + life_time;
    double decay_x = mother->x + mother->px/mother->E*life_time;
    double decay_y = mother->y + mother->py/mother->E*life_time;
    double decay_z = mother->z + mother->pz/mother->E*life_time;

    // compute the momentum of decay daughters
    double tp2_lrf_x = p2_lrf*sqrt(1. - cos12_lrf*cos12_lrf);
    double tp2_lrf_z = p2_lrf*cos12_lrf;
    double tp3_lrf_x = - tp2_lrf_x;
    double tp3_lrf_z = - (p1_lrf + tp2_lrf_z); 
    double phi = 2.*M_PI*drand48();
    double ksi = 2.*M_PI*drand48();
    double cos_theta = 2.*drand48() - 1.0;

    double sin_phi = sin(phi);
    double cos_phi = cos(phi);
    double sin_ksi = sin(ksi);
    double cos_ksi = cos(ksi);
    double sin_theta = sqrt(1. - cos_theta*cos_theta);

    double p1_lrf_x = - p1_lrf*sin_theta*cos_ksi;
    double p1_lrf_y = p1_lrf*sin_theta*sin_ksi;
    double p1_lrf_z = p1_lrf*cos_theta;
    E1_lrf = sqrt(m1*m1 + p1_lrf_x*p1_lrf_x
                  + p1_lrf_y*p1_lrf_y + p1_lrf_z*p1_lrf_z);
    double p2_lrf_x = (tp2_lrf_x*(cos_phi*cos_theta*cos_ksi - sin_phi*sin_ksi)
                       - tp2_lrf_z*sin_theta*cos_ksi);
    double p2_lrf_y = (tp2_lrf_x*(-cos_phi*cos_theta*sin_ksi - sin_phi*cos_ksi)
                       + tp2_lrf_z*sin_theta*sin_ksi);
    double p2_lrf_z = tp2_lrf_x*(cos_phi*sin_theta) + tp2_lrf_z*cos_theta;
    E2_lrf = sqrt(m2*m2 + p2_lrf_x*p2_lrf_x
                  + p2_lrf_y*p2_lrf_y + p2_lrf_z*p2_lrf_z);
    double p3_lrf_x = (tp3_lrf_x*(cos_phi*cos_theta*cos_ksi - sin_phi*sin_ksi)
                       - tp3_lrf_z*sin_theta*cos_ksi);
    double p3_lrf_y = (tp3_lrf_x*(-cos_phi*cos_theta*sin_ksi - sin_phi*cos_ksi)
                       + tp3_lrf_z*(sin_theta*sin_ksi));
    double p3_lrf_z = tp3_lrf_x*cos_phi*sin_theta + tp3_lrf_z*cos_theta;
    E3_lrf = sqrt(m3*m3 + p3_lrf_x*p3_lrf_x
                  + p3_lrf_y*p3_lrf_y + p3_lrf_z*p3_lrf_z);

    double vx = mother->px/mother->E;
    double vy = mother->py/mother->E;
    double vz = mother->pz/mother->E;
    double v2 = vx*vx + vy*vy + vz*vz;
    double gamma = 1./sqrt(1. - v2);
    double gamma_m_1 = gamma - 1.;
    double vp1 = vx*p1_lrf_x + vy*p1_lrf_y + vz*p1_lrf_z;
    double vp2 = vx*p2_lrf_x + vy*p2_lrf_y + vz*p2_lrf_z;
    double vp3 = vx*p3_lrf_x + vy*p3_lrf_y + vz*p3_lrf_z;

    daughter1->E = gamma*(E1_lrf - vp1);
    daughter1->px = p1_lrf_x + (gamma_m_1*vp1/v2 - gamma*E1_lrf)*vx;
    daughter1->py = p1_lrf_y + (gamma_m_1*vp1/v2 - gamma*E1_lrf)*vy;
    daughter1->pz = p1_lrf_z + (gamma_m_1*vp1/v2 - gamma*E1_lrf)*vz;
    daughter1->t = decay_time;
    daughter1->x = decay_x;
    daughter1->y = decay_y;
    daughter1->z = decay_z;

    daughter2->E = gamma*(E2_lrf - vp2);
    daughter2->px = p2_lrf_x + (gamma_m_1*vp2/v2 - gamma*E2_lrf)*vx;
    daughter2->py = p2_lrf_y + (gamma_m_1*vp2/v2 - gamma*E2_lrf)*vy;
    daughter2->pz = p2_lrf_z + (gamma_m_1*vp2/v2 - gamma*E2_lrf)*vz;
    daughter2->t = decay_time;
    daughter2->x = decay_x;
    daughter2->y = decay_y;
    daughter2->z = decay_z;

    daughter3->E = gamma*(E3_lrf - vp3);
    daughter3->px = p3_lrf_x + (gamma_m_1*vp3/v2 - gamma*E3_lrf)*vx;
    daughter3->py = p3_lrf_y + (gamma_m_1*vp3/v2 - gamma*E3_lrf)*vy;
    daughter3->pz = p3_lrf_z + (gamma_m_1*vp3/v2 - gamma*E3_lrf)*vz;
    daughter3->t = decay_time;
    daughter3->x = decay_x;
    daughter3->y = decay_y;
    daughter3->z = decay_z;

    return;
}

double particle_decay::sample_breit_wigner(double mass, double width,
                                           double M_min) {
    // this function sample the Breit Wigner distribution for the mass
    // particle mass sampled is >= M_min
    // From Wiki: https://en.wikipedia.org/wiki/Cauchy_distribution
    // The CDF of BreitWigner (or Cauchy) distribution is
    // y = 1/pi arctan((x-x0)/gamma) + 1/2
    // The inverse CDF is 
    // x = x0 + gamma*tan(pi*(y - 1/2))

    // compute the minimum probability given by M_min
    double p_min = atan2((M_min - mass), width/2.)/M_PI + 0.5;
    // generate a random number from p_min to 1
    double y = (1. - p_min)*drand48() + p_min;
    // compute the corresponding sampled mass
    double mass_sampled = mass + width/2.*tan(M_PI*(y - 0.5));
    return(mass_sampled);
}
