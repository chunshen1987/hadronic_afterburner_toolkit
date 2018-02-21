// Copyright Chun Shen @ 2016
#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <cstdlib>
#include <vector>

#include "zlib.h"
#include "./particleSamples.h"

using namespace std;

particleSamples::particleSamples(ParameterReader* paraRdr_in, string path_in) {
    paraRdr = paraRdr_in;
    path = path_in;

    echo_level = paraRdr->getVal("echo_level");
    event_buffer_size = paraRdr->getVal("event_buffer_size");
    read_in_mode = paraRdr->getVal("read_in_mode");
    run_mode = paraRdr->getVal("run_mode");
    if (run_mode == 1) {
        reject_decay_flag = paraRdr->getVal("reject_decay_flag");
        tau_reject = paraRdr->getVal("tau_reject");
    } else {
        reject_decay_flag = 0;
        tau_reject = 10000.;
    }

    if (run_mode == 2) {
        net_particle_flag = paraRdr->getVal("net_particle_flag");
        if (net_particle_flag == 1) {
            anti_particle_list = new vector< vector<particle_info>* >;
        }
    } else {
        net_particle_flag = 0;
    }
    
    // read in particle Monte-Carlo number
    particle_monval = paraRdr->getVal("particle_monval");
    flag_isospin = paraRdr->getVal("distinguish_isospin");
    if (read_in_mode == 1 || read_in_mode == 2
            || read_in_mode == 3 || read_in_mode == 4) {
        get_UrQMD_id(particle_monval);
    }

    resonance_feed_down_flag = paraRdr->getVal("resonance_feed_down_flag");
    select_resonances_flag = 0;
    if (resonance_feed_down_flag == 1) {
        resonance_list = new vector< vector<particle_info>* >;
        decayer_ptr = new particle_decay;
        select_resonances_flag = paraRdr->getVal("select_resonances_flag");
        initialize_selected_resonance_list();
    }

    particle_list = new vector< vector<particle_info>* >;
    particle_list_mixed_event = new vector< vector<particle_info>* >;
    if (abs(particle_monval) == 3122) {
        // for Lambda and anti-Lambda
        resonance_weak_feed_down_flag = paraRdr->getVal(
                                            "resonance_weak_feed_down_flag");
        if (resonance_weak_feed_down_flag == 1) {
            // include Sigma0 feed down to Lambda
            resonance_list = new vector< vector<particle_info>* >;
            decayer_ptr = new particle_decay;
        }
    } else {
        resonance_weak_feed_down_flag = 0;
    }

    if (particle_monval == 333) {
        // for phi(1020) meson, we need reconstruct them from (K^+ K^-) pairs
        reconst_flag = 1;
        reconst_list_1 = new vector< vector<particle_info>* >;
        reconst_list_2 = new vector< vector<particle_info>* >;
    } else {
        reconst_flag = 0;
    }

    flag_charge_dependence = 0;
    if (particle_monval == 9999) {
        flag_charge_dependence = paraRdr->getVal("flag_charge_dependence");
        if (flag_charge_dependence == 1) {
            positive_charge_hadron_list = new vector< vector<particle_info>* >;
            negative_charge_hadron_list = new vector< vector<particle_info>* >;
        }
    }

    ostringstream filename;
    ostringstream filename_mixed_event;
    if (read_in_mode == 0) {
        filename << path << "/OSCAR.DAT";
        filename_mixed_event << path << "/OSCAR_mixed_event.DAT";
    } else if (read_in_mode == 1 || read_in_mode == 2 || read_in_mode == 3
               || read_in_mode == 4 || read_in_mode == 5) {
        filename << path << "/particle_list.dat";
        filename_mixed_event << path << "/particle_list_mixed_event.dat";
    } else if (read_in_mode == 10) {
        filename << path << "/particle_samples.gz";
        filename_mixed_event << path << "/particle_samples_mixed_event.gz";
    }

    if (read_in_mode != 2 && read_in_mode != 10) {
        inputfile.open(filename.str().c_str());
        if (!inputfile.is_open()) {
            cout << "particleSamples:: Error: input file: " << filename.str() 
                 << " can not open!" << endl;
            exit(1);
        }
        if (run_mode == 1) {
            inputfile_mixed_event.open(filename_mixed_event.str().c_str());
            if (!inputfile_mixed_event.is_open()) {
                cout << "particleSamples:: Error: input file: " 
                     << filename_mixed_event.str() << " can not open!" << endl;
                exit(1);
            }
        }
    } else {
        inputfile_gz = gzopen(filename.str().c_str(), "rb");
        if (!inputfile_gz) {
            cout << "particleSamples:: Error: input file: " << filename.str() 
                 << " can not open!" << endl;
            exit(1);
        }
        if (run_mode == 1) {
            inputfile_mixed_event_gz =
                            gzopen(filename_mixed_event.str().c_str(), "rb");
            if (!inputfile_mixed_event_gz) {
                cout << "particleSamples:: Error: input file: " 
                     << filename_mixed_event.str() << " can not open!" << endl;
                exit(1);
            }
        }
    }

    // skip the header file for OSCAR
    string temp;
    if (read_in_mode == 0) {
        getline(inputfile, temp);
        getline(inputfile, temp);
        getline(inputfile, temp);
        if (run_mode == 1) {
            getline(inputfile_mixed_event, temp);
            getline(inputfile_mixed_event, temp);
            getline(inputfile_mixed_event, temp);
        }
    }
    
    // skip the header in JAM
    if (read_in_mode == 5) {
        getline(inputfile, temp);
        if (run_mode == 1) {
            getline(inputfile_mixed_event, temp);
        }
    }

    initialize_charged_hadron_pdg_list();
    initialize_charged_hadron_urqmd_id_list();
}

particleSamples::~particleSamples() {
    if (read_in_mode != 2 && read_in_mode != 10) {
        inputfile.close();
    } else {
        gzclose(inputfile_gz);
    }
    if (run_mode == 1) {
        if (read_in_mode != 2 && read_in_mode != 10) {
            inputfile_mixed_event.close();
        } else {
            gzclose(inputfile_mixed_event_gz);
        }
    }
    
    for (unsigned int i = 0; i < particle_list->size(); i++)
        (*particle_list)[i]->clear();
    particle_list->clear();
    delete particle_list;
    for (unsigned int i = 0; i < particle_list_mixed_event->size(); i++)
        (*particle_list_mixed_event)[i]->clear();
    particle_list_mixed_event->clear();
    delete particle_list_mixed_event;

    if (net_particle_flag == 1) {
        for (unsigned int i = 0; i < anti_particle_list->size(); i++) {
            (*anti_particle_list)[i]->clear();
        }
        anti_particle_list->clear();
        delete anti_particle_list;
    }

    if (resonance_weak_feed_down_flag == 1) {
        for (unsigned int i = 0; i < resonance_list->size(); i++)
            (*resonance_list)[i]->clear();
        resonance_list->clear();
        delete resonance_list;
        delete decayer_ptr;
    }

    if (resonance_feed_down_flag == 1) {
        for (unsigned int i = 0; i < resonance_list->size(); i++)
            (*resonance_list)[i]->clear();
        resonance_list->clear();
        delete resonance_list;
        delete decayer_ptr;
        if (select_resonances_flag == 1)
            select_resonances_list.clear();
    }

    if (reconst_flag == 1) {
        for (unsigned int i = 0; i < reconst_list_1->size(); i++)
            (*reconst_list_1)[i]->clear();
        reconst_list_1->clear();
        delete reconst_list_1;
        for (unsigned int i = 0; i < reconst_list_2->size(); i++)
            (*reconst_list_2)[i]->clear();
        reconst_list_2->clear();
        delete reconst_list_2;
    }

    if (flag_charge_dependence == 1) {
        for (unsigned int i = 0; i < positive_charge_hadron_list->size(); i++)
            (*positive_charge_hadron_list)[i]->clear();
        positive_charge_hadron_list->clear();
        delete positive_charge_hadron_list;
        for (unsigned int i = 0; i < negative_charge_hadron_list->size(); i++)
            (*negative_charge_hadron_list)[i]->clear();
        negative_charge_hadron_list->clear();
        delete negative_charge_hadron_list;
    }
}

void particleSamples::initialize_charged_hadron_pdg_list() {
    charged_hadron_pdg_list[0] = 211;       // pion
    charged_hadron_pdg_list[1] = 321;       // kaon
    charged_hadron_pdg_list[2] = 2212;      // proton
    charged_hadron_pdg_list[3] = 3222;      // Sigma^+
    charged_hadron_pdg_list[4] = 3112;      // Sigma^-
    charged_hadron_pdg_list[5] = 3312;      // Xi^-
}

void particleSamples::initialize_charged_hadron_urqmd_id_list() {
    charged_hadron_urqmd_id_list[0] = 101;   // pion
    charged_hadron_urqmd_id_list[1] = 106;   // kaon
    charged_hadron_urqmd_id_list[2] = 1;     // proton
    charged_hadron_urqmd_id_list[3] = 40;    // Sigma^+ and Sigma^-
    charged_hadron_urqmd_id_list[4] = 49;    // Xi^-
}

void particleSamples::initialize_baryon_urqmd_id_list() {
    baryon_urqmd_id_list[0] = 1;     // proton and neturon
    baryon_urqmd_id_list[1] = 40;    // Sigma^+, Sigma^0, and Sigma^-
    baryon_urqmd_id_list[2] = 49;    // Xi^-
    baryon_urqmd_id_list[3] = 27;    // Lambda
    baryon_urqmd_id_list[4] = 55;    // Omega
}

void particleSamples::initialize_selected_resonance_list() {
    ostringstream filename;
    filename << "EOS/selected_resonances_list.dat";
    ifstream reso_file(filename.str().c_str());
    if (!reso_file.is_open()) {
        cout << "[Error]:particleSamples::initialize_selected_resonance_list:"
             << "Can not find the selected resonance list file: "
             << filename.str() << endl;
        exit(1);
    }
    int temp_id;
    reso_file >> temp_id;
    while (!reso_file.eof()) {
        select_resonances_list.push_back(temp_id);
        reso_file >> temp_id;
    }
    reso_file.close();
    if (echo_level > 8) {
        cout << "[Info]:particleSamples::initialize_selected_resonance_list:"
             << "selected resonance list: " << endl;
        for (unsigned int ireso = 0; ireso < select_resonances_list.size();
                ireso++) {
            cout << select_resonances_list[ireso] << endl;
        }
    }
}

void particleSamples::get_UrQMD_id(int monval) {
    // find the corresponding UrQMD id number
    if (monval == 211) {
        // pion^+
        particle_urqmd_id = 101;
        particle_urqmd_isospin = 2;
    } else if (monval == -211) {
        // pion^-
        particle_urqmd_id = 101;
        particle_urqmd_isospin = -2;
    } else if (monval == 111) {
        // pion^0
        particle_urqmd_id = 101;
        particle_urqmd_isospin = 0;
    } else if (monval == 321) {
        // Kaon^+
        particle_urqmd_id = 106;
        particle_urqmd_isospin = 1;
    } else if (monval == 311) {
        // Kaon^0
        particle_urqmd_id = 106;
        particle_urqmd_isospin = -1;
    } else if (monval == -321) {
        // Kaon^-
        particle_urqmd_id = -106;
        particle_urqmd_isospin = -1;
    } else if (monval == -311) {
        // anti-Kaon^0
        particle_urqmd_id = -106;
        particle_urqmd_isospin = 1;
    } else if (monval == 2212) {
        // proton
        particle_urqmd_id = 1;
        particle_urqmd_isospin = 1;
    } else if (monval == -2212) {
        // anti-proton
        particle_urqmd_id = -1;
        particle_urqmd_isospin = -1;
    } else if (monval == 3222) {
        // Sigma^+
        particle_urqmd_id = 40;
        particle_urqmd_isospin = 2;
    } else if (monval == -3222) {
        // anti-Sigma^+
        particle_urqmd_id = -40;
        particle_urqmd_isospin = -2;
    } else if (monval == 3212) {
        // Sigma^0
        particle_urqmd_id = 40;
        particle_urqmd_isospin = 0;
    } else if (monval == -3212) {
        // anti-Sigma^0
        particle_urqmd_id = -40;
        particle_urqmd_isospin = 0;
    } else if (monval == 3112) {
        // Sigma^-
        particle_urqmd_id = 40;
        particle_urqmd_isospin = -2;
    } else if (monval == -3112) {
        // anti-Sigma^-
        particle_urqmd_id = -40;
        particle_urqmd_isospin = 2;
    } else if (monval == 3312) {
        // Xi^-
        particle_urqmd_id = 49;
        particle_urqmd_isospin = -1;
    } else if (monval == -3312) {
        // anti-Xi^-
        particle_urqmd_id = -49;
        particle_urqmd_isospin = 1;
    } else if (monval == 3122) {
        // Lambda
        particle_urqmd_id = 27;
        particle_urqmd_isospin = 0;
    } else if (monval == -3122) {
        // anti-Lambda
        particle_urqmd_id = -27;
        particle_urqmd_isospin = 0;
    } else if (monval == 3334) {
        // Omega
        particle_urqmd_id = 55;
        particle_urqmd_isospin = 0;
    } else if (monval == -3334) {
        // anti-Omega
        particle_urqmd_id = -55;
        particle_urqmd_isospin = 0;
    } else if (monval == 333) {
        // phi(1020) meson
        particle_urqmd_id = 109;
        particle_urqmd_isospin = 0;
    } else if (monval == 221) {
        // eta meson
        particle_urqmd_id = 102;
        particle_urqmd_isospin = 0;
    } else if (monval == 22) {
        // photons
        particle_urqmd_id = 100;
        particle_urqmd_isospin = 0;
    } else if (monval == 9996) {
        // all strange hadrons
        particle_urqmd_id = 9996;
        particle_urqmd_isospin = 0;
    } else if (monval == -9996) {
        // all anti-strange hadrons
        particle_urqmd_id = -9996;
        particle_urqmd_isospin = 0;
    } else if (monval == 9997) {
        // all baryons
        particle_urqmd_id = 9997;
        particle_urqmd_isospin = 0;
    } else if (monval == -9997) {
        // all anti-baryons
        particle_urqmd_id = -9997;
        particle_urqmd_isospin = 0;
    } else if (monval == 9998) {
        // all positive charged hadrons
        particle_urqmd_id = 9998;
        particle_urqmd_isospin = 0;
    } else if (monval == -9998) {
        // all negative charged hadrons
        particle_urqmd_id = -9998;
        particle_urqmd_isospin = 0;
    } else if (monval == 9999) {
        // all charged hadrons
        particle_urqmd_id = 9999;
        particle_urqmd_isospin = 0;
    } else if (monval == -1) {
        // all hadrons
        particle_urqmd_id = -1;
        particle_urqmd_isospin = 0;
    }
}

int particleSamples::get_pdg_id(int urqmd_id, int urqmd_isospin) {
    int monval = 0;
    if (urqmd_id == 40 && urqmd_isospin == 0) {
        monval = 3212;
    } else if (urqmd_id == -40 && urqmd_isospin == 0) {
        monval = -3212;
    }
    return(monval);
}

int particleSamples::read_in_particle_samples() {
    if (read_in_mode == 0) {
        read_in_particle_samples_OSCAR();
        resonance_weak_feed_down_flag = 0;
        reconst_flag = 0;
    } else if (read_in_mode == 1) {
        read_in_particle_samples_UrQMD();
    } else if (read_in_mode == 2) {
        read_in_particle_samples_UrQMD_zipped();
    } else if (read_in_mode == 3) {
        read_in_particle_samples_Sangwook();
    } else if (read_in_mode == 4) {
        read_in_particle_samples_UrQMD_3p3();
    } else if (read_in_mode == 5) {
        read_in_particle_samples_JAM();
    } else if (read_in_mode == 10) {
        read_in_particle_samples_gzipped();
        resonance_weak_feed_down_flag = 0;
        reconst_flag = 0;
    }

    if (resonance_weak_feed_down_flag == 1) {
        perform_weak_resonance_feed_down();
    }

    if (resonance_feed_down_flag == 1) {
        perform_resonance_feed_down(particle_list);
    }

    if (reconst_flag == 1) {
        perform_particle_reconstruction();
    }

    return(0);
}

int particleSamples::read_in_particle_samples_mixed_event() {
    if (read_in_mode == 0) {
        read_in_particle_samples_OSCAR_mixed_event();
        resonance_weak_feed_down_flag = 0;
    } else if (read_in_mode == 1) {
        read_in_particle_samples_UrQMD_mixed_event();
    } else if (read_in_mode == 2) {
        read_in_particle_samples_UrQMD_mixed_event_zipped();
    } else if (read_in_mode == 3) {
        read_in_particle_samples_mixed_event_Sangwook();
    } else if (read_in_mode == 4) {
        read_in_particle_samples_UrQMD_3p3_mixed_event();
    } else if (read_in_mode == 5) {
        read_in_particle_samples_JAM_mixed_event();
    } else if (read_in_mode == 10) {
        read_in_particle_samples_mixed_event_gzipped();
    }

    if (resonance_weak_feed_down_flag == 1) {
        perform_weak_resonance_feed_down();
    }
    
    if (resonance_feed_down_flag == 1) {
        perform_resonance_feed_down(particle_list_mixed_event);
    }

    return(0);
}

int particleSamples::decide_to_pick_UrQMD(int pid, int iso3, int charge,
                                          int parent_proc_type) {
    int pick_flag = 0;
    if (resonance_feed_down_flag == 1) {  // all hadrons
        pick_flag = 1;
    } else if (particle_urqmd_id == 9999) {  // charged hadrons
        int in_flag = 0;
        for (int i = 0; i < 5; i++) {
            if (abs(pid) == charged_hadron_urqmd_id_list[i]) {
                in_flag = 1;
                break;
            }
        }
        if (in_flag == 1 && charge != 0)
            pick_flag = 1;
    } else if (particle_urqmd_id == 9998) {
        // all positive charged hadrons
        int in_flag = 0;
        for (int i = 0; i < 5; i++) {
            if (abs(pid) == charged_hadron_urqmd_id_list[i]) {
                in_flag = 1;
                break;
            }
        }
        if (in_flag == 1 && charge > 0)
            pick_flag = 1;
    } else if (particle_urqmd_id == -9998) {
        // all negative charged hadrons
        int in_flag = 0;
        for (int i = 0; i < 5; i++) {
            if (abs(pid) == charged_hadron_urqmd_id_list[i]) {
                in_flag = 1;
                break;
            }
        }
        if (in_flag == 1 && charge < 0)
            pick_flag = 1;
    } else if (particle_urqmd_id == 9997) {
        // all baryons
        int in_flag = 0;
        for (int i = 0; i < 5; i++) {
            if (pid == baryon_urqmd_id_list[i]) {
                in_flag = 1;
                break;
            }
        }
        if (in_flag == 1)
            pick_flag = 1;
    } else if (particle_urqmd_id == -9997) {
        // all anti-baryons
        int in_flag = 0;
        for (int i = 0; i < 5; i++) {
            if (pid == -baryon_urqmd_id_list[i]) {
                in_flag = 1;
                break;
            }
        }
        if (in_flag == 1)
            pick_flag = 1;
    } else {
        // identified particle
        if (flag_isospin == 0) {
            if (pid == particle_urqmd_id) {
                if (reject_decay_flag == 1 && parent_proc_type == 20) {
                    pick_flag = 0;
                } else {
                    pick_flag = 1;
                }
            }
        } else {
            if (pid == particle_urqmd_id && iso3 == particle_urqmd_isospin) {
                if (reject_decay_flag == 1 && parent_proc_type == 20) {
                    pick_flag = 0;
                } else {
                    pick_flag = 1;
                }
            }
        }
    }
    return(pick_flag);
}

int particleSamples::decide_to_pick_UrQMD_anti_particles(int pid, int iso3,
                                                         int charge) {
    // this function judge whether particle is the anti-particle of the
    // particle of interest
    int pick_flag = 0;
    if (particle_urqmd_id == 9998) {
        // anti-particles for all positive charged hadrons
        // pick all negative charged hadrons
        int in_flag = 0;
        for (int i = 0; i < 5; i++) {
            if (abs(pid) == charged_hadron_urqmd_id_list[i]) {
                in_flag = 1;
                break;
            }
        }
        if (in_flag == 1 && charge < 0)
            pick_flag = 1;
    } else if (particle_urqmd_id == 9997) {
        // anti-particles for all baryons:
        // pick all anti-baryons
        int in_flag = 0;
        for (int i = 0; i < 5; i++) {
            if (pid == -baryon_urqmd_id_list[i]) {
                in_flag = 1;
                break;
            }
        }
        if (in_flag == 1)
            pick_flag = 1;
    } else {
        // for identified particle: pick its anti-particle
        if (pid == -particle_urqmd_id && iso3 == -particle_urqmd_isospin) {
            pick_flag = 1;
        }
    }
    return(pick_flag);
}

int particleSamples::decide_to_pick_UrQMD_resonance(int pid, int iso3,
                                                    int charge) {
    int pick_flag = 0;
    if (particle_urqmd_id == 27) {
        // Lambda particles, we will consider feed down from Sigma^0
        if (pid == 40 && iso3 == 0 && charge == 0) {
            pick_flag = 1;
        }
    } else if (particle_urqmd_id == -27) {
        // Anti-Lambda particles, we will consider feed down from Anti-Sigma^0
        if (pid == -40 && iso3 == 0 && charge == 0) {
            pick_flag = 1;
        }
    } else {
        // Other particles do not have weak decay contributions
        pick_flag = 0;
    }
    return(pick_flag);
}

void particleSamples::decide_to_pick_UrQMD_reconst(
                    int pid, int iso3, int charge, int parent_proc_type,
                    int *flag1, int *flag2) {
    *flag1 = 0;
    *flag2 = 0;
    if (particle_monval == 333) {
        // particle of interest is phi(1020)
        // collect (K^+, K^-) pairs for reconstruction
        if (pid == 106 && iso3 == 1
            && charge == 1 && parent_proc_type == 20) {
            // current particle is K^+ from a decay
            *flag1 = 1;
        }
        if (pid == -106 && iso3 == -1
            && charge == -1 && parent_proc_type == 20) {
            // current particle is K^- from a decay
            *flag2 = 1;
        }
    }
    return;
}

int particleSamples::decide_to_pick_JAM(int pid, int *charge_flag) {
    int pick_flag = 0;
    for (int i = 0; i < 6; i++) {
        if (abs(pid) == charged_hadron_pdg_list[i]) {
            // includes anti-particles
            if (particle_monval == 9999) {
                pick_flag = 1;
            }
            if (pid > 0) {
                if (i < 4) {
                    *charge_flag = 1;
                } else {
                    *charge_flag = -1;
                }
            } else {
                if (i < 4) {
                    *charge_flag = -1;
                } else {
                    *charge_flag = 1;
                }
            }
            if (particle_monval == 9998 && *charge_flag == 1) {
                pick_flag = 1;
            }
            if (particle_monval == -9998 && *charge_flag == -1) {
                pick_flag = 1;
            }
            break;
        }
    }
    return(pick_flag);
}

int particleSamples::decide_to_pick_from_OSCAR_file(int temp_monval) {
    int pick_flag = 0;
    if (resonance_feed_down_flag == 1) {
        if (select_resonances_flag == 0) {
            // pick everything
            pick_flag = 1;
        } else {
            for (unsigned int ireso = 0; ireso < select_resonances_list.size();
                    ireso++) {
                if (temp_monval == select_resonances_list[ireso]) {
                    pick_flag = 1;
                    break;
                }
            }
        }
    } else if (flag_isospin == 0) {
        // do not distinguish the isospin of the particle
        if (abs(temp_monval) == particle_monval) {
            pick_flag = 1;
        }
    } else {
        // distinguish the isospin of the particle
        if (temp_monval == particle_monval) {
            pick_flag = 1;
        }
        if (net_particle_flag == 1) {
            if (temp_monval == -particle_monval) {
                pick_flag = 2;  // set to 2 to trick anti_particle_pick_flag
            }
        }
    }
    return(pick_flag);
}

int particleSamples::decide_to_pick_OSCAR(int monval) {
    int pick_flag = 0;
    if (particle_monval == 9999) {
        int charge = decayer_ptr->get_particle_charge(monval);
        if (charge != 0) {
            pick_flag = 1;
        }
    } else if (particle_monval == 9998) {
        int charge = decayer_ptr->get_particle_charge(monval);
        if (charge > 0) {
            pick_flag = 1;
        }
    } else if (particle_monval == -9998) {
        int charge = decayer_ptr->get_particle_charge(monval);
        if (charge < 0) {
            pick_flag = 1;
        }
    } else if (particle_monval == 9997) {
        int baryon = decayer_ptr->get_particle_baryon_number(monval);
        if (baryon > 0) {
            pick_flag = 1;
        }
    } else if (particle_monval == -9997) {
        int baryon = decayer_ptr->get_particle_baryon_number(monval);
        if (baryon < 0) {
            pick_flag = 1;
        }
    } else if (particle_monval == 9996) {
        int strange = decayer_ptr->get_particle_strange_number(monval);
        if (strange > 0) {
            pick_flag = 1;
        }
    } else if (particle_monval == -9996) {
        int strange = decayer_ptr->get_particle_strange_number(monval);
        if (strange < 0) {
            pick_flag = 1;
        }
    } else {
        if (monval == particle_monval) {
            pick_flag = 1;
        }
    }
    return(pick_flag);
}

int particleSamples::read_in_particle_samples_OSCAR() {
    // clean out the previous record
    for (unsigned int i = 0; i < particle_list->size(); i++)
        (*particle_list)[i]->clear();
    particle_list->clear();
    if (net_particle_flag == 1) {
        for (unsigned int i = 0; i < anti_particle_list->size(); i++) {
            (*anti_particle_list)[i]->clear();
        }
        anti_particle_list->clear();
    }
    
    string temp_string;
    int event_id, n_particle, dummy;
    int ievent;
    int temp_monval;
    for (ievent = 0; ievent < event_buffer_size; ievent++) {
        getline(inputfile, temp_string);
        stringstream temp1(temp_string);
        temp1 >> event_id >> n_particle;
        if (!inputfile.eof()) {
            particle_list->push_back(new vector<particle_info> );
            if (net_particle_flag == 1) {
                anti_particle_list->push_back(new vector<particle_info> );
            }
            int idx = ievent;

            int pick_flag = 0;
            for (int ipart = 0; ipart < n_particle; ipart++) {
                getline(inputfile, temp_string);
                stringstream temp2(temp_string);
                temp2 >> dummy >> temp_monval;
                pick_flag = decide_to_pick_from_OSCAR_file(temp_monval);
                if (pick_flag != 0) {
                     particle_info *temp_particle_info = new particle_info;
                     temp_particle_info->monval = temp_monval;
                     temp2 >> temp_particle_info->px 
                           >> temp_particle_info->py
                           >> temp_particle_info->pz 
                           >> temp_particle_info->E
                           >> temp_particle_info->mass 
                           >> temp_particle_info->x 
                           >> temp_particle_info->y
                           >> temp_particle_info->z 
                           >> temp_particle_info->t;
                     if (pick_flag == 1) {
                        (*particle_list)[idx]->push_back(*temp_particle_info);
                     } else if (pick_flag == 2) {
                        (*anti_particle_list)[idx]->push_back(
                                                        *temp_particle_info);
                     }
                     delete temp_particle_info;
                }
            }
        } else {
            break;
        }
    }
    return(0);
}

int particleSamples::read_in_particle_samples_JAM() {
    // clean out the previous record
    for (unsigned int i = 0; i < particle_list->size(); i++)
        (*particle_list)[i]->clear();
    particle_list->clear();
    
    if (resonance_weak_feed_down_flag == 1) {
        for (unsigned int i = 0; i < resonance_list->size(); i++)
            (*resonance_list)[i]->clear();
        resonance_list->clear();
    }
    
    if (reconst_flag == 1) {
        for (unsigned int i = 0; i < reconst_list_1->size(); i++)
            (*reconst_list_1)[i]->clear();
        reconst_list_1->clear();
        for (unsigned int i = 0; i < reconst_list_2->size(); i++)
            (*reconst_list_2)[i]->clear();
        reconst_list_2->clear();
    }
    
    
    if (flag_charge_dependence == 1) {
        for (unsigned int i = 0; i < positive_charge_hadron_list->size(); i++)
            (*positive_charge_hadron_list)[i]->clear();
        positive_charge_hadron_list->clear();
        for (unsigned int i = 0; i < negative_charge_hadron_list->size(); i++)
            (*negative_charge_hadron_list)[i]->clear();
        negative_charge_hadron_list->clear();
    }
    
    string temp_string;
    int event_id, n_particle;
    char cdummy;
    int ievent;
    int temp_monval;
    for (ievent = 0; ievent < event_buffer_size; ievent++) {
        getline(inputfile, temp_string);
        stringstream temp1(temp_string);
        temp1 >> cdummy >> event_id >> n_particle;
        if (!inputfile.eof()) {
            particle_list->push_back(new vector<particle_info> );

            if (resonance_weak_feed_down_flag == 1) {
                resonance_list->push_back(new vector<particle_info>);
            }

            if (reconst_flag == 1) {
                reconst_list_1->push_back(new vector<particle_info>);
                reconst_list_2->push_back(new vector<particle_info>);
            }
            
            if (flag_charge_dependence == 1) {
                positive_charge_hadron_list->push_back(
                                                new vector<particle_info>);
                negative_charge_hadron_list->push_back(
                                                new vector<particle_info>);
            }

            int idx = ievent;
            (*particle_list)[idx]->clear(); // clean out the previous record
            if (resonance_weak_feed_down_flag == 1) {
                (*resonance_list)[idx]->clear();
            }
            if (reconst_flag == 1) {
                (*reconst_list_1)[idx]->clear();
                (*reconst_list_2)[idx]->clear();
            }
            if (flag_charge_dependence == 1) {
                (*positive_charge_hadron_list)[idx]->clear();
                (*negative_charge_hadron_list)[idx]->clear();
            }

            int pick_flag = 0;
            int charge_flag = 0;
            int resonance_pick_flag = 0;
            int reconst_1_pick_flag = 0;
            int reconst_2_pick_flag = 0;
            for (int ipart = 0; ipart < n_particle; ipart++) {
                getline(inputfile, temp_string);
                stringstream temp2(temp_string);
                temp2 >> temp_monval;

                // particle of interest
                if (flag_isospin == 0) {
                    if (abs(temp_monval) == particle_monval)
                        pick_flag = 1;
                    else
                        pick_flag = 0;
                } else {
                    if (temp_monval == particle_monval)
                        pick_flag = 1;
                    else
                        pick_flag = 0;
                }

                if (fabs(particle_monval) > 9990) {
                    pick_flag = decide_to_pick_JAM(temp_monval, &charge_flag);
                }

                if (resonance_weak_feed_down_flag == 1) {
                    if (temp_monval == 3122 && particle_monval == 3212) {
                        // pick Sigma^0 for Lambda
                        resonance_pick_flag = 1;
                    } else if (temp_monval == -3122
                               && particle_monval == -3212) {
                        // pick anti-Sigma^0 for anti-Lambda
                        resonance_pick_flag = 1;
                    }
                }

                if (reconst_flag == 1) {
                    if (particle_monval == 333) {
                        // particle of interest is phi(1020)
                        if (temp_monval == 321) {
                            // current particle is K^+
                            reconst_1_pick_flag = 1;
                        } else if (temp_monval == -321) {
                            // current particle is K^-
                            reconst_2_pick_flag = 1;
                        }
                    }
                }

                int trigger = (pick_flag + resonance_pick_flag
                               + reconst_1_pick_flag + reconst_2_pick_flag);
                if (trigger > 0) {
                    // pick up
                    particle_info *temp_particle_info = new particle_info;
                    temp2 >> temp_particle_info->mass
                          >> temp_particle_info->px
                          >> temp_particle_info->py
                          >> temp_particle_info->pz
                          >> temp_particle_info->x 
                          >> temp_particle_info->y
                          >> temp_particle_info->z 
                          >> temp_particle_info->t;
                    temp_particle_info->E = sqrt(
                        temp_particle_info->mass*temp_particle_info->mass
                        + temp_particle_info->px*temp_particle_info->px
                        + temp_particle_info->py*temp_particle_info->py
                        + temp_particle_info->pz*temp_particle_info->pz);
                    if (pick_flag == 1) {
                        (*particle_list)[idx]->push_back(*temp_particle_info);
                        if (flag_charge_dependence == 1) {
                            if (charge_flag > 0) {
                                (*positive_charge_hadron_list)[idx]->push_back(
                                                    *temp_particle_info);
                            } else {
                                (*negative_charge_hadron_list)[idx]->push_back(
                                                    *temp_particle_info);
                            }
                        }
                    } else if (resonance_pick_flag == 1) {
                        temp_particle_info->monval = temp_monval;
                        (*resonance_list)[idx]->push_back(*temp_particle_info);
                    } else if (reconst_1_pick_flag == 1) {
                        (*reconst_list_1)[idx]->push_back(*temp_particle_info);
                    } else if (reconst_2_pick_flag == 1) {
                        (*reconst_list_2)[idx]->push_back(*temp_particle_info);
                    }
                    delete temp_particle_info;
                }
            }
        } else {
            break;
        }
    }
    return(0);
}

int particleSamples::read_in_particle_samples_OSCAR_mixed_event() {
    // clean out the previous record
    for (unsigned int i = 0; i < particle_list_mixed_event->size(); i++)
        (*particle_list_mixed_event)[i]->clear();
    particle_list_mixed_event->clear();
    
    string temp_string;
    int event_id, n_particle, dummy;
    int ievent;
    int temp_monval;
    for (ievent = 0; ievent < event_buffer_size; ievent++) {
        getline(inputfile_mixed_event, temp_string);
        stringstream temp1(temp_string);
        temp1 >> event_id >> n_particle;
        if (!inputfile_mixed_event.eof()) {
            particle_list_mixed_event->push_back(new vector<particle_info> );
            int idx = ievent;

            int pick_flag = 0;
            for (int ipart = 0; ipart < n_particle; ipart++) {
                getline(inputfile_mixed_event, temp_string);
                stringstream temp2(temp_string);
                temp2 >> dummy >> temp_monval;
                pick_flag = decide_to_pick_from_OSCAR_file(temp_monval);
                if (pick_flag == 1) {
                     particle_info *temp_particle_info = new particle_info;
                     temp_particle_info->monval = temp_monval;
                     temp2 >> temp_particle_info->px 
                           >> temp_particle_info->py
                           >> temp_particle_info->pz 
                           >> temp_particle_info->E
                           >> temp_particle_info->mass 
                           >> temp_particle_info->x 
                           >> temp_particle_info->y
                           >> temp_particle_info->z 
                           >> temp_particle_info->t;
                     (*particle_list_mixed_event)[idx]->push_back(
                                                        *temp_particle_info);
                     delete temp_particle_info;
                }
            }
        } else {
            break;
        }
    }
    return(0);
}

int particleSamples::read_in_particle_samples_JAM_mixed_event() {
    // clean out the previous record
    for (unsigned int i = 0; i < particle_list_mixed_event->size(); i++)
        (*particle_list_mixed_event)[i]->clear();
    particle_list_mixed_event->clear();
    
    if (resonance_weak_feed_down_flag == 1) {
        for (unsigned int i = 0; i < resonance_list->size(); i++)
            (*resonance_list)[i]->clear();
        resonance_list->clear();
    }
    
    string temp_string;
    int event_id, n_particle;
    char cdummy;
    int ievent;
    int temp_monval;
    for (ievent = 0; ievent < event_buffer_size; ievent++) {
        getline(inputfile_mixed_event, temp_string);
        stringstream temp1(temp_string);
        temp1 >> cdummy >> event_id >> n_particle;
        if (!inputfile_mixed_event.eof()) {
            particle_list_mixed_event->push_back(new vector<particle_info> );
            if (resonance_weak_feed_down_flag == 1) {
                resonance_list->push_back(new vector<particle_info>);
            }
            int idx = ievent;

            int pick_flag = 0;
            int charge_flag = 0;
            int resonance_pick_flag = 0;
            for (int ipart = 0; ipart < n_particle; ipart++) {
                getline(inputfile_mixed_event, temp_string);
                stringstream temp2(temp_string);
                temp2 >> temp_monval;
                if (flag_isospin == 0) {
                    if (abs(temp_monval) == particle_monval)
                        pick_flag = 1;
                    else
                        pick_flag = 0;
                } else {
                    if (temp_monval == particle_monval)
                        pick_flag = 1;
                    else
                        pick_flag = 0;
                }
                if (particle_monval == 9999) {
                    pick_flag = decide_to_pick_JAM(temp_monval, &charge_flag);
                }

                if (resonance_weak_feed_down_flag == 1) {
                    // pick up Sigma^0 for Lambda
                    if (temp_monval == 3122 && particle_monval == 3212) {
                        resonance_pick_flag = 1;
                    } else if (temp_monval == -3122
                               && particle_monval == -3212) {
                        resonance_pick_flag = 1;
                    }
                }

                int trigger = pick_flag + resonance_pick_flag;
                if (trigger > 0) {
                     particle_info *temp_particle_info = new particle_info;
                     temp2 >> temp_particle_info->mass
                           >> temp_particle_info->px
                           >> temp_particle_info->py
                           >> temp_particle_info->pz
                           >> temp_particle_info->x
                           >> temp_particle_info->y
                           >> temp_particle_info->z
                           >> temp_particle_info->t;
                     temp_particle_info->E = sqrt(
                         temp_particle_info->mass*temp_particle_info->mass
                         + temp_particle_info->px*temp_particle_info->px
                         + temp_particle_info->py*temp_particle_info->py
                         + temp_particle_info->pz*temp_particle_info->pz);
                     if (pick_flag == 1) {
                        (*particle_list_mixed_event)[idx]->push_back(
                                                        *temp_particle_info);
                     } else if (resonance_pick_flag == 1) {
                         temp_particle_info->monval = temp_monval;
                        (*resonance_list)[idx]->push_back(*temp_particle_info);
                     }
                     delete temp_particle_info;
                }
            }
        } else {
            break;
        }
    }
    return(0);
}

int particleSamples::read_in_particle_samples_UrQMD() {
    // clean out the previous record
    for (unsigned int i = 0; i < particle_list->size(); i++)
        (*particle_list)[i]->clear();
    particle_list->clear();
    
    if (resonance_weak_feed_down_flag == 1) {
        for (unsigned int i = 0; i < resonance_list->size(); i++)
            (*resonance_list)[i]->clear();
        resonance_list->clear();
    }

    if (net_particle_flag == 1) {
        for (unsigned int i = 0; i < anti_particle_list->size(); i++) {
            (*anti_particle_list)[i]->clear();
        }
        anti_particle_list->clear();
    }

    if (reconst_flag == 1) {
        for (unsigned int i = 0; i < reconst_list_1->size(); i++)
            (*reconst_list_1)[i]->clear();
        reconst_list_1->clear();
        for (unsigned int i = 0; i < reconst_list_2->size(); i++)
            (*reconst_list_2)[i]->clear();
        reconst_list_2->clear();
    }

    string temp_string;
    int n_particle;
    double dummy;
    int parent_proc_type;
    int ievent;
    int urqmd_pid, urqmd_iso3, urqmd_charge;
    double temp_mass;
    for (ievent = 0; ievent < event_buffer_size; ievent++) {
        getline(inputfile, temp_string);
        if (!inputfile.eof()) {
            // create one event
            particle_list->push_back(new vector<particle_info> );
            if (resonance_weak_feed_down_flag == 1) {
                resonance_list->push_back(new vector<particle_info>);
            }
            if (reconst_flag == 1) {
                reconst_list_1->push_back(new vector<particle_info>);
                reconst_list_2->push_back(new vector<particle_info>);
            }
            if (net_particle_flag == 1) {
                anti_particle_list->push_back(new vector<particle_info>);
            }

            // first skip the header
            for (int i = 0; i < 16; i++)
                getline(inputfile, temp_string);
            // then get number of particles within the event
            getline(inputfile, temp_string);
            stringstream temp1(temp_string);
            temp1 >> n_particle;
            getline(inputfile, temp_string);  // then get one useless line

            int idx = ievent;
            (*particle_list)[idx]->clear(); // clean out the previous record
            if (resonance_weak_feed_down_flag == 1) {
                (*resonance_list)[idx]->clear();
            }
            if (net_particle_flag == 1) {
                (*anti_particle_list)[idx]->clear();
            }
            if (reconst_flag == 1) {
                (*reconst_list_1)[idx]->clear();
                (*reconst_list_2)[idx]->clear();
            }

            int pick_flag = 0;
            int resonance_pick_flag = 0;
            int reconst_1_pick_flag = 0;
            int reconst_2_pick_flag = 0;
            int anti_particle_pick_flag = 0;
            for (int ipart = 0; ipart < n_particle; ipart++) {
                getline(inputfile, temp_string);
                stringstream temp2(temp_string);
                temp2 >> dummy >> dummy >> dummy >> dummy
                      >> dummy >> dummy >> dummy >> dummy
                      >> temp_mass >> urqmd_pid >> urqmd_iso3 >> urqmd_charge
                      >> dummy >> dummy >> parent_proc_type;

                pick_flag = decide_to_pick_UrQMD(
                        urqmd_pid, urqmd_iso3, urqmd_charge, parent_proc_type);

                if (resonance_weak_feed_down_flag == 1) {
                    resonance_pick_flag = decide_to_pick_UrQMD_resonance(
                                        urqmd_pid, urqmd_iso3, urqmd_charge);
                }

                if (net_particle_flag == 1) {
                    anti_particle_pick_flag =
                        decide_to_pick_UrQMD_anti_particles(
                                        urqmd_pid, urqmd_iso3, urqmd_charge);
                }

                if (reconst_flag == 1) {
                    decide_to_pick_UrQMD_reconst(
                        urqmd_pid, urqmd_iso3, urqmd_charge, parent_proc_type,
                        &reconst_1_pick_flag, &reconst_2_pick_flag);
                }

                int trigger = (pick_flag + resonance_pick_flag
                               + net_particle_flag
                               + reconst_1_pick_flag + reconst_2_pick_flag);
                if (trigger > 0) {
                    particle_info *temp_particle_info = new particle_info;
                    temp2 >> temp_particle_info->t
                          >> temp_particle_info->x 
                          >> temp_particle_info->y
                          >> temp_particle_info->z 
                          >> temp_particle_info->E
                          >> temp_particle_info->px 
                          >> temp_particle_info->py
                          >> temp_particle_info->pz;
                    temp_particle_info->mass = temp_mass;
                    if (pick_flag == 1) {
                        if (reject_decay_flag == 2 && parent_proc_type == 20) {
                            double tau = sqrt(
                                temp_particle_info->t*temp_particle_info->t
                                - temp_particle_info->z*temp_particle_info->z);
                            if (tau > tau_reject) {
                                pick_flag = 0;
                            } else {
                                pick_flag = 1;
                            }
                        }
                        if (pick_flag == 1) {
                            (*particle_list)[idx]->push_back(
                                                        *temp_particle_info);
                        }
                    } else if (resonance_pick_flag == 1) {
                        temp_particle_info->monval = get_pdg_id(urqmd_pid,
                                                                urqmd_iso3);
                        (*resonance_list)[idx]->push_back(*temp_particle_info);
                    } else if (reconst_1_pick_flag == 1) {
                        (*reconst_list_1)[idx]->push_back(*temp_particle_info);
                    } else if (reconst_2_pick_flag == 1) {
                        (*reconst_list_2)[idx]->push_back(*temp_particle_info);
                    } else if (anti_particle_pick_flag == 1) {
                        (*anti_particle_list)[idx]->push_back(
                                                        *temp_particle_info);
                    }
                    delete temp_particle_info;
                }
            }
        } else {
            break;
        }
    }
    return(0);
}

int particleSamples::read_in_particle_samples_UrQMD_zipped() {
    // clean out the previous record
    for (unsigned int i = 0; i < particle_list->size(); i++)
        (*particle_list)[i]->clear();
    particle_list->clear();
    
    if (resonance_weak_feed_down_flag == 1) {
        for (unsigned int i = 0; i < resonance_list->size(); i++)
            (*resonance_list)[i]->clear();
        resonance_list->clear();
    }

    if (net_particle_flag == 1) {
        for (unsigned int i = 0; i < anti_particle_list->size(); i++) {
            (*anti_particle_list)[i]->clear();
        }
        anti_particle_list->clear();
    }

    if (reconst_flag == 1) {
        for (unsigned int i = 0; i < reconst_list_1->size(); i++)
            (*reconst_list_1)[i]->clear();
        reconst_list_1->clear();
        for (unsigned int i = 0; i < reconst_list_2->size(); i++)
            (*reconst_list_2)[i]->clear();
        reconst_list_2->clear();
    }

    if (flag_charge_dependence == 1) {
        for (unsigned int i = 0; i < positive_charge_hadron_list->size(); i++)
            (*positive_charge_hadron_list)[i]->clear();
        positive_charge_hadron_list->clear();
        for (unsigned int i = 0; i < negative_charge_hadron_list->size(); i++)
            (*negative_charge_hadron_list)[i]->clear();
        negative_charge_hadron_list->clear();
    }

    string temp_string;
    int n_particle;
    double dummy;
    int parent_proc_type;
    int ievent;
    int urqmd_pid, urqmd_iso3, urqmd_charge;
    for (ievent = 0; ievent < event_buffer_size; ievent++) {
        temp_string = gz_readline(inputfile_gz);
        if (!gzeof(inputfile_gz)) {
            // create one event
            particle_list->push_back(new vector<particle_info> );
            if (resonance_weak_feed_down_flag == 1) {
                resonance_list->push_back(new vector<particle_info>);
            }
            if (reconst_flag == 1) {
                reconst_list_1->push_back(new vector<particle_info>);
                reconst_list_2->push_back(new vector<particle_info>);
            }
            if (net_particle_flag == 1) {
                anti_particle_list->push_back(new vector<particle_info>);
            }
            if (flag_charge_dependence == 1) {
                positive_charge_hadron_list->push_back(
                                                new vector<particle_info>);
                negative_charge_hadron_list->push_back(
                                                new vector<particle_info>);
            }

            // get number of particles within the event
            stringstream temp1(temp_string);
            temp1 >> n_particle;
            // cout << "first line: " << temp_string << endl;
            // cout << "check n_particle = " << n_particle << endl;
            temp_string = gz_readline(inputfile_gz);  // then get one useless line

            int idx = ievent;
            (*particle_list)[idx]->clear(); // clean out the previous record
            if (resonance_weak_feed_down_flag == 1) {
                (*resonance_list)[idx]->clear();
            }
            if (net_particle_flag == 1) {
                (*anti_particle_list)[idx]->clear();
            }
            if (reconst_flag == 1) {
                (*reconst_list_1)[idx]->clear();
                (*reconst_list_2)[idx]->clear();
            }
            if (flag_charge_dependence == 1) {
                (*positive_charge_hadron_list)[idx]->clear();
                (*negative_charge_hadron_list)[idx]->clear();
            }

            int pick_flag = 0;
            int resonance_pick_flag = 0;
            int reconst_1_pick_flag = 0;
            int reconst_2_pick_flag = 0;
            int anti_particle_pick_flag = 0;
            for (int ipart = 0; ipart < n_particle; ipart++) {
                temp_string = gz_readline(inputfile_gz);
                stringstream temp2(temp_string);
                temp2 >> urqmd_pid >> urqmd_iso3 >> urqmd_charge
                      >> dummy >> dummy >> parent_proc_type;

                pick_flag = decide_to_pick_UrQMD(
                        urqmd_pid, urqmd_iso3, urqmd_charge, parent_proc_type);

                if (resonance_weak_feed_down_flag == 1) {
                    resonance_pick_flag = decide_to_pick_UrQMD_resonance(
                                        urqmd_pid, urqmd_iso3, urqmd_charge);
                }

                if (net_particle_flag == 1) {
                    anti_particle_pick_flag =
                        decide_to_pick_UrQMD_anti_particles(
                                        urqmd_pid, urqmd_iso3, urqmd_charge);
                }

                if (reconst_flag == 1) {
                    decide_to_pick_UrQMD_reconst(
                        urqmd_pid, urqmd_iso3, urqmd_charge, parent_proc_type,
                        &reconst_1_pick_flag, &reconst_2_pick_flag);
                }

                int trigger = (pick_flag + resonance_pick_flag
                               + net_particle_flag
                               + reconst_1_pick_flag + reconst_2_pick_flag);
                if (trigger > 0) {
                    particle_info *temp_particle_info = new particle_info;
                    temp2 >> temp_particle_info->mass
                          >> temp_particle_info->t
                          >> temp_particle_info->x 
                          >> temp_particle_info->y
                          >> temp_particle_info->z 
                          >> temp_particle_info->E
                          >> temp_particle_info->px 
                          >> temp_particle_info->py
                          >> temp_particle_info->pz;
                    temp_particle_info->monval = get_pdg_id(urqmd_pid,
                                                            urqmd_iso3);

                    if (pick_flag == 1) {
                        if (reject_decay_flag == 2 && parent_proc_type == 20) {
                            double tau = sqrt(
                                temp_particle_info->t*temp_particle_info->t
                                - temp_particle_info->z*temp_particle_info->z);
                            if (tau > tau_reject) {
                                pick_flag = 0;
                            } else {
                                pick_flag = 1;
                            }
                        }
                        if (pick_flag == 1) {
                            (*particle_list)[idx]->push_back(
                                                        *temp_particle_info);
                            if (flag_charge_dependence == 1) {
                                if (urqmd_charge > 0) {
                                    (*positive_charge_hadron_list)[idx]->push_back(
                                                        *temp_particle_info);
                                } else {
                                    (*negative_charge_hadron_list)[idx]->push_back(
                                                        *temp_particle_info);
                                }
                            }
                        }
                    } else if (resonance_pick_flag == 1) {
                        (*resonance_list)[idx]->push_back(*temp_particle_info);
                    } else if (reconst_1_pick_flag == 1) {
                        (*reconst_list_1)[idx]->push_back(*temp_particle_info);
                    } else if (reconst_2_pick_flag == 1) {
                        (*reconst_list_2)[idx]->push_back(*temp_particle_info);
                    } else if (anti_particle_pick_flag == 1) {
                        (*anti_particle_list)[idx]->push_back(
                                                        *temp_particle_info);
                    }
                    delete temp_particle_info;
                }
            }
        } else {
            break;
        }
    }
    return(0);
}


int particleSamples::read_in_particle_samples_gzipped() {
    // clean out the previous record
    for (unsigned int i = 0; i < particle_list->size(); i++) {
        (*particle_list)[i]->clear();
    }
    particle_list->clear();
    
    string temp_string;
    int n_particle;
    int temp_monval;
    for (int ievent = 0; ievent < event_buffer_size; ievent++) {
        temp_string = gz_readline(inputfile_gz);
        if (!gzeof(inputfile_gz)) {
            // create one event
            particle_list->push_back(new vector<particle_info> );

            // get number of particles within the event
            stringstream temp1(temp_string);
            temp1 >> n_particle;
            //cout << "first line: " << temp_string << endl;
            //cout << "check n_particle = " << n_particle << endl;

            int idx = ievent;
            (*particle_list)[idx]->clear(); // clean out the previous record

            int pick_flag = 0;
            for (int ipart = 0; ipart < n_particle; ipart++) {
                temp_string = gz_readline(inputfile_gz);
                stringstream temp2(temp_string);
                temp2 >> temp_monval;
                pick_flag = decide_to_pick_from_OSCAR_file(temp_monval);
                if (pick_flag != 0) {
                    particle_info *temp_particle_info = new particle_info;
                    temp_particle_info->monval = temp_monval;
                    temp2 >> temp_particle_info->mass
                          >> temp_particle_info->t
                          >> temp_particle_info->x 
                          >> temp_particle_info->y
                          >> temp_particle_info->z 
                          >> temp_particle_info->E
                          >> temp_particle_info->px 
                          >> temp_particle_info->py
                          >> temp_particle_info->pz;

                    (*particle_list)[idx]->push_back(*temp_particle_info);
                    delete temp_particle_info;
                }
            }
        } else {
            break;
        }
    }
    return(0);
}

int particleSamples::read_in_particle_samples_UrQMD_3p3() {
    // clean out the previous record
    for (unsigned int i = 0; i < particle_list->size(); i++)
        (*particle_list)[i]->clear();
    particle_list->clear();
    
    if (resonance_weak_feed_down_flag == 1) {
        for (unsigned int i = 0; i < resonance_list->size(); i++)
            (*resonance_list)[i]->clear();
        resonance_list->clear();
    }
    
    if (reconst_flag == 1) {
        for (unsigned int i = 0; i < reconst_list_1->size(); i++)
            (*reconst_list_1)[i]->clear();
        reconst_list_1->clear();
        for (unsigned int i = 0; i < reconst_list_2->size(); i++)
            (*reconst_list_2)[i]->clear();
        reconst_list_2->clear();
    }

    string temp_string;
    int n_particle;
    double dummy;
    int parent_proc_type;
    int ievent;
    int urqmd_pid, urqmd_iso3, urqmd_charge;
    double temp_mass;
    for (ievent = 0; ievent < event_buffer_size; ievent++) {
        getline(inputfile, temp_string);
        if (!inputfile.eof()) {
            // create one event
            particle_list->push_back(new vector<particle_info> );
            if (resonance_weak_feed_down_flag == 1) {
                resonance_list->push_back(new vector<particle_info>);
            }
            if (reconst_flag == 1) {
                reconst_list_1->push_back(new vector<particle_info>);
                reconst_list_2->push_back(new vector<particle_info>);
            }

            // first skip the header
            for (int i = 0; i < 13; i++)
                getline(inputfile, temp_string);
            // then get number of particles within the event
            getline(inputfile, temp_string);
            stringstream temp1(temp_string);
            temp1 >> n_particle;
            getline(inputfile, temp_string);  // then get one useless line

            int idx = ievent;
            (*particle_list)[idx]->clear(); // clean out the previous record
            if (resonance_weak_feed_down_flag == 1) {
                (*resonance_list)[idx]->clear();
            }
            if (reconst_flag == 1) {
                (*reconst_list_1)[idx]->clear();
                (*reconst_list_2)[idx]->clear();
            }

            int pick_flag = 0;
            int resonance_pick_flag = 0;
            int reconst_1_pick_flag = 0;
            int reconst_2_pick_flag = 0;
            for (int ipart = 0; ipart < n_particle; ipart++) {
                getline(inputfile, temp_string);
                stringstream temp2(temp_string);
                temp2 >> dummy >> dummy >> dummy >> dummy
                      >> dummy >> dummy >> dummy >> dummy
                      >> temp_mass >> urqmd_pid >> urqmd_iso3 >> urqmd_charge
                      >> dummy >> dummy >> parent_proc_type;

                pick_flag = decide_to_pick_UrQMD(
                        urqmd_pid, urqmd_iso3, urqmd_charge, parent_proc_type);

                if (resonance_weak_feed_down_flag == 1) {
                    resonance_pick_flag = decide_to_pick_UrQMD_resonance(
                                        urqmd_pid, urqmd_iso3, urqmd_charge);
                }

                if (reconst_flag == 1) {
                    decide_to_pick_UrQMD_reconst(
                        urqmd_pid, urqmd_iso3, urqmd_charge, parent_proc_type,
                        &reconst_1_pick_flag, &reconst_2_pick_flag);
                }
                
                int trigger = (pick_flag + resonance_pick_flag
                               + reconst_1_pick_flag + reconst_2_pick_flag);
                if (trigger > 0) {
                    particle_info *temp_particle_info = new particle_info;
                    temp2 >> temp_particle_info->t
                          >> temp_particle_info->x 
                          >> temp_particle_info->y
                          >> temp_particle_info->z 
                          >> temp_particle_info->E
                          >> temp_particle_info->px 
                          >> temp_particle_info->py
                          >> temp_particle_info->pz;
                    temp_particle_info->mass = temp_mass;
                    temp_particle_info->monval = get_pdg_id(urqmd_pid,
                                                            urqmd_iso3);
                    if (pick_flag == 1) {
                        if (reject_decay_flag == 2 && parent_proc_type == 20) {
                            double tau = sqrt(
                                temp_particle_info->t*temp_particle_info->t
                                - temp_particle_info->z*temp_particle_info->z);
                            if (tau > tau_reject) {
                                pick_flag = 0;
                            } else {
                                pick_flag = 1;
                            }
                        }
                        if (pick_flag == 1) {
                            (*particle_list)[idx]->push_back(
                                                        *temp_particle_info);
                        }
                    } else if (resonance_pick_flag == 1) {
                        (*resonance_list)[idx]->push_back(*temp_particle_info);
                    } else if (reconst_1_pick_flag == 1) {
                        (*reconst_list_1)[idx]->push_back(*temp_particle_info);
                    } else if (reconst_2_pick_flag == 1) {
                        (*reconst_list_2)[idx]->push_back(*temp_particle_info);
                    }
                    delete temp_particle_info;
                }
            }
        } else {
            break;
        }
    }
    return(0);
}

int particleSamples::read_in_particle_samples_UrQMD_mixed_event() {
    // clean out the previous record
    for (unsigned int i = 0; i < particle_list_mixed_event->size(); i++)
        (*particle_list_mixed_event)[i]->clear();
    particle_list_mixed_event->clear();
    
    if (resonance_weak_feed_down_flag == 1) {
        for (unsigned int i = 0; i < resonance_list->size(); i++)
            (*resonance_list)[i]->clear();
        resonance_list->clear();
    }

    string temp_string;
    int n_particle;
    double dummy;
    int parent_proc_type;
    int ievent;
    int urqmd_pid, urqmd_iso3, urqmd_charge;
    double temp_mass;
    for (ievent = 0; ievent < event_buffer_size; ievent++) {
        getline(inputfile_mixed_event, temp_string);
        if (!inputfile_mixed_event.eof()) {
            particle_list_mixed_event->push_back(new vector<particle_info> );
            if (resonance_weak_feed_down_flag == 1) {
                resonance_list->push_back(new vector<particle_info>);
            }

            // first skip the header
            for(int i = 0; i < 16; i++)
                getline(inputfile_mixed_event, temp_string);
            // then get number of particles within the event
            getline(inputfile_mixed_event, temp_string);
            stringstream temp1(temp_string);
            temp1 >> n_particle;
            // then get one useless line
            getline(inputfile_mixed_event, temp_string);  

            // clean out the previous record
            int idx = ievent;
            (*particle_list_mixed_event)[idx]->clear(); 
            if (resonance_weak_feed_down_flag == 1) {
                (*resonance_list)[idx]->clear();
            }

            int pick_flag = 0;
            int resonance_pick_flag = 0;
            for(int ipart = 0; ipart < n_particle; ipart++)
            {
                getline(inputfile_mixed_event, temp_string);
                stringstream temp2(temp_string);
                temp2 >> dummy >> dummy >> dummy >> dummy
                      >> dummy >> dummy >> dummy >> dummy
                      >> temp_mass >> urqmd_pid >> urqmd_iso3 >> urqmd_charge
                      >> dummy >> dummy >> parent_proc_type;

                pick_flag = decide_to_pick_UrQMD(
                        urqmd_pid, urqmd_iso3, urqmd_charge, parent_proc_type);

                if (resonance_weak_feed_down_flag == 1) {
                    resonance_pick_flag = decide_to_pick_UrQMD_resonance(
                                        urqmd_pid, urqmd_iso3, urqmd_charge);
                }
                
                int trigger = pick_flag + resonance_pick_flag;
                if (trigger > 0) {
                    particle_info *temp_particle_info = new particle_info;
                    temp2 >> temp_particle_info->t
                          >> temp_particle_info->x 
                          >> temp_particle_info->y
                          >> temp_particle_info->z 
                          >> temp_particle_info->E
                          >> temp_particle_info->px 
                          >> temp_particle_info->py
                          >> temp_particle_info->pz ;
                    temp_particle_info->mass = temp_mass;
                    temp_particle_info->monval = get_pdg_id(urqmd_pid,
                                                            urqmd_iso3);
                    if (pick_flag == 1) {
                        if (reject_decay_flag == 2 && parent_proc_type == 20) {
                            double tau = sqrt(
                                temp_particle_info->t*temp_particle_info->t
                                - temp_particle_info->z*temp_particle_info->z);
                            if (tau > tau_reject) {
                                pick_flag = 0;
                            } else {
                                pick_flag = 1;
                            }
                        }
                        if (pick_flag == 1) {
                            (*particle_list_mixed_event)[idx]->push_back(
                                                           *temp_particle_info);
                        }
                    } else if (resonance_pick_flag == 1) {
                        (*resonance_list)[idx]->push_back(*temp_particle_info);
                    }
                    delete temp_particle_info;
                }
            }
        } else {
            break;
        }
    }
    return(0);
}

int particleSamples::read_in_particle_samples_UrQMD_mixed_event_zipped() {
    // clean out the previous record
    for (unsigned int i = 0; i < particle_list_mixed_event->size(); i++)
        (*particle_list_mixed_event)[i]->clear();
    particle_list_mixed_event->clear();
    
    if (resonance_weak_feed_down_flag == 1) {
        for (unsigned int i = 0; i < resonance_list->size(); i++)
            (*resonance_list)[i]->clear();
        resonance_list->clear();
    }

    string temp_string;
    int n_particle;
    double dummy;
    int parent_proc_type;
    int urqmd_pid, urqmd_iso3, urqmd_charge;
    for (int ievent = 0; ievent < event_buffer_size; ievent++) {
        temp_string = gz_readline(inputfile_mixed_event_gz);
        if (!gzeof(inputfile_mixed_event_gz)) {
            particle_list_mixed_event->push_back(new vector<particle_info> );
            if (resonance_weak_feed_down_flag == 1) {
                resonance_list->push_back(new vector<particle_info>);
            }

            // get number of particles within the event
            stringstream temp1(temp_string);
            temp1 >> n_particle;
            // then get one useless line
            temp_string = gz_readline(inputfile_mixed_event_gz);  

            // clean out the previous record
            int idx = ievent;
            (*particle_list_mixed_event)[idx]->clear(); 
            if (resonance_weak_feed_down_flag == 1) {
                (*resonance_list)[idx]->clear();
            }

            int pick_flag = 0;
            int resonance_pick_flag = 0;
            for (int ipart = 0; ipart < n_particle; ipart++) {
                temp_string = gz_readline(inputfile_mixed_event_gz);
                stringstream temp2(temp_string);
                temp2 >> urqmd_pid >> urqmd_iso3 >> urqmd_charge
                      >> dummy >> dummy >> parent_proc_type;

                pick_flag = decide_to_pick_UrQMD(
                        urqmd_pid, urqmd_iso3, urqmd_charge, parent_proc_type);

                if (resonance_weak_feed_down_flag == 1) {
                    resonance_pick_flag = decide_to_pick_UrQMD_resonance(
                                        urqmd_pid, urqmd_iso3, urqmd_charge);
                }
                
                int trigger = pick_flag + resonance_pick_flag;
                if (trigger > 0) {
                    particle_info *temp_particle_info = new particle_info;
                    temp2 >> temp_particle_info->mass
                          >> temp_particle_info->t
                          >> temp_particle_info->x
                          >> temp_particle_info->y
                          >> temp_particle_info->z
                          >> temp_particle_info->E
                          >> temp_particle_info->px
                          >> temp_particle_info->py
                          >> temp_particle_info->pz;
                    temp_particle_info->monval = get_pdg_id(urqmd_pid,
                                                            urqmd_iso3);
                    if (pick_flag == 1) {
                        if (reject_decay_flag == 2 && parent_proc_type == 20) {
                            double tau = sqrt(
                                temp_particle_info->t*temp_particle_info->t
                                - temp_particle_info->z*temp_particle_info->z);
                            if (tau > tau_reject) {
                                pick_flag = 0;
                            } else {
                                pick_flag = 1;
                            }
                        }
                        if (pick_flag == 1) {
                            (*particle_list_mixed_event)[idx]->push_back(
                                                           *temp_particle_info);
                        }
                    } else if (resonance_pick_flag == 1) {
                        (*resonance_list)[idx]->push_back(*temp_particle_info);
                    }
                    delete temp_particle_info;
                }
            }
        } else {
            break;
        }
    }
    return(0);
}

int particleSamples::read_in_particle_samples_mixed_event_gzipped() {
    // clean out the previous record
    for (unsigned int i = 0; i < particle_list_mixed_event->size(); i++)
        (*particle_list_mixed_event)[i]->clear();
    particle_list_mixed_event->clear();
    
    string temp_string;
    int n_particle;
    int temp_monval;
    for (int ievent = 0; ievent < event_buffer_size; ievent++) {
        temp_string = gz_readline(inputfile_mixed_event_gz);
        if (!gzeof(inputfile_mixed_event_gz)) {
            particle_list_mixed_event->push_back(new vector<particle_info> );
            // get number of particles within the event
            stringstream temp1(temp_string);
            temp1 >> n_particle;
            // clean out the previous record
            int idx = ievent;
            (*particle_list_mixed_event)[idx]->clear(); 

            int pick_flag = 0;
            for (int ipart = 0; ipart < n_particle; ipart++) {
                temp_string = gz_readline(inputfile_mixed_event_gz);
                stringstream temp2(temp_string);
                temp2 >> temp_monval;

                pick_flag = decide_to_pick_OSCAR(temp_monval);

                if (pick_flag != 0) {
                    particle_info *temp_particle_info = new particle_info;
                    temp_particle_info->monval = temp_monval;
                    temp2 >> temp_particle_info->mass
                          >> temp_particle_info->t
                          >> temp_particle_info->x
                          >> temp_particle_info->y
                          >> temp_particle_info->z
                          >> temp_particle_info->E
                          >> temp_particle_info->px
                          >> temp_particle_info->py
                          >> temp_particle_info->pz;

                    (*particle_list_mixed_event)[idx]->push_back(
                                                           *temp_particle_info);
                    delete temp_particle_info;
                }
            }
        } else {
            break;
        }
    }
    return(0);
}

int particleSamples::read_in_particle_samples_UrQMD_3p3_mixed_event() {
    // clean out the previous record
    for (unsigned int i = 0; i < particle_list_mixed_event->size(); i++)
        (*particle_list_mixed_event)[i]->clear();
    particle_list_mixed_event->clear();

    if (resonance_weak_feed_down_flag == 1) {
        for (unsigned int i = 0; i < resonance_list->size(); i++)
            (*resonance_list)[i]->clear();
        resonance_list->clear();
    }

    string temp_string;
    int n_particle;
    double dummy;
    int parent_proc_type;
    int ievent;
    int urqmd_pid, urqmd_iso3, urqmd_charge;
    double temp_mass;
    for (ievent = 0; ievent < event_buffer_size; ievent++) {
        getline(inputfile_mixed_event, temp_string);
        if (!inputfile_mixed_event.eof()) {
            particle_list_mixed_event->push_back(new vector<particle_info> );
            if (resonance_weak_feed_down_flag == 1) {
                resonance_list->push_back(new vector<particle_info>);
            }

            // first skip the header
            for (int i = 0; i < 13; i++)
                getline(inputfile_mixed_event, temp_string);
            // then get number of particles within the event
            getline(inputfile_mixed_event, temp_string);
            stringstream temp1(temp_string);
            temp1 >> n_particle;
            // then get one useless line
            getline(inputfile_mixed_event, temp_string);

            // clean out the previous record
            int idx = ievent;
            (*particle_list_mixed_event)[idx]->clear();
            if (resonance_weak_feed_down_flag == 1) {
                (*resonance_list)[idx]->clear();
            }

            int pick_flag = 0;
            int resonance_pick_flag = 0;
            for (int ipart = 0; ipart < n_particle; ipart++) {
                getline(inputfile_mixed_event, temp_string);
                stringstream temp2(temp_string);
                temp2 >> dummy >> dummy >> dummy >> dummy
                      >> dummy >> dummy >> dummy >> dummy
                      >> temp_mass >> urqmd_pid >> urqmd_iso3 >> urqmd_charge
                      >> dummy >> dummy >> parent_proc_type;

                pick_flag = decide_to_pick_UrQMD(
                        urqmd_pid, urqmd_iso3, urqmd_charge, parent_proc_type);

                if (resonance_weak_feed_down_flag == 1) {
                    resonance_pick_flag = decide_to_pick_UrQMD_resonance(
                                        urqmd_pid, urqmd_iso3, urqmd_charge);
                }
                
                int trigger = pick_flag + resonance_pick_flag;
                if (trigger > 0) {
                    particle_info *temp_particle_info = new particle_info;
                    temp2 >> temp_particle_info->t
                          >> temp_particle_info->x 
                          >> temp_particle_info->y
                          >> temp_particle_info->z 
                          >> temp_particle_info->E
                          >> temp_particle_info->px 
                          >> temp_particle_info->py
                          >> temp_particle_info->pz ;
                    temp_particle_info->mass = temp_mass;
                    temp_particle_info->monval = get_pdg_id(urqmd_pid,
                                                            urqmd_iso3);
                    if (pick_flag == 1) {
                        if (reject_decay_flag == 2 && parent_proc_type == 20) {
                            double tau = sqrt(
                                temp_particle_info->t*temp_particle_info->t
                                - temp_particle_info->z*temp_particle_info->z);
                            if (tau > tau_reject) {
                                pick_flag = 0;
                            } else {
                                pick_flag = 1;
                            }
                        }
                        if (pick_flag == 1) {
                            (*particle_list_mixed_event)[idx]->push_back(
                                                           *temp_particle_info);
                        }
                    } else if (resonance_pick_flag == 1) {
                        (*resonance_list)[idx]->push_back(*temp_particle_info);
                    }
                    delete temp_particle_info;
                }
            }
        } else {
            break;
        }
    }
    return(0);
}

int particleSamples::read_in_particle_samples_Sangwook() {
    // clean out the previous record
    for (unsigned int i = 0; i < particle_list->size(); i++)
        (*particle_list)[i]->clear();
    particle_list->clear();

    if (resonance_weak_feed_down_flag == 1) {
        for (unsigned int i = 0; i < resonance_list->size(); i++)
            (*resonance_list)[i]->clear();
        resonance_list->clear();
    }

    if (reconst_flag == 1) {
        for (unsigned int i = 0; i < reconst_list_1->size(); i++)
            (*reconst_list_1)[i]->clear();
        reconst_list_1->clear();
        for (unsigned int i = 0; i < reconst_list_2->size(); i++)
            (*reconst_list_2)[i]->clear();
        reconst_list_2->clear();
    }

    string temp_string;
    int n_particle;
    double dummy;
    int parent_proc_type;
    int ievent;
    int urqmd_pid, urqmd_iso3, urqmd_charge;
    double temp_mass;
    for (ievent = 0; ievent < event_buffer_size; ievent++) {
        getline(inputfile, temp_string);
        if (!inputfile.eof()) {
            // create one event
            particle_list->push_back(new vector<particle_info> );
            if (resonance_weak_feed_down_flag == 1) {
                resonance_list->push_back(new vector<particle_info>);
            }
            if (reconst_flag == 1) {
                reconst_list_1->push_back(new vector<particle_info>);
                reconst_list_2->push_back(new vector<particle_info>);
            }

            // get number of particles within the event
            stringstream temp1(temp_string);
            temp1 >> n_particle;
            getline(inputfile, temp_string);  // then get one useless line

            int idx = ievent;
            (*particle_list)[idx]->clear();  // clean out the previous record
            if (resonance_weak_feed_down_flag == 1) {
                (*resonance_list)[idx]->clear();
            }
            if (reconst_flag == 1) {
                (*reconst_list_1)[idx]->clear();
                (*reconst_list_2)[idx]->clear();
            }

            int pick_flag = 0;
            int resonance_pick_flag = 0;
            int reconst_1_pick_flag = 0;
            int reconst_2_pick_flag = 0;
            for (int ipart = 0; ipart < n_particle; ipart++) {
                getline(inputfile, temp_string);
                stringstream temp2(temp_string);
                temp2 >> dummy >> dummy >> dummy >> dummy
                      >> dummy >> dummy >> dummy >> dummy
                      >> temp_mass >> urqmd_pid >> urqmd_iso3 >> urqmd_charge
                      >> dummy >> dummy >> parent_proc_type;

                pick_flag = decide_to_pick_UrQMD(
                        urqmd_pid, urqmd_iso3, urqmd_charge, parent_proc_type);

                if (resonance_weak_feed_down_flag == 1) {
                    resonance_pick_flag = decide_to_pick_UrQMD_resonance(
                                        urqmd_pid, urqmd_iso3, urqmd_charge);
                }

                if (reconst_flag == 1) {
                    decide_to_pick_UrQMD_reconst(
                        urqmd_pid, urqmd_iso3, urqmd_charge, parent_proc_type,
                        &reconst_1_pick_flag, &reconst_2_pick_flag);
                }

                int trigger = (pick_flag + resonance_pick_flag
                               + reconst_1_pick_flag + reconst_2_pick_flag);
                if (trigger > 0) {
                    particle_info *temp_particle_info = new particle_info;
                    temp2 >> temp_particle_info->t
                          >> temp_particle_info->x
                          >> temp_particle_info->y
                          >> temp_particle_info->z
                          >> temp_particle_info->E
                          >> temp_particle_info->px
                          >> temp_particle_info->py
                          >> temp_particle_info->pz;
                    temp_particle_info->mass = temp_mass;
                    temp_particle_info->monval = get_pdg_id(urqmd_pid,
                                                            urqmd_iso3);
                    if (pick_flag == 1) {
                        if (reject_decay_flag == 2 && parent_proc_type == 20) {
                            double tau = sqrt(
                                temp_particle_info->t*temp_particle_info->t
                                - temp_particle_info->z*temp_particle_info->z);
                            if (tau > tau_reject) {
                                pick_flag = 0;
                            } else {
                                pick_flag = 1;
                            }
                        }
                        if (pick_flag == 1) {
                            (*particle_list)[idx]->push_back(
                                                        *temp_particle_info);
                        }
                    } else if (resonance_pick_flag == 1) {
                        (*resonance_list)[idx]->push_back(*temp_particle_info);
                    } else if (reconst_1_pick_flag == 1) {
                        (*reconst_list_1)[idx]->push_back(*temp_particle_info);
                    } else if (reconst_2_pick_flag == 1) {
                        (*reconst_list_2)[idx]->push_back(*temp_particle_info);
                    }
                    delete temp_particle_info;
                }
            }
        } else {
            break;
        }
    }
    return(0);
}

int particleSamples::read_in_particle_samples_mixed_event_Sangwook() {
    // clean out the previous record
    for (unsigned int i = 0; i < particle_list_mixed_event->size(); i++)
        (*particle_list_mixed_event)[i]->clear();
    particle_list_mixed_event->clear();

    if (resonance_weak_feed_down_flag == 1) {
        for (unsigned int i = 0; i < resonance_list->size(); i++)
            (*resonance_list)[i]->clear();
        resonance_list->clear();
    }

    string temp_string;
    int n_particle;
    double dummy;
    int parent_proc_type;
    int ievent;
    int urqmd_pid, urqmd_iso3, urqmd_charge;
    double temp_mass;
    for (ievent = 0; ievent < event_buffer_size; ievent++) {
        getline(inputfile_mixed_event, temp_string);
        if (!inputfile_mixed_event.eof()) {
            particle_list_mixed_event->push_back(new vector<particle_info> );
            if (resonance_weak_feed_down_flag == 1) {
                resonance_list->push_back(new vector<particle_info>);
            }

            // get number of particles within the event
            stringstream temp1(temp_string);
            temp1 >> n_particle;
            // then get one useless line
            getline(inputfile_mixed_event, temp_string);

            // clean out the previous record
            int idx = ievent;
            (*particle_list_mixed_event)[idx]->clear();
            if (resonance_weak_feed_down_flag == 1) {
                (*resonance_list)[idx]->clear();
            }

            int pick_flag = 0;
            int resonance_pick_flag = 0;
            for (int ipart = 0; ipart < n_particle; ipart++) {
                getline(inputfile_mixed_event, temp_string);
                stringstream temp2(temp_string);
                temp2 >> dummy >> dummy >> dummy >> dummy
                      >> dummy >> dummy >> dummy >> dummy
                      >> temp_mass >> urqmd_pid >> urqmd_iso3 >> urqmd_charge
                      >> dummy >> dummy >> parent_proc_type;

                pick_flag = decide_to_pick_UrQMD(
                        urqmd_pid, urqmd_iso3, urqmd_charge, parent_proc_type);

                if (resonance_weak_feed_down_flag == 1) {
                    resonance_pick_flag = decide_to_pick_UrQMD_resonance(
                                        urqmd_pid, urqmd_iso3, urqmd_charge);
                }

                int trigger = pick_flag + resonance_pick_flag;
                if (trigger > 0) {
                    particle_info *temp_particle_info = new particle_info;
                    temp2 >> temp_particle_info->t
                          >> temp_particle_info->x
                          >> temp_particle_info->y
                          >> temp_particle_info->z
                          >> temp_particle_info->E
                          >> temp_particle_info->px
                          >> temp_particle_info->py
                          >> temp_particle_info->pz;
                    temp_particle_info->mass = temp_mass;
                    temp_particle_info->monval = get_pdg_id(urqmd_pid,
                                                            urqmd_iso3);
                    if (pick_flag == 1) {
                        if (reject_decay_flag == 2 && parent_proc_type == 20) {
                            double tau = sqrt(
                                temp_particle_info->t*temp_particle_info->t
                                - temp_particle_info->z*temp_particle_info->z);
                            if (tau > tau_reject) {
                                pick_flag = 0;
                            } else {
                                pick_flag = 1;
                            }
                        }
                        if (pick_flag == 1) {
                            (*particle_list_mixed_event)[idx]->push_back(
                                                        *temp_particle_info);
                        }
                    } else if (resonance_pick_flag == 1) {
                        (*resonance_list)[idx]->push_back(*temp_particle_info);
                    }
                    delete temp_particle_info;
                }
            }
        } else {
            break;
        }
    }
    return(0);
}

void particleSamples::perform_resonance_feed_down(
                    vector< vector<particle_info>* >* input_particle_list) {
    cout << "perform resonance decays... " << endl;
    // loop over events
    unsigned int nev = input_particle_list->size();
    for (unsigned int ievent = 0; ievent < nev; ievent++) {
        // create a temporary particle list
        vector<particle_info> temp_list;
        // copy all particles into the temp list
        unsigned int Npart = (*input_particle_list)[ievent]->size();
        for (unsigned int ipart = 0; ipart < Npart; ipart++) {
            temp_list.push_back((*(*input_particle_list)[ievent])[ipart]);
        }
        (*input_particle_list)[ievent]->clear();
        // perform resonance decays
        for (unsigned int ipart = 0; ipart < temp_list.size(); ipart++) {
            vector<particle_info> *daughter_list = new vector<particle_info>;
            decayer_ptr->perform_decays(&temp_list[ipart], daughter_list);
            for (unsigned int idaughter = 0; idaughter < daughter_list->size();
                    idaughter++) {
                if (decayer_ptr->check_particle_stable(
                                        &(*daughter_list)[idaughter]) == 1) {
                    int flag = decide_to_pick_OSCAR(
                                (*daughter_list)[idaughter].monval);
                    if (flag == 1) {
                        (*input_particle_list)[ievent]->push_back(
                                            (*daughter_list)[idaughter]);
                    }
                } else {
                    temp_list.push_back((*daughter_list)[idaughter]);
                }
            }
            daughter_list->clear();
            delete daughter_list;
        }
        temp_list.clear();
    }
}

void particleSamples::perform_weak_resonance_feed_down() {
    if (particle_monval == 3122) {
        // consider Sigma^0 feed down to Lambda
        for (unsigned int iev = 0; iev < resonance_list->size(); iev++) {
            for (unsigned int i = 0; i < (*resonance_list)[iev]->size(); i++) {
                particle_info *daughter1 = new particle_info;
                particle_info *daughter2 = new particle_info;
                daughter1->mass = 1.116;  // mass of Lambda
                daughter2->mass = 0.0;    // mass of photon
                decayer_ptr->perform_two_body_decay(
                        &(*(*resonance_list)[iev])[i],
                        daughter1, daughter2);
                (*particle_list)[iev]->push_back(*daughter1);
                delete daughter2;  // discard the photon
            }
        }
    }
}

void particleSamples::perform_particle_reconstruction() {
    if (particle_monval == 333) {
        // particle of interest is phi(1020)
        double particle_mass = 1.019;
        double particle_width = 0.00443;

        // now we loop over events
        for (unsigned int iev = 0; iev < reconst_list_1->size(); iev++) {
            // we first perfrom the (K^+, K^-) pair
            for (unsigned int i = 0; i < (*reconst_list_1)[iev]->size(); i++) {
                double E_1 = (*(*reconst_list_1)[iev])[i].E;
                double px_1 = (*(*reconst_list_1)[iev])[i].px;
                double py_1 = (*(*reconst_list_1)[iev])[i].py;
                double pz_1 = (*(*reconst_list_1)[iev])[i].pz;
                double t_1 = (*(*reconst_list_1)[iev])[i].t;
                double x_1 = (*(*reconst_list_1)[iev])[i].x;
                double y_1 = (*(*reconst_list_1)[iev])[i].y;
                double z_1 = (*(*reconst_list_1)[iev])[i].z;
                for (unsigned int j = 0; j < (*reconst_list_2)[iev]->size();
                        j++) {
                    double E_2 = (*(*reconst_list_2)[iev])[j].E;
                    double px_2 = (*(*reconst_list_2)[iev])[j].px;
                    double py_2 = (*(*reconst_list_2)[iev])[j].py;
                    double pz_2 = (*(*reconst_list_2)[iev])[j].pz;
                    double t_2 = (*(*reconst_list_2)[iev])[j].t;
                    double x_2 = (*(*reconst_list_2)[iev])[j].x;
                    double y_2 = (*(*reconst_list_2)[iev])[j].y;
                    double z_2 = (*(*reconst_list_2)[iev])[j].z;

                    // compute the invariant mass
                    double E = E_1 + E_2;
                    double px = px_1 + px_2;
                    double py = py_1 + py_2;
                    double pz = pz_1 + pz_2;
                    double invariant_mass = sqrt(E*E - px*px - py*py - pz*pz);
                    if (fabs(invariant_mass - particle_mass)
                        < particle_width) {
                        // phi(1020) resonance found
                        double spatial_distance = sqrt(
                            (t_1 - t_2)*(t_1 - t_2)
                            + (x_1 - x_2)*(x_1 - x_2)
                            + (y_1 - y_2)*(y_1 - y_2)
                            + (z_1 - z_2)*(z_1 - z_2));
                        if (spatial_distance < 0.01) {
                            // this is a real phi(1020)
                            particle_info *mother = new particle_info;
                            mother->mass = particle_mass;
                            mother->E = E;
                            mother->px = px;
                            mother->py = py;
                            mother->pz = pz;
                            mother->t = (t_1 + t_2)/2.;
                            mother->x = (x_1 + x_2)/2.;
                            mother->y = (y_1 + y_2)/2.;
                            mother->z = (z_1 + z_2)/2.;
                            (*particle_list)[iev]->push_back(*mother);
                        }
                    }
                }
            }
        }
    }
    return;
}

string particleSamples::gz_readline(gzFile gzfp) {
    stringstream line;
    char buffer[1];
    int len = gzread(gzfp, buffer, 1);
    while (len == 1 && buffer[0] != '\n') {
        line << buffer[0];
        len = gzread(gzfp, buffer, 1);
    }
    return(line.str());
}
