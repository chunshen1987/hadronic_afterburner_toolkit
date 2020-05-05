// Copyright Chun Shen @ 2016
#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <cstdlib>
#include <vector>
#include <list>

#include "string.h"
#include "zlib.h"
#include "particleSamples.h"

using namespace std;

particleSamples::particleSamples(ParameterReader &paraRdr, std::string path,
                                 std::shared_ptr<RandomUtil::Random> ran_gen) :
        paraRdr_(paraRdr), path_(path) {

    echo_level_       = paraRdr_.getVal("echo_level");
    event_buffer_size = paraRdr_.getVal("event_buffer_size");
    read_in_mode_     = paraRdr_.getVal("read_in_mode");
    run_mode_         = paraRdr_.getVal("run_mode");

    read_mixed_events = false;
    if (run_mode_ == 1 || run_mode_ == 3) {
        read_mixed_events = true;
    }

    rap_type_ = paraRdr_.getVal("rap_type");

    if (run_mode_ == 2) {
        net_particle_flag = paraRdr_.getVal("net_particle_flag");
        if (net_particle_flag == 1) {
            anti_particle_list = new vector< vector<particle_info>* >;
        }
    } else {
        net_particle_flag = 0;
    }

    // read in particle Monte-Carlo number
    particle_monval = paraRdr_.getVal("particle_monval");
    flag_isospin    = paraRdr_.getVal("distinguish_isospin");

    if (run_mode_ == 3) {
        particle_monval_a    = paraRdr_.getVal("particle_alpha");
        particle_monval_abar = -particle_monval_a;
        particle_monval_b    = paraRdr_.getVal("particle_beta");
        particle_monval_bbar = -particle_monval_b;
        balance_function_particle_a    = new vector< vector<particle_info>* >;
        balance_function_particle_abar = new vector< vector<particle_info>* >;
        balance_function_particle_b    = new vector< vector<particle_info>* >;
        balance_function_particle_bbar = new vector< vector<particle_info>* >;
        balance_function_particle_a_mixed_event    = new vector< vector<particle_info>* >;
        balance_function_particle_abar_mixed_event = new vector< vector<particle_info>* >;
        balance_function_particle_b_mixed_event    = new vector< vector<particle_info>* >;
        balance_function_particle_bbar_mixed_event = new vector< vector<particle_info>* >;
    }

    resonance_weak_feed_down_flag = (
        paraRdr_.getVal("resonance_weak_feed_down_flag"));
    decayer_ptr = new particle_decay(ran_gen, resonance_weak_feed_down_flag);
    resonance_feed_down_flag = paraRdr_.getVal("resonance_feed_down_flag");
    if (resonance_weak_feed_down_flag == 1) {
        resonance_feed_down_flag = 1;
    }
    select_resonances_flag = 0;
    if (resonance_feed_down_flag == 1) {
        select_resonances_flag = paraRdr_.getVal("select_resonances_flag");
        if (select_resonances_flag == 1) {
            initialize_selected_resonance_list();
        }
    }

    full_particle_list             = new vector< vector<particle_info>* >;
    full_particle_list_mixed_event = new vector< vector<particle_info>* >;
    particle_list                  = new vector< vector<particle_info>* >;
    particle_list_mixed_event      = new vector< vector<particle_info>* >;
    if (std::abs(particle_monval) == 3122) {
        // for Lambda and anti-Lambda
        resonance_weak_feed_down_Sigma_to_Lambda_flag = paraRdr_.getVal(
                            "resonance_weak_feed_down_Sigma_to_Lambda_flag");
        if (resonance_weak_feed_down_Sigma_to_Lambda_flag == 1) {
            // include Sigma0 feed down to Lambda
            resonance_list_Sigma0 = new vector< vector<particle_info>* >;
        }
    } else {
        resonance_weak_feed_down_Sigma_to_Lambda_flag = 0;
    }

    if (particle_monval == 333) {
        // for phi(1020) meson, we need reconstruct them from (K^+ K^-) pairs
        reconst_flag   = 1;
        reconst_list_1 = new vector< vector<particle_info>* >;
        reconst_list_2 = new vector< vector<particle_info>* >;
    } else {
        reconst_flag = 0;
    }

    flag_charge_dependence = paraRdr_.getVal("flag_charge_dependence");
    if (flag_charge_dependence == 1) {
        positive_charge_hadron_list = new vector< vector<particle_info>* >;
        negative_charge_hadron_list = new vector< vector<particle_info>* >;
    }

    ostringstream filename;
    ostringstream filename_mixed_event;
    if (read_in_mode_ == 0) {
        filename << path_ << "/OSCAR.DAT";
        filename_mixed_event << path_ << "/OSCAR_mixed_event.DAT";
    } else if (read_in_mode_ == 1 || read_in_mode_ == 2 || read_in_mode_ == 3
               || read_in_mode_ == 4 || read_in_mode_ == 5
               || read_in_mode_ == 7) {
        filename << path_ << "/particle_list.dat";
        filename_mixed_event << path_ << "/particle_list_mixed_event.dat";
    } else if (read_in_mode_ == 10) {
        filename << path_ << "/particle_samples.gz";
        filename_mixed_event << path_ << "/particle_samples_mixed_event.gz";
    } else if (read_in_mode_ == 8) {
        filename << path_ << "/particles_binary.bin";
        filename_mixed_event << path_ << "/particles_binary_mixed_event.bin";
    }

    if (read_in_mode_ == 8) {
       open_SMASH_binary(filename.str(), inputfile);
       if (read_mixed_events) {
           open_SMASH_binary(filename_mixed_event.str(), inputfile_mixed_event);
       }
    } else if (read_in_mode_ != 2 && read_in_mode_ != 7 && read_in_mode_ != 10) {
        inputfile.open(filename.str().c_str());
        if (!inputfile.is_open()) {
            messager << "particleSamples:: Error: input file: "
                     << filename.str() << " can not open!";
            messager.flush("error");
            exit(1);
        }
        if (read_mixed_events) {
            inputfile_mixed_event.open(filename_mixed_event.str().c_str());
            if (!inputfile_mixed_event.is_open()) {
                messager << "particleSamples:: Error: input file: " 
                         << filename_mixed_event.str() << " can not open!";
                messager.flush("error");
                exit(1);
            }
        }
    } else {
        inputfile_gz = gzopen(filename.str().c_str(), "rb");
        if (!inputfile_gz) {
            messager << "particleSamples:: Error: input file: "
                     << filename.str() << " can not open!";
            messager.flush("error");
            exit(1);
        }
        if (read_mixed_events) {
            inputfile_mixed_event_gz =
                            gzopen(filename_mixed_event.str().c_str(), "rb");
            if (!inputfile_mixed_event_gz) {
                messager << "particleSamples:: Error: input file: " 
                         << filename_mixed_event.str() << " can not open!";
                messager.flush("error");
                exit(1);
            }
        }
    }

    // skip the header file for OSCAR
    string temp;
    if (read_in_mode_ == 0) {
        getline(inputfile, temp);
        getline(inputfile, temp);
        getline(inputfile, temp);
        if (read_mixed_events) {
            getline(inputfile_mixed_event, temp);
            getline(inputfile_mixed_event, temp);
            getline(inputfile_mixed_event, temp);
        }
    }

    // skip the header in JAM
    if (read_in_mode_ == 5) {
        getline(inputfile, temp);
        if (read_mixed_events) {
            getline(inputfile_mixed_event, temp);
        }
    }

    build_map_urqmd_to_pdg_id();
}

particleSamples::~particleSamples() {
    if (read_in_mode_ != 2 && read_in_mode_ != 7 && read_in_mode_ != 10) {
        inputfile.close();
    } else {
        gzclose(inputfile_gz);
    }
    if (read_mixed_events) {
        if (read_in_mode_ != 2 && read_in_mode_ !=7 && read_in_mode_ != 10) {
            inputfile_mixed_event.close();
        } else {
            gzclose(inputfile_mixed_event_gz);
        }
    }

    clear_out_previous_record(full_particle_list);
    delete full_particle_list;

    clear_out_previous_record(full_particle_list_mixed_event);
    delete full_particle_list_mixed_event;

    clear_out_previous_record(particle_list);
    delete particle_list;

    clear_out_previous_record(particle_list_mixed_event);
    delete particle_list_mixed_event;

    if (net_particle_flag == 1) {
        clear_out_previous_record(anti_particle_list);
        delete anti_particle_list;
    }

    delete decayer_ptr;
    if (resonance_feed_down_flag == 1) {
        if (select_resonances_flag == 1)
            select_resonances_list.clear();
    } else if (resonance_weak_feed_down_Sigma_to_Lambda_flag == 1) {
        clear_out_previous_record(resonance_list_Sigma0);
        delete resonance_list_Sigma0;
    }

    if (reconst_flag == 1) {
        clear_out_previous_record(reconst_list_1);
        delete reconst_list_1;
        clear_out_previous_record(reconst_list_2);
        delete reconst_list_2;
    }

    if (flag_charge_dependence == 1) {
        clear_out_previous_record(positive_charge_hadron_list);
        delete positive_charge_hadron_list;
        clear_out_previous_record(negative_charge_hadron_list);
        delete negative_charge_hadron_list;
    }

    if (run_mode_ == 3) {
        clear_out_previous_record(balance_function_particle_a);
        clear_out_previous_record(balance_function_particle_abar);
        clear_out_previous_record(balance_function_particle_b);
        clear_out_previous_record(balance_function_particle_bbar);
        clear_out_previous_record(balance_function_particle_a_mixed_event);
        clear_out_previous_record(balance_function_particle_abar_mixed_event);
        clear_out_previous_record(balance_function_particle_b_mixed_event);
        clear_out_previous_record(balance_function_particle_bbar_mixed_event);
    }
}

void particleSamples::open_SMASH_binary(const std::string &SMASH_filename,
                                        std::ifstream &SMASH_inputfile) {
   SMASH_inputfile.open(SMASH_filename.c_str(), ios::binary | ios::in);
   if (!SMASH_inputfile.is_open()) {
       throw std::runtime_error("Can't open file " + SMASH_filename);
   }

   // Read header
   char magic_number[5], smash_version[20];
   uint16_t format_variant;
   uint32_t len;
   SMASH_inputfile.read(&magic_number[0], 4 * sizeof(char));
   SMASH_inputfile.read(reinterpret_cast<char *>(&smash_format_version_),
                  sizeof(std::uint16_t));
   SMASH_inputfile.read(reinterpret_cast<char *>(&format_variant),
                  sizeof(std::uint16_t));
   SMASH_inputfile.read(reinterpret_cast<char *>(&len), sizeof(std::uint32_t));
   SMASH_inputfile.read(&smash_version[0], sizeof(char) * len);

   if (strcmp(magic_number, "SMSH") != 0) {
       throw std::runtime_error(
           SMASH_filename + " is likely not a SMASH binary:" +
           " magic number does not match. magic_number = "
           + magic_number);
   }

   if (format_variant != 1) {
       throw std::runtime_error(SMASH_filename + " is not a file of" +
                                " extended SMASH binary format.");
   }

   std::cout << magic_number << " " << smash_format_version_ << " "
             << format_variant << " " << smash_version << std::endl;
 
}

void particleSamples::build_map_urqmd_to_pdg_id() {
    // mesorns
    urqmd_to_pdg[std::make_pair( 101,  2)] =  211;  // pi^+
    urqmd_to_pdg[std::make_pair( 101,  0)] =  111;  // pi^0
    urqmd_to_pdg[std::make_pair( 101, -2)] = -211;  // pi^-
    urqmd_to_pdg[std::make_pair( 106,  1)] =  321;  // K^+
    urqmd_to_pdg[std::make_pair( 106, -1)] =  311;  // K^0
    urqmd_to_pdg[std::make_pair(-106,  1)] = -311;  // anti K^0
    urqmd_to_pdg[std::make_pair(-106, -1)] = -321;  // K^-
    urqmd_to_pdg[std::make_pair( 109,  0)] =  333;  // phi(1020)
    urqmd_to_pdg[std::make_pair( 102,  0)] =  221;  // eta
    urqmd_to_pdg[std::make_pair( 100,  0)] =   22;  // photon

    // baryons
    urqmd_to_pdg[std::make_pair(  1,  1)] =  2212;  // p
    urqmd_to_pdg[std::make_pair(  1, -1)] =  2112;  // n
    urqmd_to_pdg[std::make_pair( -1, -1)] = -2212;  // anti p
    urqmd_to_pdg[std::make_pair( -1,  1)] = -2112;  // anti n
    urqmd_to_pdg[std::make_pair( 40,  2)] =  3222;  // Sigma^+
    urqmd_to_pdg[std::make_pair(-40, -2)] = -3222;  // anti Sigma^-
    urqmd_to_pdg[std::make_pair( 40,  0)] =  3212;  // Sigma^0
    urqmd_to_pdg[std::make_pair(-40,  0)] = -3212;  // anti Sigma^0
    urqmd_to_pdg[std::make_pair( 40, -2)] =  3112;  // Sigma^-
    urqmd_to_pdg[std::make_pair(-40,  2)] = -3112;  // anti Sigma^+
    urqmd_to_pdg[std::make_pair( 49,  1)] =  3322;  // Xi^0
    urqmd_to_pdg[std::make_pair(-49, -1)] = -3322;  // anti Xi^0
    urqmd_to_pdg[std::make_pair( 49, -1)] =  3312;  // Xi^-
    urqmd_to_pdg[std::make_pair(-49,  1)] = -3312;  // anti Xi^+
    urqmd_to_pdg[std::make_pair( 27,  0)] =  3122;  // Labmda
    urqmd_to_pdg[std::make_pair(-27,  0)] = -3122;  // anti Labmda
    urqmd_to_pdg[std::make_pair( 55,  0)] =  3334;  // Omega
    urqmd_to_pdg[std::make_pair(-55,  0)] = -3334;  // anti Omega
}


void particleSamples::initialize_selected_resonance_list() {
    ostringstream filename;
    filename << "EOS/selected_resonances_list.dat";
    ifstream reso_file(filename.str().c_str());
    if (!reso_file.is_open()) {
        messager << "particleSamples::initialize_selected_resonance_list:"
                 << "Can not find the selected resonance list file: "
                 << filename.str();
        messager.flush("error");
        exit(1);
    }
    int temp_id;
    reso_file >> temp_id;
    while (!reso_file.eof()) {
        select_resonances_list.push_back(temp_id);
        reso_file >> temp_id;
    }
    reso_file.close();
    if (select_resonances_flag == 1 && echo_level_ > 8) {
        messager << "particleSamples::initialize_selected_resonance_list:"
                 << "selected resonance list: ";
        messager.flush("info");
        for (unsigned int ireso = 0; ireso < select_resonances_list.size();
                ireso++) {
            messager << select_resonances_list[ireso];
            messager.flush("info");
        }
    }
}


int particleSamples::get_pdg_id(int urqmd_id, int urqmd_isospin) {
    int monval = urqmd_to_pdg[std::make_pair(urqmd_id, urqmd_isospin)];
    return(monval);
}


int particleSamples::read_in_particle_samples() {
    if (read_in_mode_ == 0) {
        read_in_particle_samples_OSCAR();
        resonance_weak_feed_down_flag = 0;
        reconst_flag = 0;
    } else if (read_in_mode_ == 1) {
        read_in_particle_samples_UrQMD();
    } else if (read_in_mode_ == 2) {
        read_in_particle_samples_UrQMD_zipped();
    } else if (read_in_mode_ == 3) {
        read_in_particle_samples_Sangwook();
    } else if (read_in_mode_ == 4) {
        read_in_particle_samples_UrQMD_3p3();
    } else if (read_in_mode_ == 5) {
        read_in_particle_samples_JAM();
    } else if (read_in_mode_ == 7) {
        read_in_particle_samples_SMASH_gzipped();
    } else if (read_in_mode_ == 8) {
        read_in_particle_samples_SMASH_binary(inputfile, full_particle_list);
    } else if (read_in_mode_ == 10) {
        read_in_particle_samples_gzipped();
        resonance_weak_feed_down_flag = 0;
        reconst_flag = 0;
    }

    for (auto &ev_i: (*full_particle_list)) {
        for (auto &part_i: (*ev_i)) {
            part_i.pT    = sqrt(part_i.px*part_i.px + part_i.py*part_i.py);
            part_i.phi_p = atan2(part_i.py, part_i.px);
            if (rap_type_ == 1) {
                part_i.rap_y = 0.5*log((part_i.E + part_i.pz)
                                       /(part_i.E - part_i.pz));
            } else {
                double p_mag = sqrt(part_i.pT*part_i.pT + part_i.pz*part_i.pz);
                part_i.rap_y = 0.5*log((p_mag + part_i.pz)
                                       /(p_mag - part_i.pz));
            }
        }
    }

    if (resonance_feed_down_flag == 1)
        perform_resonance_feed_down(full_particle_list);
    return(full_particle_list->size());
}


void particleSamples::read_in_particle_samples_and_filter() {
    read_in_particle_samples();

    filter_particles_from_events(particle_monval);
    //filter_particles_into_lists(full_particle_list);

    //// reconst phi(1020) from K^+ and K^- pair
    //if (reconst_flag == 1)
    //    perform_particle_reconstruction();
}


int particleSamples::read_in_particle_samples_mixed_event() {
    if (read_in_mode_ == 0) {
        read_in_particle_samples_OSCAR_mixed_event();
        resonance_weak_feed_down_flag = 0;
    } else if (read_in_mode_ == 1) {
        read_in_particle_samples_UrQMD_mixed_event();
    } else if (read_in_mode_ == 2) {
        read_in_particle_samples_UrQMD_mixed_event_zipped();
    } else if (read_in_mode_ == 3) {
        read_in_particle_samples_mixed_event_Sangwook();
    } else if (read_in_mode_ == 4) {
        read_in_particle_samples_UrQMD_3p3_mixed_event();
    } else if (read_in_mode_ == 5) {
        read_in_particle_samples_JAM_mixed_event();
    } else if (read_in_mode_ == 7) {
        read_in_particle_samples_SMASH_mixed_event_gzipped();
    } else if (read_in_mode_ == 8) {
        read_in_particle_samples_SMASH_binary(inputfile_mixed_event,
                                              full_particle_list_mixed_event);
    } else if (read_in_mode_ == 10) {
        read_in_particle_samples_mixed_event_gzipped();
    }

    for (auto &ev_i: (*full_particle_list_mixed_event)) {
        for (auto &part_i: (*ev_i)) {
            part_i.pT    = sqrt(part_i.px*part_i.px + part_i.py*part_i.py);
            part_i.phi_p = atan2(part_i.py, part_i.px);
            if (rap_type_ == 1) {
                part_i.rap_y = 0.5*log((part_i.E + part_i.pz)
                                       /(part_i.E - part_i.pz));
            } else {
                double p_mag = sqrt(part_i.pT*part_i.pT + part_i.pz*part_i.pz);
                part_i.rap_y = 0.5*log((p_mag + part_i.pz)
                                       /(p_mag - part_i.pz));
            }
        }
    }

    if (resonance_feed_down_flag == 1)
        perform_resonance_feed_down(full_particle_list_mixed_event);

    return(full_particle_list_mixed_event->size());
}


void particleSamples::read_in_particle_samples_mixed_event_and_filter() {
    read_in_particle_samples_mixed_event();

    filter_particles(particle_monval, full_particle_list_mixed_event,
                     particle_list_mixed_event);
    if (run_mode_ == 3) {
        filter_particles(particle_monval_a, full_particle_list_mixed_event,
                         balance_function_particle_a_mixed_event);
        filter_particles(particle_monval_abar, full_particle_list_mixed_event,
                         balance_function_particle_abar_mixed_event);
        filter_particles(particle_monval_b, full_particle_list_mixed_event,
                         balance_function_particle_b_mixed_event);
        filter_particles(particle_monval_bbar, full_particle_list_mixed_event,
                         balance_function_particle_bbar_mixed_event);
    }
}


int particleSamples::decide_to_pick_anti_particles(int monval) {
    // this function judge whether particle is the anti-particle of the
    // particle of interest
    int pick_flag = 0;
    if (particle_monval == 9998) {
        // anti-particles for all positive charged hadrons
        // pick all negative charged hadrons
        int charge = decayer_ptr->get_particle_charge(monval);
        if (charge < 0) {
            pick_flag = 1;
        }
    } else if (particle_monval == 9997) {
        // anti-particles for all baryons:
        // pick all anti-baryons
        int baryon = decayer_ptr->get_particle_baryon_number(monval);
        if (baryon < 0) {
            pick_flag = 1;
        }
    } else if (particle_monval == 9996) {
        // anti-particles for all strangness particles:
        // pick all anti-strangness particles
        int strange = decayer_ptr->get_particle_strange_number(monval);
        if (strange > 0) {
            pick_flag = 1;
        }
    } else {
        // for identified particle: pick its anti-particle
        if (particle_monval == - monval) {
            pick_flag = 1;
        }
    }
    return(pick_flag);
}



int particleSamples::decide_to_pick_charge(int monval) {
    int pick_flag = 0;
    int charge = decayer_ptr->get_particle_charge(monval);
    if (charge < 0)
        pick_flag = -1;
    else if (charge > 0)
        pick_flag = 1;
    return(pick_flag);
}



int particleSamples::decide_to_pick_resonance(int monval) {
    int pick_flag = 0;
    if (particle_monval == 3122) {
        // Lambda particles, we will consider feed down from Sigma^0
        if (monval == 3212)
            pick_flag = 1;
    } else if (particle_monval == -3122) {
        // Anti-Lambda particles, we will consider feed down from Anti-Sigma^0
        if (monval == -3212)
            pick_flag = 1;
    }
    return(pick_flag);
}


int particleSamples::decide_to_pick_reconst(int monval) {
    int flag = 0;
    if (particle_monval == 333) {
        // particle of interest is phi(1020)
        // collect (K^+, K^-) pairs for reconstruction
        if (monval == 321) {   // current particle is K^+ from a decay
            flag = 1;
        }
        if (monval == -321) {  // current particle is K^- from a decay
            flag = 2;
        }
    }
    return(flag);
}


bool particleSamples::decide_to_pick_OSCAR(int POI, int monval) {
    bool pick_flag = false;
    if (POI == 9999) {
        int charge = decayer_ptr->get_particle_charge(monval);
        if (charge != 0) {
            pick_flag = true;
        }
    } else if (POI == 9998) {
        int charge = decayer_ptr->get_particle_charge(monval);
        if (charge > 0) {
            pick_flag = true;
        }
    } else if (POI == -9998) {
        int charge = decayer_ptr->get_particle_charge(monval);
        if (charge < 0) {
            pick_flag = true;
        }
    } else if (POI == 9997) {
        int baryon = decayer_ptr->get_particle_baryon_number(monval);
        if (baryon > 0) {
            pick_flag = true;
        }
    } else if (POI == -9997) {
        int baryon = decayer_ptr->get_particle_baryon_number(monval);
        if (baryon < 0) {
            pick_flag = true;
        }
    } else if (POI == 9996) {
        int strange = decayer_ptr->get_particle_strange_number(monval);
        if (strange > 0) {
            pick_flag = true;
        }
    } else if (POI == -9996) {
        int strange = decayer_ptr->get_particle_strange_number(monval);
        if (strange < 0) {
            pick_flag = true;
        }
    } else {
        if (monval == POI) {
            pick_flag = true;
        }
    }
    return(pick_flag);
}

int particleSamples::read_in_particle_samples_OSCAR() {
    // clean out the previous record
    clear_out_previous_record(full_particle_list);

    std::string temp_string;
    int event_id, n_particle, dummy;
    int ievent = 0;
    int temp_monval;
    int num_particles = 0;
    while (num_particles < event_buffer_size) {
        getline(inputfile, temp_string);
        std::stringstream temp1(temp_string);
        temp1 >> event_id >> n_particle;

        if (inputfile.eof()) break;

        full_particle_list->push_back(new vector<particle_info> );
        for (int ipart = 0; ipart < n_particle; ipart++) {
            getline(inputfile, temp_string);
            std::stringstream temp2(temp_string);
            temp2 >> dummy >> temp_monval;
            particle_info temp_particle_info;
            temp_particle_info.monval = temp_monval;
            temp2 >> temp_particle_info.px
                  >> temp_particle_info.py
                  >> temp_particle_info.pz
                  >> temp_particle_info.E
                  >> temp_particle_info.mass
                  >> temp_particle_info.x
                  >> temp_particle_info.y
                  >> temp_particle_info.z
                  >> temp_particle_info.t;
            (*full_particle_list)[ievent]->push_back(temp_particle_info);
        }
        num_particles += n_particle;
        ievent++;
    }
    return(0);
}

int particleSamples::read_in_particle_samples_JAM() {
    // clean out the previous record
    clear_out_previous_record(full_particle_list);

    std::string temp_string;
    int event_id, n_particle;
    char cdummy;
    int ievent = 0;
    int temp_monval;
    int num_particles = 0;
    while (num_particles < event_buffer_size) {
        getline(inputfile, temp_string);
        std::stringstream temp1(temp_string);
        temp1 >> cdummy >> event_id >> n_particle;

        if (inputfile.eof()) break;

        full_particle_list->push_back(new vector<particle_info> );
        for (int ipart = 0; ipart < n_particle; ipart++) {
            getline(inputfile, temp_string);
            std::stringstream temp2(temp_string);
            temp2 >> temp_monval;

            particle_info temp_particle_info;
            temp2 >> temp_particle_info.mass
                  >> temp_particle_info.px
                  >> temp_particle_info.py
                  >> temp_particle_info.pz
                  >> temp_particle_info.x 
                  >> temp_particle_info.y
                  >> temp_particle_info.z 
                  >> temp_particle_info.t;
            temp_particle_info.E = sqrt(
                temp_particle_info.mass*temp_particle_info.mass
                + temp_particle_info.px*temp_particle_info.px
                + temp_particle_info.py*temp_particle_info.py
                + temp_particle_info.pz*temp_particle_info.pz);
            (*full_particle_list)[ievent]->push_back(temp_particle_info);
        }
        num_particles += n_particle;
        ievent++;
    }
    return(0);
}

int particleSamples::read_in_particle_samples_OSCAR_mixed_event() {
    // clean out the previous record
    clear_out_previous_record(full_particle_list_mixed_event);

    std::string temp_string;
    int event_id, n_particle, dummy;
    int ievent = 0;
    int temp_monval;
    int num_particles = 0;
    while (num_particles < event_buffer_size) {
        getline(inputfile_mixed_event, temp_string);
        std::stringstream temp1(temp_string);
        temp1 >> event_id >> n_particle;

        if (inputfile_mixed_event.eof()) break;

        full_particle_list_mixed_event->push_back(
                                            new vector<particle_info> );
        for (int ipart = 0; ipart < n_particle; ipart++) {
            getline(inputfile_mixed_event, temp_string);
            std::stringstream temp2(temp_string);
            temp2 >> dummy >> temp_monval;
            particle_info temp_particle_info;
            temp_particle_info.monval = temp_monval;
            temp2 >> temp_particle_info.px
                  >> temp_particle_info.py
                  >> temp_particle_info.pz
                  >> temp_particle_info.E
                  >> temp_particle_info.mass
                  >> temp_particle_info.x
                  >> temp_particle_info.y
                  >> temp_particle_info.z
                  >> temp_particle_info.t;
            (*full_particle_list_mixed_event)[ievent]->push_back(
                                                    temp_particle_info);
        }
        num_particles += n_particle;
        ievent++;
    }
    return(0);
}

int particleSamples::read_in_particle_samples_JAM_mixed_event() {
    // clean out the previous record
    clear_out_previous_record(full_particle_list_mixed_event);

    std::string temp_string;
    int event_id, n_particle;
    char cdummy;
    int ievent = 0;
    int temp_monval;
    int num_particles = 0;
    while (num_particles < event_buffer_size) {
        getline(inputfile_mixed_event, temp_string);
        std::stringstream temp1(temp_string);
        temp1 >> cdummy >> event_id >> n_particle;

        if (inputfile_mixed_event.eof()) break;

        full_particle_list_mixed_event->push_back(
                                            new vector<particle_info> );

        for (int ipart = 0; ipart < n_particle; ipart++) {
            getline(inputfile_mixed_event, temp_string);
            std::stringstream temp2(temp_string);
            temp2 >> temp_monval;
            particle_info temp_particle_info;
            temp2 >> temp_particle_info.mass
                  >> temp_particle_info.px
                  >> temp_particle_info.py
                  >> temp_particle_info.pz
                  >> temp_particle_info.x
                  >> temp_particle_info.y
                  >> temp_particle_info.z
                  >> temp_particle_info.t;
            temp_particle_info.E = sqrt(
                temp_particle_info.mass*temp_particle_info.mass
                + temp_particle_info.px*temp_particle_info.px
                + temp_particle_info.py*temp_particle_info.py
                + temp_particle_info.pz*temp_particle_info.pz);
            (*full_particle_list_mixed_event)[ievent]->push_back(
                                                    temp_particle_info);
        }
        num_particles += n_particle;
        ievent++;
    }
    return(0);
}

int particleSamples::read_in_particle_samples_UrQMD() {
    // clean out the previous record
    clear_out_previous_record(full_particle_list);

    std::string temp_string;
    int n_particle;
    double dummy;
    int parent_proc_type;
    int ievent = 0;
    int urqmd_pid, urqmd_iso3, urqmd_charge;
    double temp_mass;
    int num_particles = 0;
    while (num_particles < event_buffer_size) {
        getline(inputfile, temp_string);

        if (inputfile.eof()) break;

        // create one event
        full_particle_list->push_back(new vector<particle_info> );

        // first skip the header
        for (int i = 0; i < 16; i++)
            getline(inputfile, temp_string);

        // then get number of particles within the event
        getline(inputfile, temp_string);

        std::stringstream temp1(temp_string);
        temp1 >> n_particle;
        getline(inputfile, temp_string);  // then get one useless line

        for (int ipart = 0; ipart < n_particle; ipart++) {
            getline(inputfile, temp_string);
            std::stringstream temp2(temp_string);
            temp2 >> dummy >> dummy >> dummy >> dummy
                  >> dummy >> dummy >> dummy >> dummy
                  >> temp_mass >> urqmd_pid >> urqmd_iso3 >> urqmd_charge
                  >> dummy >> dummy >> parent_proc_type;

            particle_info temp_particle_info;
            temp2 >> temp_particle_info.t
                  >> temp_particle_info.x
                  >> temp_particle_info.y
                  >> temp_particle_info.z
                  >> temp_particle_info.E
                  >> temp_particle_info.px
                  >> temp_particle_info.py
                  >> temp_particle_info.pz;
            temp_particle_info.mass = temp_mass;
            temp_particle_info.monval = get_pdg_id(urqmd_pid, urqmd_iso3);
            (*full_particle_list)[ievent]->push_back(temp_particle_info);
        }
        num_particles += n_particle;
        ievent++;
    }
    return(0);
}


int particleSamples::read_in_particle_samples_UrQMD_zipped() {
    // clean out the previous record
    clear_out_previous_record(full_particle_list);

    std::string temp_string;
    int n_particle;
    double dummy;
    int parent_proc_type;
    int ievent = 0;
    int urqmd_pid, urqmd_iso3, urqmd_charge;
    int num_particles = 0;
    while (num_particles < event_buffer_size) {
        temp_string = gz_readline(inputfile_gz);

        if (gzeof(inputfile_gz)) break;

        full_particle_list->push_back(new vector<particle_info> );

        // get number of particles within the event
        std::stringstream temp1(temp_string);
        temp1 >> n_particle;

        // then get one useless line
        temp_string = gz_readline(inputfile_gz);

        for (int ipart = 0; ipart < n_particle; ipart++) {
            temp_string = gz_readline(inputfile_gz);
            std::stringstream temp2(temp_string);
            temp2 >> urqmd_pid >> urqmd_iso3 >> urqmd_charge
                  >> dummy >> dummy >> parent_proc_type;

            particle_info temp_particle_info;
            temp2 >> temp_particle_info.mass
                  >> temp_particle_info.t
                  >> temp_particle_info.x
                  >> temp_particle_info.y
                  >> temp_particle_info.z
                  >> temp_particle_info.E
                  >> temp_particle_info.px
                  >> temp_particle_info.py
                  >> temp_particle_info.pz;
            temp_particle_info.monval = get_pdg_id(urqmd_pid, urqmd_iso3);
            (*full_particle_list)[ievent]->push_back(temp_particle_info);
        }
        num_particles += n_particle;
        ievent++;
    }
    return(0);
}


int particleSamples::read_in_particle_samples_SMASH_gzipped() {
    // clean out the previous record
    clear_out_previous_record(full_particle_list);

    std::string temp_string;
    int n_particle;
    int parent_proc_type;
    int pdg_mother1, pdg_mother2;
    int pdg_id, charge;
    int ievent = 0;
    int num_particles = 0;
    while (num_particles < event_buffer_size) {
        temp_string = gz_readline(inputfile_gz);

        if (gzeof(inputfile_gz)) break;

        full_particle_list->push_back(new vector<particle_info> );

        // get number of particles within the event
        std::stringstream temp1(temp_string);
        temp1 >> n_particle;

        for (int ipart = 0; ipart < n_particle; ipart++) {
            temp_string = gz_readline(inputfile_gz);
            std::stringstream temp2(temp_string);
            temp2 >> pdg_id >> charge >> parent_proc_type
                  >> pdg_mother1 >> pdg_mother2;

            particle_info temp_particle_info;
            temp2 >> temp_particle_info.mass
                  >> temp_particle_info.t
                  >> temp_particle_info.x
                  >> temp_particle_info.y
                  >> temp_particle_info.z
                  >> temp_particle_info.E
                  >> temp_particle_info.px
                  >> temp_particle_info.py
                  >> temp_particle_info.pz;
            temp_particle_info.monval = pdg_id;
            (*full_particle_list)[ievent]->push_back(temp_particle_info);
        }
        num_particles += n_particle;
        ievent++;
    }
    return(0);
}

int particleSamples::read_in_particle_samples_SMASH_binary(
    std::ifstream &SMASH_inputfile,
    std::vector< std::vector<particle_info>* >* resulting_particle_list) {
    // clean out the previous record
    clear_out_previous_record(resulting_particle_list);
    int ievent = 0, num_particles = 0;
    while (num_particles < event_buffer_size) {
      char block_type;
      SMASH_inputfile.read(&block_type, sizeof(char));
      if (!SMASH_inputfile) {
        break;
      }
      if (block_type == 'f') {
        uint32_t ev;
        double impact_parameter;
        char empty;
        SMASH_inputfile.read(reinterpret_cast<char *>(&ev), sizeof(std::uint32_t));
        SMASH_inputfile.read(reinterpret_cast<char *>(&impact_parameter), sizeof(double));
        if (smash_format_version_ > 6) {
          SMASH_inputfile.read(&empty, sizeof(char));
        }
        ievent++;
        continue;
      }
      if (block_type != 'p') {
        break;
      }
      resulting_particle_list->push_back(new vector<particle_info> );

      uint32_t n_part_lines;
      SMASH_inputfile.read(reinterpret_cast<char *>(&n_part_lines), sizeof(std::uint32_t));
      for (size_t i = 0; i < n_part_lines; i++) {
        double t, x, y, z, m, p0, px, py, pz,
               form_time, xsecfac, time_last_coll;
        int32_t pdg, id, charge, ncoll, proc_id_origin,
                proc_type_origin, pdg_mother1, pdg_mother2;
        SMASH_inputfile.read(reinterpret_cast<char *>(&t), sizeof(double));
        SMASH_inputfile.read(reinterpret_cast<char *>(&x), sizeof(double));
        SMASH_inputfile.read(reinterpret_cast<char *>(&y), sizeof(double));
        SMASH_inputfile.read(reinterpret_cast<char *>(&z), sizeof(double));
        SMASH_inputfile.read(reinterpret_cast<char *>(&m), sizeof(double));
        SMASH_inputfile.read(reinterpret_cast<char *>(&p0), sizeof(double));
        SMASH_inputfile.read(reinterpret_cast<char *>(&px), sizeof(double));
        SMASH_inputfile.read(reinterpret_cast<char *>(&py), sizeof(double));
        SMASH_inputfile.read(reinterpret_cast<char *>(&pz), sizeof(double));
        SMASH_inputfile.read(reinterpret_cast<char *>(&pdg), sizeof(std::int32_t));
        SMASH_inputfile.read(reinterpret_cast<char *>(&id),  sizeof(std::int32_t));
        SMASH_inputfile.read(reinterpret_cast<char *>(&charge), sizeof(std::int32_t));
        SMASH_inputfile.read(reinterpret_cast<char *>(&ncoll),  sizeof(std::int32_t));
        SMASH_inputfile.read(reinterpret_cast<char *>(&form_time), sizeof(double));
        SMASH_inputfile.read(reinterpret_cast<char *>(&xsecfac), sizeof(double));
        SMASH_inputfile.read(reinterpret_cast<char *>(&proc_id_origin), sizeof(std::int32_t));
        SMASH_inputfile.read(reinterpret_cast<char *>(&proc_type_origin), sizeof(std::int32_t));
        SMASH_inputfile.read(reinterpret_cast<char *>(&time_last_coll), sizeof(double));
        SMASH_inputfile.read(reinterpret_cast<char *>(&pdg_mother1), sizeof(std::int32_t));
        SMASH_inputfile.read(reinterpret_cast<char *>(&pdg_mother2), sizeof(std::int32_t));

        double origin_t = time_last_coll;
        double vx = px/p0, vy = py/p0, vz = pz/p0;
        double dt = t - origin_t;
        double origin_x = x - vx * dt,
               origin_y = y - vy * dt,
               origin_z = z - vz * dt;

        particle_info temp_particle_info;
        temp_particle_info.monval = pdg;
        temp_particle_info.mass = m;
        temp_particle_info.t = origin_t;
        temp_particle_info.x = origin_x;
        temp_particle_info.y = origin_y;
        temp_particle_info.z = origin_z;
        temp_particle_info.E  = p0;
        temp_particle_info.px = px;
        temp_particle_info.py = py;
        temp_particle_info.pz = pz;

        (*resulting_particle_list)[ievent]->push_back(temp_particle_info);
      }
      num_particles += n_part_lines;
      std::cout << "Read in " << n_part_lines << " particles" << std::endl;
    }
    return(0);
}


int particleSamples::read_in_particle_samples_gzipped() {
    // clean out the previous record
    clear_out_previous_record(full_particle_list);

    std::string temp_string;
    int n_particle;
    int temp_monval;
    int num_particles = 0;
    int ievent = 0;
    while (num_particles < event_buffer_size) {
        temp_string = gz_readline(inputfile_gz);

        if (gzeof(inputfile_gz)) break;

        // create one event
        full_particle_list->push_back(new vector<particle_info> );

        // get number of particles within the event
        std::stringstream temp1(temp_string);
        temp1 >> n_particle;

        for (int ipart = 0; ipart < n_particle; ipart++) {
            temp_string = gz_readline(inputfile_gz);
            std::stringstream temp2(temp_string);
            temp2 >> temp_monval;
            particle_info temp_particle_info;
            temp_particle_info.monval = temp_monval;
            temp2 >> temp_particle_info.mass
                  >> temp_particle_info.t
                  >> temp_particle_info.x
                  >> temp_particle_info.y
                  >> temp_particle_info.z
                  >> temp_particle_info.E
                  >> temp_particle_info.px
                  >> temp_particle_info.py
                  >> temp_particle_info.pz;

            (*full_particle_list)[ievent]->push_back(temp_particle_info);
        }
        num_particles += n_particle;
        ievent++;
    }
    return(0);
}


void particleSamples::filter_particles_from_events(const int PoI_monval) {
    filter_particles(PoI_monval, full_particle_list, particle_list);

    if (resonance_weak_feed_down_Sigma_to_Lambda_flag == 1) {
        if (PoI_monval == 3122) {
            filter_particles(3212, full_particle_list, resonance_list_Sigma0);
        } else if (PoI_monval == -3122) {
            filter_particles(-3212, full_particle_list, resonance_list_Sigma0);
        }
    }

    if (net_particle_flag == 1) {
        filter_particles(-PoI_monval, full_particle_list, anti_particle_list);
    }

    if (reconst_flag == 1 && PoI_monval == 333) {
        filter_particles(321, full_particle_list, reconst_list_1);
        filter_particles(-321, full_particle_list, reconst_list_2);
        perform_particle_reconstruction();
    }

    if (flag_charge_dependence == 1 && PoI_monval == 9999) {
        filter_particles(9998, full_particle_list,
                         positive_charge_hadron_list);
        filter_particles(-9998, full_particle_list,
                         negative_charge_hadron_list);
    }

    if (run_mode_ == 3) {
        filter_particles(particle_monval_a, full_particle_list,
                         balance_function_particle_a);
        filter_particles(particle_monval_abar, full_particle_list,
                         balance_function_particle_abar);
        filter_particles(particle_monval_b, full_particle_list,
                         balance_function_particle_b);
        filter_particles(particle_monval_bbar, full_particle_list,
                         balance_function_particle_bbar);
    }
}


void particleSamples::filter_particles(
                    const int PoI_monval,
                    vector< vector<particle_info>* >* full_list,
                    vector< vector<particle_info>* >* filted_list) {
    // clean out the previous record
    clear_out_previous_record(filted_list);

    int i = 0;
    for (auto &ev_i: (*full_list)) {
        filted_list->push_back(new vector<particle_info> );
        for (auto &part_i: (*ev_i)) {
            bool pick_flag = decide_to_pick_OSCAR(PoI_monval,
                                                  part_i.monval);
            if (pick_flag)
                (*filted_list)[i]->push_back(part_i);
        }
        i++;
    }
}


void particleSamples::clear_out_previous_record(
                    vector< vector<particle_info>* >* plist) {
    for (auto &ev_i: (*plist)) {
        ev_i->clear();
        delete ev_i;
    }
    plist->clear();
}


void particleSamples::filter_particles_into_lists(
                    vector< vector<particle_info>* >* full_list) {
    // clean out the previous record
    clear_out_previous_record(particle_list);

    if (resonance_weak_feed_down_Sigma_to_Lambda_flag == 1)
        clear_out_previous_record(resonance_list_Sigma0);

    if (net_particle_flag == 1)
        clear_out_previous_record(anti_particle_list);

    if (reconst_flag == 1) {
        clear_out_previous_record(reconst_list_1);
        clear_out_previous_record(reconst_list_2);
    }

    if (flag_charge_dependence == 1) {
        clear_out_previous_record(positive_charge_hadron_list);
        clear_out_previous_record(negative_charge_hadron_list);
    }

    if (run_mode_ == 3) {
        clear_out_previous_record(balance_function_particle_a);
        clear_out_previous_record(balance_function_particle_abar);
        clear_out_previous_record(balance_function_particle_b);
        clear_out_previous_record(balance_function_particle_bbar);
    }

    int iev = 0;
    for (auto &ev_i: (*full_list)) {
        particle_list->push_back(new vector<particle_info> );
        if (resonance_weak_feed_down_Sigma_to_Lambda_flag == 1)
            resonance_list_Sigma0->push_back(new vector<particle_info>);

        if (net_particle_flag == 1)
            anti_particle_list->push_back(new vector<particle_info>);

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

        if (run_mode_ == 3) {
            balance_function_particle_a->push_back(
                                            new vector<particle_info>);
            balance_function_particle_abar->push_back(
                                            new vector<particle_info>);
            balance_function_particle_b->push_back(
                                            new vector<particle_info>);
            balance_function_particle_bbar->push_back(
                                            new vector<particle_info>);
        }

        for (auto &part_i: (*ev_i)) {
            bool pick_flag = decide_to_pick_OSCAR(particle_monval,
                                                  part_i.monval);
            if (pick_flag)
                (*particle_list)[iev]->push_back(part_i);

            if (resonance_weak_feed_down_Sigma_to_Lambda_flag == 1) {
                int resonance_pick_flag = decide_to_pick_resonance(
                                                            part_i.monval);
                if (resonance_pick_flag == 1)
                    (*resonance_list_Sigma0)[iev]->push_back(part_i);
            }

            if (reconst_flag == 1) {
                int reconst_pick_flag = decide_to_pick_reconst(part_i.monval);
                if (reconst_pick_flag == 1)
                    (*reconst_list_1)[iev]->push_back(part_i);
                else if (reconst_pick_flag == 2)
                    (*reconst_list_2)[iev]->push_back(part_i);
            }
            if (net_particle_flag == 1) {
                int anti_particle_pick_flag = decide_to_pick_anti_particles(
                                                                part_i.monval);
                if (anti_particle_pick_flag == 1)
                    (*anti_particle_list)[iev]->push_back(part_i);
            }

            if (flag_charge_dependence == 1) {
                int charge_flag = decide_to_pick_charge(part_i.monval);
                if (charge_flag == 1)
                    (*positive_charge_hadron_list)[iev]->push_back(part_i);
                else if (charge_flag == -1)
                    (*negative_charge_hadron_list)[iev]->push_back(part_i);
            }

            if (run_mode_ == 3) {
                if (decide_to_pick_OSCAR(particle_monval_a, part_i.monval))
                    (*balance_function_particle_a)[iev]->push_back(part_i);
                if (decide_to_pick_OSCAR(particle_monval_abar, part_i.monval))
                    (*balance_function_particle_abar)[iev]->push_back(part_i);
                if (decide_to_pick_OSCAR(particle_monval_b, part_i.monval))
                    (*balance_function_particle_b)[iev]->push_back(part_i);
                if (decide_to_pick_OSCAR(particle_monval_bbar, part_i.monval))
                    (*balance_function_particle_bbar)[iev]->push_back(part_i);
            }
        }
        iev++;
    }
}


int particleSamples::read_in_particle_samples_UrQMD_3p3() {
    // clean out the previous record
    clear_out_previous_record(full_particle_list);

    std::string temp_string;
    int n_particle;
    double dummy;
    int parent_proc_type;
    int ievent = 0;
    int urqmd_pid, urqmd_iso3, urqmd_charge;
    double temp_mass;
    int num_particles = 0;
    while (num_particles < event_buffer_size) {
        getline(inputfile, temp_string);
        if (inputfile.eof()) break;
        // create one event
        full_particle_list->push_back(new vector<particle_info> );

        // first skip the header
        for (int i = 0; i < 13; i++)
            getline(inputfile, temp_string);
        // then get number of particles within the event
        getline(inputfile, temp_string);
        std::stringstream temp1(temp_string);
        temp1 >> n_particle;
        getline(inputfile, temp_string);  // then get one useless line

        for (int ipart = 0; ipart < n_particle; ipart++) {
            getline(inputfile, temp_string);
            std::stringstream temp2(temp_string);
            temp2 >> dummy >> dummy >> dummy >> dummy
                  >> dummy >> dummy >> dummy >> dummy
                  >> temp_mass >> urqmd_pid >> urqmd_iso3 >> urqmd_charge
                  >> dummy >> dummy >> parent_proc_type;

            particle_info temp_particle_info;
            temp2 >> temp_particle_info.t
                  >> temp_particle_info.x
                  >> temp_particle_info.y
                  >> temp_particle_info.z
                  >> temp_particle_info.E
                  >> temp_particle_info.px
                  >> temp_particle_info.py
                  >> temp_particle_info.pz;
            temp_particle_info.mass = temp_mass;
            temp_particle_info.monval = get_pdg_id(urqmd_pid, urqmd_iso3);
            (*full_particle_list)[ievent]->push_back(temp_particle_info);
        }
        num_particles += n_particle;
        ievent++;
    }
    return(0);
}


int particleSamples::read_in_particle_samples_UrQMD_mixed_event() {
    // clean out the previous record
    clear_out_previous_record(full_particle_list_mixed_event);

    std::string temp_string;
    int n_particle;
    double dummy;
    int parent_proc_type;
    int ievent = 0;
    int urqmd_pid, urqmd_iso3, urqmd_charge;
    double temp_mass;
    int num_particles = 0;
    while (num_particles < event_buffer_size) {
        getline(inputfile_mixed_event, temp_string);

        if (inputfile_mixed_event.eof()) break;

        full_particle_list_mixed_event->push_back(
                                            new vector<particle_info> );

        // first skip the header
        for(int i = 0; i < 16; i++)
            getline(inputfile_mixed_event, temp_string);
        // then get number of particles within the event
        getline(inputfile_mixed_event, temp_string);
        std::stringstream temp1(temp_string);
        temp1 >> n_particle;
        // then get one useless line
        getline(inputfile_mixed_event, temp_string);  

        for (int ipart = 0; ipart < n_particle; ipart++) {
            getline(inputfile_mixed_event, temp_string);
            std::stringstream temp2(temp_string);
            temp2 >> dummy >> dummy >> dummy >> dummy
                  >> dummy >> dummy >> dummy >> dummy
                  >> temp_mass >> urqmd_pid >> urqmd_iso3 >> urqmd_charge
                  >> dummy >> dummy >> parent_proc_type;

            particle_info temp_particle_info;
            temp2 >> temp_particle_info.t
                  >> temp_particle_info.x
                  >> temp_particle_info.y
                  >> temp_particle_info.z
                  >> temp_particle_info.E
                  >> temp_particle_info.px
                  >> temp_particle_info.py
                  >> temp_particle_info.pz;
            temp_particle_info.mass = temp_mass;
            temp_particle_info.monval = get_pdg_id(urqmd_pid, urqmd_iso3);
            (*full_particle_list_mixed_event)[ievent]->push_back(
                                                    temp_particle_info);
        }
        num_particles += n_particle;
        ievent++;
    }
    return(0);
}


int particleSamples::read_in_particle_samples_UrQMD_mixed_event_zipped() {
    // clean out the previous record
    clear_out_previous_record(full_particle_list_mixed_event);

    std::string temp_string;
    int n_particle;
    double dummy;
    int parent_proc_type;
    int urqmd_pid, urqmd_iso3, urqmd_charge;
    int ievent = 0;
    int num_particles = 0;
    while (num_particles < event_buffer_size) {
        temp_string = gz_readline(inputfile_mixed_event_gz);

        if (gzeof(inputfile_mixed_event_gz)) break;

        full_particle_list_mixed_event->push_back(
                                            new vector<particle_info> );

        // get number of particles within the event
        std::stringstream temp1(temp_string);
        temp1 >> n_particle;
        // then get one useless line
        temp_string = gz_readline(inputfile_mixed_event_gz);  

        for (int ipart = 0; ipart < n_particle; ipart++) {
            temp_string = gz_readline(inputfile_mixed_event_gz);
            std::stringstream temp2(temp_string);
            temp2 >> urqmd_pid >> urqmd_iso3 >> urqmd_charge
                  >> dummy >> dummy >> parent_proc_type;

            particle_info temp_particle_info;
            temp2 >> temp_particle_info.mass
                  >> temp_particle_info.t
                  >> temp_particle_info.x
                  >> temp_particle_info.y
                  >> temp_particle_info.z
                  >> temp_particle_info.E
                  >> temp_particle_info.px
                  >> temp_particle_info.py
                  >> temp_particle_info.pz;
            temp_particle_info.monval = get_pdg_id(urqmd_pid,
                                                    urqmd_iso3);
            (*full_particle_list_mixed_event)[ievent]->push_back(
                                                    temp_particle_info);
        }
        num_particles += n_particle;
        ievent++;
    }
    return(0);
}


int particleSamples::read_in_particle_samples_SMASH_mixed_event_gzipped() {
    // clean out the previous record
    clear_out_previous_record(full_particle_list_mixed_event);

    std::string temp_string;
    int n_particle;
    int parent_proc_type;
    int pdg_mother1, pdg_mother2;
    int pdg_id, charge;
    int ievent = 0;
    int num_particles = 0;
    while (num_particles < event_buffer_size) {
        temp_string = gz_readline(inputfile_mixed_event_gz);

        if (gzeof(inputfile_mixed_event_gz)) break;

        full_particle_list_mixed_event->push_back(
                                            new vector<particle_info> );

        // get number of particles within the event
        std::stringstream temp1(temp_string);
        temp1 >> n_particle;
        for (int ipart = 0; ipart < n_particle; ipart++) {
            temp_string = gz_readline(inputfile_mixed_event_gz);
            std::stringstream temp2(temp_string);
            temp2 >> pdg_id >> charge >> parent_proc_type
                  >> pdg_mother1 >> pdg_mother2;

            particle_info temp_particle_info;
            temp2 >> temp_particle_info.mass
                  >> temp_particle_info.t
                  >> temp_particle_info.x
                  >> temp_particle_info.y
                  >> temp_particle_info.z
                  >> temp_particle_info.E
                  >> temp_particle_info.px
                  >> temp_particle_info.py
                  >> temp_particle_info.pz;
            temp_particle_info.monval = pdg_id;
            (*full_particle_list_mixed_event)[ievent]->push_back(
                                                    temp_particle_info);
        }
        num_particles += n_particle;
        ievent++;
    }
    return(0);
}


int particleSamples::read_in_particle_samples_mixed_event_gzipped() {
    // clean out the previous record
    for (unsigned int i = 0; i < full_particle_list_mixed_event->size(); i++)
        (*full_particle_list_mixed_event)[i]->clear();
    full_particle_list_mixed_event->clear();

    string temp_string;
    int n_particle;
    int temp_monval;
    int ievent = 0;
    int num_particles = 0;
    while (num_particles < event_buffer_size) {
        temp_string = gz_readline(inputfile_mixed_event_gz);
        if (gzeof(inputfile_mixed_event_gz)) break;
        full_particle_list_mixed_event->push_back(
                                            new vector<particle_info> );
        // get number of particles within the event
        stringstream temp1(temp_string);
        temp1 >> n_particle;

        for (int ipart = 0; ipart < n_particle; ipart++) {
            temp_string = gz_readline(inputfile_mixed_event_gz);
            stringstream temp2(temp_string);
            temp2 >> temp_monval;

            particle_info temp_particle_info;
            temp_particle_info.monval = temp_monval;
            temp2 >> temp_particle_info.mass
                  >> temp_particle_info.t
                  >> temp_particle_info.x
                  >> temp_particle_info.y
                  >> temp_particle_info.z
                  >> temp_particle_info.E
                  >> temp_particle_info.px
                  >> temp_particle_info.py
                  >> temp_particle_info.pz;

            (*full_particle_list_mixed_event)[ievent]->push_back(
                                                   temp_particle_info);
        }
        num_particles += n_particle;
        ievent++;
    }
    return(0);
}


int particleSamples::read_in_particle_samples_UrQMD_3p3_mixed_event() {
    // clean out the previous record
    clear_out_previous_record(full_particle_list_mixed_event);

    std::string temp_string;
    int n_particle;
    double dummy;
    int parent_proc_type;
    int ievent = 0;
    int urqmd_pid, urqmd_iso3, urqmd_charge;
    double temp_mass;
    int num_particles = 0;
    while (num_particles < event_buffer_size) {
        getline(inputfile_mixed_event, temp_string);

        if (inputfile_mixed_event.eof()) break;

        full_particle_list_mixed_event->push_back(
                                            new vector<particle_info> );

        // first skip the header
        for (int i = 0; i < 13; i++)
            getline(inputfile_mixed_event, temp_string);
        // then get number of particles within the event
        getline(inputfile_mixed_event, temp_string);
        std::stringstream temp1(temp_string);
        temp1 >> n_particle;
        // then get one useless line
        getline(inputfile_mixed_event, temp_string);

        for (int ipart = 0; ipart < n_particle; ipart++) {
            getline(inputfile_mixed_event, temp_string);
            std::stringstream temp2(temp_string);
            temp2 >> dummy >> dummy >> dummy >> dummy
                  >> dummy >> dummy >> dummy >> dummy
                  >> temp_mass >> urqmd_pid >> urqmd_iso3 >> urqmd_charge
                  >> dummy >> dummy >> parent_proc_type;

            particle_info temp_particle_info;
            temp2 >> temp_particle_info.t
                  >> temp_particle_info.x
                  >> temp_particle_info.y
                  >> temp_particle_info.z
                  >> temp_particle_info.E
                  >> temp_particle_info.px
                  >> temp_particle_info.py
                  >> temp_particle_info.pz;
            temp_particle_info.mass = temp_mass;
            temp_particle_info.monval = get_pdg_id(urqmd_pid, urqmd_iso3);
            (*full_particle_list_mixed_event)[ievent]->push_back(
                                                    temp_particle_info);
        }
        num_particles += n_particle;
        ievent++;
    }
    return(0);
}


int particleSamples::read_in_particle_samples_Sangwook() {
    // clean out the previous record
    clear_out_previous_record(full_particle_list);

    std::string temp_string;
    int n_particle;
    double dummy;
    int parent_proc_type;
    int ievent = 0;
    int urqmd_pid, urqmd_iso3, urqmd_charge;
    double temp_mass;
    int num_particles = 0;
    while (num_particles < event_buffer_size) {
        getline(inputfile, temp_string);

        if (inputfile.eof()) break;

        // create one event
        full_particle_list->push_back(new vector<particle_info> );

        // get number of particles within the event
        std::stringstream temp1(temp_string);
        temp1 >> n_particle;
        getline(inputfile, temp_string);  // then get one useless line

        for (int ipart = 0; ipart < n_particle; ipart++) {
            getline(inputfile, temp_string);
            std::stringstream temp2(temp_string);
            temp2 >> dummy >> dummy >> dummy >> dummy
                  >> dummy >> dummy >> dummy >> dummy
                  >> temp_mass >> urqmd_pid >> urqmd_iso3 >> urqmd_charge
                  >> dummy >> dummy >> parent_proc_type;

            particle_info temp_particle_info;
            temp2 >> temp_particle_info.t
                  >> temp_particle_info.x
                  >> temp_particle_info.y
                  >> temp_particle_info.z
                  >> temp_particle_info.E
                  >> temp_particle_info.px
                  >> temp_particle_info.py
                  >> temp_particle_info.pz;
            temp_particle_info.mass = temp_mass;
            temp_particle_info.monval = get_pdg_id(urqmd_pid, urqmd_iso3);
            (*full_particle_list)[ievent]->push_back(temp_particle_info);
        }
        num_particles += n_particle;
        ievent++;
    }
    return(0);
}


int particleSamples::read_in_particle_samples_mixed_event_Sangwook() {
    // clean out the previous record
    clear_out_previous_record(full_particle_list_mixed_event);

    std::string temp_string;
    int n_particle;
    double dummy;
    int parent_proc_type;
    int ievent = 0;
    int urqmd_pid, urqmd_iso3, urqmd_charge;
    double temp_mass;
    int num_particles = 0;
    while (num_particles < event_buffer_size) {
        getline(inputfile_mixed_event, temp_string);

        if (inputfile_mixed_event.eof()) break;

        full_particle_list_mixed_event->push_back(
                                            new vector<particle_info> );

        // get number of particles within the event
        std::stringstream temp1(temp_string);
        temp1 >> n_particle;
        // then get one useless line
        getline(inputfile_mixed_event, temp_string);

        for (int ipart = 0; ipart < n_particle; ipart++) {
            getline(inputfile_mixed_event, temp_string);
            std::stringstream temp2(temp_string);
            temp2 >> dummy >> dummy >> dummy >> dummy
                  >> dummy >> dummy >> dummy >> dummy
                  >> temp_mass >> urqmd_pid >> urqmd_iso3 >> urqmd_charge
                  >> dummy >> dummy >> parent_proc_type;

            particle_info temp_particle_info;
            temp2 >> temp_particle_info.t
                  >> temp_particle_info.x
                  >> temp_particle_info.y
                  >> temp_particle_info.z
                  >> temp_particle_info.E
                  >> temp_particle_info.px
                  >> temp_particle_info.py
                  >> temp_particle_info.pz;
            temp_particle_info.mass = temp_mass;
            temp_particle_info.monval = get_pdg_id(urqmd_pid, urqmd_iso3);
            (*full_particle_list_mixed_event)[ievent]->push_back(
                                                    temp_particle_info);
        }
        num_particles += n_particle;
        ievent++;
    }
    return(0);
}


void particleSamples::perform_resonance_feed_down(
                    vector< vector<particle_info>* >* input_particle_list) {
    messager.info("perform resonance decays... ");
    // loop over events
    int ievent = 0;
    for (auto &ev_i: (*input_particle_list)) {
        // create a temporary particle list
        std::list<particle_info> current_ev_list;

        // copy all particles into the temp list
        for (auto &part_i: (*ev_i)) {
            current_ev_list.push_back(part_i);
        }

        (*input_particle_list)[ievent]->clear();

        // perform resonance decays
        while (current_ev_list.size() > 0) {
            auto part_i = current_ev_list.front();
            current_ev_list.pop_front();
            vector<particle_info> daughter_list;
            decayer_ptr->perform_decays(part_i, daughter_list);
            for (auto &daughter_i: daughter_list) {
                if (decayer_ptr->check_particle_stable(daughter_i) == 1) {
                    (*input_particle_list)[ievent]->push_back(daughter_i);
                } else {
                    current_ev_list.push_back(daughter_i);
                }
            }
            daughter_list.clear();
        }
        current_ev_list.clear();
        ievent++;
    }
}


void particleSamples::perform_weak_resonance_feed_down_Sigma_to_Lambda() {
    if (std::abs(particle_monval) == 3122) {
        // consider Sigma^0 feed down to Lambda
        for (unsigned int iev = 0; iev < resonance_list_Sigma0->size();
             iev++) {
            for (unsigned int i = 0;
                 i < (*resonance_list_Sigma0)[iev]->size(); i++) {
                particle_info daughter1;
                particle_info daughter2;
                daughter1.mass = 1.116;  // mass of Lambda
                daughter2.mass = 0.0;    // mass of photon
                decayer_ptr->perform_two_body_decay(
                        (*(*resonance_list_Sigma0)[iev])[i],
                        daughter1, daughter2);
                (*particle_list)[iev]->push_back(daughter1);
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
                    if (std::abs(invariant_mass - particle_mass)
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
