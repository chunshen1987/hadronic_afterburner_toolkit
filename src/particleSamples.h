// Copyright Chun Shen @ 2016
#ifndef SRC_particleSamples_h_
#define SRC_particleSamples_h_

#include <fstream>
#include <vector>
#include <string>
#include <map>
#include <utility>

#include "zlib.h"

#include "ParameterReader.h"
#include "particle_info.h"
#include "particle_decay.h"
#include "Random.h"
#include "pretty_ostream.h"

class particleSamples {
 private:
    const ParameterReader paraRdr_;
    pretty_ostream messager;

    //! path for results folder
    const std::string path_;

    std::ifstream inputfile;
    std::ifstream inputfile_mixed_event;
    uint16_t smash_format_version_;
    gzFile inputfile_gz;
    gzFile inputfile_mixed_event_gz;
    int event_buffer_size;

    int echo_level_;
    int read_in_mode_;
    bool read_mixed_events;

    bool analyze_flow;
    bool analyze_HBT;
    bool analyze_BF;
    bool analyze_ebedis;

    int rap_type_;

    //! the monte-carlo number of the particle of interest
    int particle_monval;

    //! perform resonance feed down for all unstable particles
    int resonance_feed_down_flag;
    int resonance_weak_feed_down_flag;
    //! perform resonance decays for only selected particle species
    int select_resonances_flag;
    //! store the monval of the selected resonance states
    std::vector<int> select_resonances_list;

    //! include Sigma0 feed down of Lambda
    int resonance_weak_feed_down_Sigma_to_Lambda_flag;

    //! reconst phi meson from (K^+, K^-) pairs
    int reconst_flag;

    //! flag to collect net particle distribution (for run_mode == 2)
    int net_particle_flag;

    //! flag to distinguish particle's isospin
    int flag_isospin;

    //! flag to collect positive and negative charged hadron seperately
    int flag_charge_dependence;

    //! particle monte-carlo number of balance function
    int particle_monval_a;
    int particle_monval_abar;
    int particle_monval_b;
    int particle_monval_bbar;

    std::map<std::pair<int, int>, int> urqmd_to_pdg;

    //! this list store all particle samples
    std::vector< std::vector<particle_info>* >* full_particle_list;
    std::vector< std::vector<particle_info>* >* full_particle_list_mixed_event;

    //! particle list to store the select particle sample
    std::vector< std::vector<particle_info>* >* particle_list;

    //! particle list to store anti-particles
    //! (used when net_particle_flag == 1)
    std::vector< std::vector<particle_info>* >* anti_particle_list;

    //! particle list to store the selected particle sample from a mix event
    //! (used when run_mode == 1 for HBT calculation)
    std::vector< std::vector<particle_info>* >* particle_list_mixed_event;

    std::vector< std::vector<particle_info>* >* resonance_list;

    //! particle list to store the resonance particles (Sigma0)
    //! (used when resonance_weak_feed_down_Sigma_to_Lambda_flag == 1)
    std::vector< std::vector<particle_info>* >* resonance_list_Sigma0;

    //! particle list to store the (K^+ and K^-) pairs
    //! (used when reconst_flag == 1)
    std::vector< std::vector<particle_info>* >* reconst_list_1;
    std::vector< std::vector<particle_info>* >* reconst_list_2;

    //! particle list to store the positive hadrons
    //! (used when flag_charge_dependence == 1)
    std::vector< std::vector<particle_info>* >* positive_charge_hadron_list;
    //! particle list to store the negative hadrons
    //! (used when flag_charge_dependence == 1)
    std::vector< std::vector<particle_info>* >* negative_charge_hadron_list;

    //! particle list to store the select particle samples for balance function
    std::vector< std::vector<particle_info>* >* balance_function_particle_a;
    std::vector< std::vector<particle_info>* >* balance_function_particle_b;
    std::vector< std::vector<particle_info>* >* balance_function_particle_abar;
    std::vector< std::vector<particle_info>* >* balance_function_particle_bbar;

    std::vector< std::vector<particle_info>* >* balance_function_particle_a_mixed_event;
    std::vector< std::vector<particle_info>* >* balance_function_particle_b_mixed_event;
    std::vector< std::vector<particle_info>* >* balance_function_particle_abar_mixed_event;
    std::vector< std::vector<particle_info>* >* balance_function_particle_bbar_mixed_event;

    //! particle decay
    particle_decay *decayer_ptr;

 public:
    particleSamples(ParameterReader &paraRdr, std::string path,
                    std::shared_ptr<RandomUtil::Random> ran_gen);
    ~particleSamples();

    int get_pdg_id(int urqmd_id, int urqmd_isospin);
    void build_map_urqmd_to_pdg_id();

    int decide_to_pick_resonance(int monval);
    int decide_to_pick_reconst(int monval);
    int decide_to_pick_anti_particles(int monval);
    int decide_to_pick_charge(int monval);
    bool decide_to_pick_OSCAR(int POI, int monval);

    void initialize_selected_resonance_list();
    void perform_resonance_feed_down(
            std::vector< std::vector<particle_info>* >* input_particle_list);
    void perform_weak_resonance_feed_down_Sigma_to_Lambda();
    void perform_particle_reconstruction();

    int read_in_particle_samples();
    void read_in_particle_samples_and_filter();
    int read_in_particle_samples_mixed_event();
    void read_in_particle_samples_mixed_event_and_filter();

    int read_in_particle_samples_OSCAR();
    int read_in_particle_samples_OSCAR_mixed_event();
    int read_in_particle_samples_JAM();
    int read_in_particle_samples_JAM_mixed_event();
    int read_in_particle_samples_UrQMD();
    int read_in_particle_samples_UrQMD_mixed_event();
    int read_in_particle_samples_UrQMD_zipped();
    int read_in_particle_samples_UrQMD_mixed_event_zipped();
    int read_in_particle_samples_UrQMD_3p3();
    int read_in_particle_samples_UrQMD_3p3_mixed_event();
    int read_in_particle_samples_Sangwook();
    int read_in_particle_samples_SMASH_binary();
    int read_in_particle_samples_mixed_event_Sangwook();
    int read_in_particle_samples_SMASH_gzipped();
    int read_in_particle_samples_SMASH_mixed_event_gzipped();
    int read_in_particle_samples_gzipped();
    int read_in_particle_samples_mixed_event_gzipped();
    void clear_out_previous_record(
                    std::vector< std::vector<particle_info>* >* plist);

    void filter_particles_from_events(const int PoI_monval);

    void filter_particles(const int PoI_monval,
            std::vector< std::vector<particle_info>* >* full_list,
            std::vector< std::vector<particle_info>* >* filted_list);
    void filter_particles_into_lists(
                    std::vector< std::vector<particle_info>* >* full_list);

    std::string gz_readline(gzFile gzfp);
    bool end_of_file() const {
        if (read_in_mode_ == 2 || read_in_mode_ == 7 || read_in_mode_ == 10) {
            return(gzeof(inputfile_gz));
        } else {
            return(inputfile.eof());
        }
    }
    bool end_of_file_mixed_event() {
        if (read_in_mode_ == 2 || read_in_mode_ == 7 || read_in_mode_ == 10) {
            return(gzeof(inputfile_mixed_event_gz));
        } else {
            return(inputfile_mixed_event.eof());
        }
    }

    int get_event_buffer_size() const {return(event_buffer_size);}

    int get_number_of_events() const {return(particle_list->size());}
    int get_number_of_events_anti_particle() {
        return(anti_particle_list->size());
    }
    int get_number_of_mixed_events() {
        return(particle_list_mixed_event->size());
    }

    int get_number_of_particles(int event_id) {
        return((*particle_list)[event_id]->size());
    }
    int get_number_of_particles_mixed_event(int event_id) {
        return((*particle_list_mixed_event)[event_id]->size());
    }
    int get_number_of_anti_particles(int event_id) {
        return((*anti_particle_list)[event_id]->size());
    }
    int get_number_of_positive_particles(int event_id) {
        return((*positive_charge_hadron_list)[event_id]->size());
    }
    int get_number_of_negative_particles(int event_id) {
        return((*negative_charge_hadron_list)[event_id]->size());
    }

    particle_info get_particle(int event_id, int part_id) {
        return((*(*particle_list)[event_id])[part_id]);
    }
    particle_info get_anti_particle(int event_id, int part_id) {
        return((*(*anti_particle_list)[event_id])[part_id]);
    }
    particle_info get_particle_from_mixed_event(int event_id, int part_id) {
        return((*(*particle_list_mixed_event)[event_id])[part_id]);
    }
    particle_info get_positive_particle(int event_id, int part_id) {
        return((*(*positive_charge_hadron_list)[event_id])[part_id]);
    }
    particle_info get_negative_particle(int event_id, int part_id) {
        return((*(*negative_charge_hadron_list)[event_id])[part_id]);
    }

    std::vector< std::vector<particle_info>* >* get_balance_function_particle_list_a() const {
        return(balance_function_particle_a);
    }

    std::vector< std::vector<particle_info>* >* get_balance_function_particle_list_abar() const {
        return(balance_function_particle_abar);
    }

    std::vector< std::vector<particle_info>* >* get_balance_function_particle_list_b() const {
        return(balance_function_particle_b);
    }

    std::vector< std::vector<particle_info>* >* get_balance_function_particle_list_bbar() const {
        return(balance_function_particle_bbar);
    }

    std::vector< std::vector<particle_info>* >* get_balance_function_particle_list_a_mixed_event() const {
        return(balance_function_particle_a_mixed_event);
    }

    std::vector< std::vector<particle_info>* >* get_balance_function_particle_list_abar_mixed_event() const {
        return(balance_function_particle_abar_mixed_event);
    }

    std::vector< std::vector<particle_info>* >* get_balance_function_particle_list_b_mixed_event() const {
        return(balance_function_particle_b_mixed_event);
    }

    std::vector< std::vector<particle_info>* >* get_balance_function_particle_list_bbar_mixed_event() const {
        return(balance_function_particle_bbar_mixed_event);
    }
};

#endif  // SRC_particleSamples_h_
