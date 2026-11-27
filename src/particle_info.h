// Copyright Chun Shen @ 2016
#ifndef SRC_particle_info_h_
#define SRC_particle_info_h_

#include <cstddef>

typedef struct {
    int monval;
    double mass;
    int baryon;
    double px, py, pz, E;
    double x, y, z, t;
    double rap_y, rap_eta, pT, phi_p;
    std::size_t event_track_id = 0;
    // true for a parent rebuilt from its decay daughters (phi -> K+ K-)
    bool reconstructed = false;
} particle_info;

#endif  // SRC_particle_info_h_
