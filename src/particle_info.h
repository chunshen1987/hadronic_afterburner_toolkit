// Copyright Chun Shen @ 2016
#ifndef SRC_particle_info_h_
#define SRC_particle_info_h_

typedef struct {
    int monval;
    double mass;
    double px, py, pz, E;
    double x, y, z, t;
    double rap_y, pT, phi_p;
} particle_info;

#endif  // SRC_particle_info_h_
