// Helpers for unit tests that build small read_in_mode = 9 fixtures.
#ifndef UNIT_TESTS_BINARY_EVENT_WRITER_H_
#define UNIT_TESTS_BINARY_EVENT_WRITER_H_

#include <sys/stat.h>

#include <cmath>
#include <fstream>
#include <string>
#include <vector>

namespace test_fixtures {

// One track of a read_in_mode = 9 record. Momenta are in GeV and
// coordinates in fm; E is set on the mass shell.
struct Track {
    int pdg;
    double mass;
    double t, x, y, z;
    double px, py, pz;
    double E() const {
        return std::sqrt(mass * mass + px * px + py * py + pz * pz);
    }
};

// Track with transverse momentum pT, momentum rapidity y and azimuth phi.
inline Track track_from_rapidity(
    int pdg, double mass, double pT, double y, double phi) {
    const double mT = std::sqrt(mass * mass + pT * pT);
    return Track {
        pdg,
        mass,
        0.0,
        0.0,
        0.0,
        0.0,
        pT * std::cos(phi),
        pT * std::sin(phi),
        mT * std::sinh(y)};
}

// Writes events in the read_in_mode = 9 layout: a native int track count,
// then per track an int PDG code and float m, t, x, y, z, E, px, py, pz.
// Creates the directory and writes <directory>/particle_list.bin.
inline void write_binary_events(
    const std::string &directory,
    const std::vector<std::vector<Track>> &events) {
    mkdir(directory.c_str(), 0755);
    std::ofstream out(
        (directory + "/particle_list.bin").c_str(), std::ios::binary);
    for (const auto &event : events) {
        const int count = static_cast<int>(event.size());
        out.write(reinterpret_cast<const char *>(&count), sizeof(int));
        for (const auto &track : event) {
            const float record[9] = {
                static_cast<float>(track.mass), static_cast<float>(track.t),
                static_cast<float>(track.x),    static_cast<float>(track.y),
                static_cast<float>(track.z),    static_cast<float>(track.E()),
                static_cast<float>(track.px),   static_cast<float>(track.py),
                static_cast<float>(track.pz)};
            out.write(reinterpret_cast<const char *>(&track.pdg), sizeof(int));
            out.write(reinterpret_cast<const char *>(record), sizeof(record));
        }
    }
}

}  // namespace test_fixtures

#endif  // UNIT_TESTS_BINARY_EVENT_WRITER_H_
