// Copyright Chun Shen @ 2016

#include "particle_decay.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>

#include "arsenal.h"

using std::cout;
using std::endl;

namespace {
constexpr double mass_margin = 1e-5;  // GeV
constexpr double broad_width = 1e-4;  // GeV
constexpr int max_daughters = 5;
}  // namespace

particle_decay::particle_decay(
    std::shared_ptr<RandomUtil::Random> ran_gen, int weak_flag) {
    ran_gen_ptr_ = ran_gen;
    weak_flag_ = weak_flag;

    read_resonances_list();
}

particle_decay::~particle_decay() {
    const long unreachable =
        unreachable_resonances_.load(std::memory_order_relaxed);
    const long fallback =
        mass_sampling_fallbacks_.load(std::memory_order_relaxed);
    if (unreachable > 0) {
        std::cerr << "[Warning] particle_decay: kept " << unreachable
                  << " resonances below all of their decay thresholds" << endl;
    }
    if (fallback > 0) {
        std::cerr
            << "[Warning] particle_decay: used fallback daughter masses in "
            << fallback << " decays" << endl;
    }
}

int particle_decay::read_resonances_list() {
    cout << " -- Read in particle resonance decay table..." << endl;
    const std::string filename =
        weak_flag_ == 0 ? "EOS/pdg.dat" : "EOS/pdg_weak.dat";
    std::ifstream table(filename);
    if (!table) throw std::runtime_error("Cannot open " + filename);

    resonance_table_.clear();
    std::vector<int> baryons;
    int monval;
    while (table >> monval) {
        particle_decay_info particle {};
        particle.monval = monval;
        if (!(table >> particle.name >> particle.mass >> particle.width
              >> particle.gspin >> particle.baryon >> particle.strange
              >> particle.charm >> particle.bottom >> particle.gisospin
              >> particle.charge >> particle.decays)
            || particle.decays < 0) {
            throw std::runtime_error("Invalid resonance entry in " + filename);
        }

        for (int i = 0; i < particle.decays; i++) {
            decay_channel_info channel {};
            int parent;
            if (!(table >> parent >> channel.decay_Npart
                  >> channel.branching_ratio >> channel.decay_part[0]
                  >> channel.decay_part[1] >> channel.decay_part[2]
                  >> channel.decay_part[3] >> channel.decay_part[4])
                || channel.decay_Npart < -max_daughters
                || channel.decay_Npart > max_daughters) {
                throw std::runtime_error(
                    "Invalid decay channel for PDG " + std::to_string(monval));
            }
            // The UrQMD-derived tables use a negative count for channels
            // that are closed at the daughters' pole masses.
            channel.decay_Npart = std::abs(channel.decay_Npart);
            particle.decay_channels.push_back(channel);
        }
        particle.stable = particle.decay_channels.empty()
                          || particle.decay_channels.front().decay_Npart == 1;
        particle.sign = particle.baryon == 0 ? -1 : 1;
        resonance_table_[monval] = particle;
        if (particle.baryon == 1) baryons.push_back(monval);
    }
    if (!table.eof() || resonance_table_.empty()) {
        throw std::runtime_error("Invalid resonance table " + filename);
    }

    // Build missing antibaryons after reading all particle quantum numbers.
    for (int id : baryons) {
        if (resonance_table_.count(-id)) continue;
        auto anti = particle_data(id);
        anti.monval = -id;
        anti.name = "Anti-" + anti.name;
        anti.baryon = -anti.baryon;
        anti.strange = -anti.strange;
        anti.charm = -anti.charm;
        anti.bottom = -anti.bottom;
        anti.charge = -anti.charge;
        for (auto &channel : anti.decay_channels) {
            for (int i = 0; i < channel.decay_Npart; i++) {
                int &daughter = channel.decay_part[i];
                const int lookup =
                    resonance_table_.count(daughter) ? daughter : -daughter;
                const auto &info = particle_data(lookup);
                if (info.baryon || info.strange || info.charge) {
                    daughter = -daughter;
                }
            }
        }
        resonance_table_.emplace(-id, anti);
    }

    prepare_channel_kinematics();
    return resonance_table_.size();
}

const particle_decay_info &particle_decay::particle_data(int monval) const {
    const auto found = resonance_table_.find(monval);
    if (found == resonance_table_.end()) {
        throw std::runtime_error(
            "PDG " + std::to_string(monval)
            + " is missing from the decay table");
    }
    return found->second;
}

bool particle_decay::is_broad(const particle_decay_info &particle) {
    return !particle.stable && particle.width > broad_width;
}

double particle_decay::compute_minimum_mass(
    int monval, std::map<int, int> &state) {
    const auto it = resonance_table_.find(monval);
    if (it == resonance_table_.end()) {
        throw std::runtime_error(
            "PDG " + std::to_string(monval)
            + " is missing from the decay table");
    }
    auto &particle = it->second;
    if (state[monval] == 2) return particle.minimum_mass;
    if (state[monval] == 1) return particle.mass;  // A cycle in the table.
    state[monval] = 1;
    particle.minimum_mass = particle.mass;
    if (is_broad(particle)) {
        double threshold = 0.0;
        for (const auto &channel : particle.decay_channels) {
            if (channel.branching_ratio <= 0.0 || channel.decay_Npart < 2)
                continue;
            double daughter_sum = 0.0;
            for (int i = 0; i < channel.decay_Npart; i++) {
                daughter_sum +=
                    compute_minimum_mass(channel.decay_part[i], state);
            }
            threshold = std::max(threshold, daughter_sum);
        }
        // One margin for reachability and one for sampling the daughters.
        particle.minimum_mass =
            std::min(particle.mass, threshold + 2 * mass_margin);
    }
    state[monval] = 2;
    return particle.minimum_mass;
}

void particle_decay::prepare_channel_kinematics() {
    std::map<int, int> state;
    for (const auto &entry : resonance_table_) {
        compute_minimum_mass(entry.first, state);
    }
    for (auto &entry : resonance_table_) {
        auto &particle = entry.second;
        if (particle.stable) continue;
        for (auto &channel : particle.decay_channels) {
            if (channel.branching_ratio <= 0.0 || channel.decay_Npart < 2)
                continue;
            channel.pole_threshold = 0.0;
            channel.minimum_threshold = 0.0;
            channel.has_broad_daughter = false;
            for (int i = 0; i < channel.decay_Npart; i++) {
                const auto &daughter = particle_data(channel.decay_part[i]);
                channel.pole_threshold += daughter.mass;
                channel.minimum_threshold += daughter.minimum_mass;
                channel.has_broad_daughter |= is_broad(daughter);
            }
        }
        if (!has_reachable_channel(particle, particle.mass)) {
            std::cerr << "[Warning] particle_decay: PDG " << particle.monval
                      << " has no channel open at its pole mass" << endl;
        } else if (!has_reachable_channel(particle, particle.minimum_mass)) {
            std::cerr << "[Warning] particle_decay: PDG " << particle.monval
                      << " has no channel open at its minimum daughter mass"
                      << endl;
        }
    }
}

bool particle_decay::channel_reachable(
    const decay_channel_info &channel, double mass) {
    return channel.decay_Npart >= 2 && channel.branching_ratio > 0.0
           && (mass > channel.pole_threshold
               || (channel.has_broad_daughter
                   && mass - mass_margin > channel.minimum_threshold));
}

bool particle_decay::has_reachable_channel(
    const particle_decay_info &particle, double mass) {
    return std::any_of(
        particle.decay_channels.begin(), particle.decay_channels.end(),
        [mass](const decay_channel_info &channel) {
            return channel_reachable(channel, mass);
        });
}

//! This is a test function to check whether the resonance table is read in
//! correctly
void particle_decay::check_resonance_table() const {
    for (auto itr_reso = resonance_table_.begin();
         itr_reso != resonance_table_.end(); itr_reso++) {
        auto ireso = itr_reso->second;
        cout << "name: " << ireso.name << endl;
        cout << "monval: " << ireso.monval << endl;
        cout << "mass: " << ireso.mass << endl;
        cout << "stable: " << ireso.stable << endl;
        cout << "# of decay channels: " << ireso.decays << endl;
        for (int j = 0; j < ireso.decays; j++) {
            cout << j
                 << ": # of daughters: " << ireso.decay_channels[j].decay_Npart
                 << endl;
            cout << j << ": branching_ratio: "
                 << ireso.decay_channels[j].branching_ratio << endl;
            cout << j << ": decay part: ";
            for (int k = 0; k < ireso.decay_channels[j].decay_Npart; k++) {
                cout << ireso.decay_channels[j].decay_part[k] << "  ";
            }
            cout << endl;
        }
    }
}

//! This function returns particle width in GeV
double particle_decay::get_particle_width(const particle_info part) const {
    double width = particle_data(part.monval).width;
    return (width);
}

//! This function checks whether the particle is stable
int particle_decay::check_particle_stable(const particle_info part) const {
    const auto &particle = particle_data(part.monval);
    if (particle.stable) return 1;
    if (has_reachable_channel(particle, part.mass)) return 0;
    unreachable_resonances_.fetch_add(1, std::memory_order_relaxed);
    return 1;
}

//! This function returns the electric charge of particle
int particle_decay::get_particle_charge(const int monval) const {
    int charge = particle_data(monval).charge;
    return (charge);
}

//! This function returns the baryon number of particle
int particle_decay::get_particle_baryon_number(const int monval) const {
    int baryon = particle_data(monval).baryon;
    return (baryon);
}

//! This function returns the strange number of particle
int particle_decay::get_particle_strange_number(const int monval) const {
    int strange = particle_data(monval).strange;
    return (strange);
}

void particle_decay::sample_daughter_masses(
    double mass, std::vector<particle_info> &daughters,
    RandomUtil::Random &ran) {
    const int count = daughters.size();
    const double limit = mass - mass_margin;
    std::array<const particle_decay_info *, max_daughters> data {};
    std::array<double, max_daughters> lower {}, angle_min {}, angle_span {};
    std::array<double, max_daughters> upper {};
    double lower_sum = 0.0;
    int broad_count = 0;
    for (int i = 0; i < count; i++) {
        const auto &info = particle_data(daughters[i].monval);
        data[i] = &info;
        lower[i] = info.minimum_mass;
        lower_sum += lower[i];
        if (is_broad(info)) broad_count++;
    }
    for (int i = 0; i < count; i++) {
        const auto &info = *data[i];
        if (!is_broad(info)) continue;
        upper[i] = limit - (lower_sum - lower[i]);
        const double half_width = 0.5 * info.width;
        angle_min[i] = std::atan2(lower[i] - info.mass, half_width);
        angle_span[i] =
            std::atan2(upper[i] - info.mass, half_width) - angle_min[i];
    }

    // Rejection sampling preserves the product of the truncated Breit-Wigner
    // densities, conditioned on the daughter masses fitting in this decay.
    for (int attempt = 0; attempt < 100000; attempt++) {
        double sum = 0.0;
        for (int i = 0; i < count; i++) {
            const auto &info = *data[i];
            if (is_broad(info)) {
                const double angle =
                    angle_min[i] + angle_span[i] * ran.rand_uniform();
                const double sampled =
                    info.mass + 0.5 * info.width * std::tan(angle);
                daughters[i].mass =
                    std::min(std::max(sampled, lower[i]), upper[i]);
            }
            sum += daughters[i].mass;
        }
        if (sum < limit) return;
    }

    // Extremely suppressed tails still give a kinematically open decay.
    mass_sampling_fallbacks_.fetch_add(1, std::memory_order_relaxed);
    const double share = 0.5 * (limit - lower_sum) / broad_count;
    for (int i = 0; i < count; i++) {
        daughters[i].mass = lower[i] + (is_broad(*data[i]) ? share : 0.0);
    }
}

void particle_decay::perform_decays(
    particle_info &mother, std::vector<particle_info> &daughter_list) {
    perform_decays(mother, daughter_list, *ran_gen_ptr_);
}

void particle_decay::perform_decays(
    particle_info &mother, std::vector<particle_info> &daughter_list,
    RandomUtil::Random &ran) {
    const auto &particle = particle_data(mother.monval);
    double total = 0.0;
    if (!particle.stable) {
        for (const auto &channel : particle.decay_channels) {
            if (channel_reachable(channel, mother.mass)) {
                total += channel.branching_ratio;
            }
        }
    }
    if (total <= 0.0) {
        // check_particle_stable recognizes these as terminal in a cascade.
        daughter_list.push_back(mother);
        return;
    }

    double target = ran.rand_uniform();
    // Retain the original channel draw when the reachable weights sum to one.
    if (std::abs(total - 1.0) > 1e-12) target *= total;
    const decay_channel_info *picked = nullptr;
    double cumulative = 0.0;
    for (const auto &channel : particle.decay_channels) {
        if (!channel_reachable(channel, mother.mass)) continue;
        picked = &channel;
        cumulative += channel.branching_ratio;
        if (cumulative > target) break;
    }

    const int count = picked->decay_Npart;
    static thread_local std::vector<particle_info> daughters;
    daughters.assign(count, particle_info {});
    for (int i = 0; i < count; i++) {
        const auto &info = particle_data(picked->decay_part[i]);
        daughters[i].monval = info.monval;
        daughters[i].mass = info.mass;
        daughters[i].baryon = info.baryon;
    }
    if (mother.mass <= picked->pole_threshold) {
        sample_daughter_masses(mother.mass, daughters, ran);
    }
    if (count == 2) {
        perform_two_body_decay(mother, daughters[0], daughters[1], ran);
    } else if (count == 3) {
        perform_three_body_decay(
            mother, daughters[0], daughters[1], daughters[2], ran);
    } else {
        perform_n_body_decay(mother, daughters, ran);
    }
    daughter_list.insert(
        daughter_list.end(), daughters.begin(), daughters.end());
}

namespace {
long double breakup_momentum(
    long double mass, long double first, long double second) {
    const long double sum = first + second;
    if (mass <= sum) return 0.0;
    const long double difference = first - second;
    return std::sqrt(
               (mass - sum) * (mass + sum) * (mass - difference)
               * (mass + difference))
           / (2 * mass);
}

struct DecayFrame {
    double E, px, py, pz, mass;
};

DecayFrame decay_frame(const particle_info &mother, long double threshold) {
    const long double E = mother.E;
    const long double px = mother.px, py = mother.py, pz = mother.pz;
    const long double momentum_squared = px * px + py * py + pz * pz;
    const long double mass_squared = E * E - momentum_squared;
    if (E > 0.0 && mass_squared >= threshold * threshold) {
        const double mass = std::sqrt(mass_squared);
        if (mass > 0.0 && mass >= threshold) {
            return {mother.E, mother.px, mother.py, mother.pz, mass};
        }
    }

    // A channel near threshold can close when a single-precision file rounds
    // the input four-vector. Repair only rounding-sized energy discrepancies.
    const double on_shell_energy = std::hypot(
        mother.mass, static_cast<double>(std::sqrt(momentum_squared)));
    const double tolerance =
        8 * std::numeric_limits<float>::epsilon() * on_shell_energy;
    if (mother.mass > 0.0 && mother.mass >= threshold && mother.E > 0.0
        && std::abs(mother.E - on_shell_energy) <= tolerance) {
        return {on_shell_energy, mother.px, mother.py, mother.pz, mother.mass};
    }
    throw std::runtime_error(
        "Stored four-momentum cannot support decay of PDG "
        + std::to_string(mother.monval));
}

void boost_from_rest_frame(const DecayFrame &frame, particle_info &daughter) {
    const double dot = frame.px * daughter.px + frame.py * daughter.py
                       + frame.pz * daughter.pz;
    const double coefficient =
        (daughter.E + dot / (frame.E + frame.mass)) / frame.mass;
    daughter.E = (frame.E * daughter.E + dot) / frame.mass;
    daughter.px += coefficient * frame.px;
    daughter.py += coefficient * frame.py;
    daughter.pz += coefficient * frame.pz;
}
}  // namespace

// Raubold-Lynch phase space, weighted by the product of the sequential
// two-body breakup momenta. Each accepted step has an isotropic direction.
void particle_decay::perform_n_body_decay(
    const particle_info &mother, std::vector<particle_info> &daughters,
    RandomUtil::Random &ran) {
    const int count = daughters.size();
    std::array<long double, max_daughters> mass_sum {}, invariant {},
        momentum {};
    mass_sum[0] = daughters[0].mass;
    for (int i = 1; i < count; i++) {
        mass_sum[i] = mass_sum[i - 1] + daughters[i].mass;
    }
    const long double threshold = mass_sum[count - 1];
    if (mother.mass <= threshold) {
        throw std::runtime_error(
            "Closed phase space for PDG " + std::to_string(mother.monval));
    }
    const auto frame = decay_frame(mother, threshold);
    const long double kinetic = frame.mass - threshold;
    long double max_weight = 1.0;
    for (int i = 1; i < count; i++) {
        max_weight *= breakup_momentum(
            kinetic + mass_sum[i], mass_sum[i - 1], daughters[i].mass);
    }

    static thread_local std::vector<double> fraction;
    fraction.resize(count);
    fraction[count - 1] = 1.0;
    while (true) {
        for (int i = 1; i < count - 1; i++) {
            fraction[i] = ran.rand_uniform();
        }
        std::sort(fraction.begin() + 1, fraction.begin() + count - 1);
        invariant[0] = mass_sum[0];
        long double weight = 1.0;
        for (int i = 1; i < count; i++) {
            invariant[i] = mass_sum[i] + fraction[i] * kinetic;
            momentum[i] = breakup_momentum(
                invariant[i], invariant[i - 1], daughters[i].mass);
            weight *= momentum[i];
        }
        if (ran.rand_uniform() * max_weight < weight) break;
    }

    for (int i = 1; i < count; i++) {
        const double cosine = 2.0 * ran.rand_uniform() - 1.0;
        const double sine = std::sqrt(1.0 - cosine * cosine);
        const double phi = 2.0 * M_PI * ran.rand_uniform();
        const double p = momentum[i];
        auto &daughter = daughters[i];
        daughter.px = p * sine * std::cos(phi);
        daughter.py = p * sine * std::sin(phi);
        daughter.pz = p * cosine;
        daughter.E = std::hypot(p, daughter.mass);
        const double subsystem_mass = invariant[i - 1];
        const double subsystem_energy = std::hypot(p, subsystem_mass);
        if (i == 1) {
            daughters[0].E = subsystem_energy;
            daughters[0].px = -daughter.px;
            daughters[0].py = -daughter.py;
            daughters[0].pz = -daughter.pz;
        } else {
            const DecayFrame subsystem {
                subsystem_energy, -daughter.px, -daughter.py, -daughter.pz,
                subsystem_mass};
            for (int j = 0; j < i; j++) {
                boost_from_rest_frame(subsystem, daughters[j]);
            }
        }
    }

    const double width = get_particle_width(mother);
    double lifetime = 1e10;
    if (width > 1e-10) {
        lifetime = -std::log(ran.rand_uniform()) * frame.E
                   / (frame.mass * width) * AfterburnerUtil::hbarc;
    }
    const double vx = frame.px / frame.E;
    const double vy = frame.py / frame.E;
    const double vz = frame.pz / frame.E;
    for (auto &daughter : daughters) {
        boost_from_rest_frame(frame, daughter);
        daughter.t = mother.t + lifetime;
        daughter.x = mother.x + vx * lifetime;
        daughter.y = mother.y + vy * lifetime;
        daughter.z = mother.z + vz * lifetime;
    }
}

//! This function returns the particle mass for a given particle id
double particle_decay::get_particle_mass(const int POI_monval) const {
    double mass = particle_data(POI_monval).mass;
    return (mass);
}

//! This function perform two body decay
void particle_decay::perform_two_body_decay(
    particle_info &mother, particle_info &daughter1, particle_info &daughter2) {
    perform_two_body_decay(mother, daughter1, daughter2, *ran_gen_ptr_);
}

void particle_decay::perform_two_body_decay(
    particle_info &mother, particle_info &daughter1, particle_info &daughter2,
    RandomUtil::Random &ran) {
    double M_pole = mother.mass;
    double M_width = get_particle_width(mother);
    double m1 = daughter1.mass;
    double m2 = daughter2.mass;
    double M_min = m1 + m2;
    if (M_pole < M_min) {
        throw std::runtime_error(
            "Listed mass is below the decay threshold for PDG "
            + std::to_string(mother.monval));
    }
    const auto frame = decay_frame(mother, M_min);
    const double p_lrf = breakup_momentum(frame.mass, m1, m2);

    // randomly pick emission angle
    double phi = ran.rand_uniform() * 2 * M_PI;
    double cos_theta = 2. * (ran.rand_uniform() - 0.5);
    double sin_theta = sqrt(1. - cos_theta * cos_theta);

    // compute daughter particles' energy and momentum in the rest frame
    double E1_lrf = sqrt(p_lrf * p_lrf + m1 * m1);
    double p1_lrf_x = p_lrf * sin_theta * cos(phi);
    double p1_lrf_y = p_lrf * sin_theta * sin(phi);
    double p1_lrf_z = p_lrf * cos_theta;
    double E2_lrf = sqrt(p_lrf * p_lrf + m2 * m2);
    double p2_lrf_x = -p1_lrf_x;
    double p2_lrf_y = -p1_lrf_y;
    double p2_lrf_z = -p1_lrf_z;

    // compute mother velocity
    double vx = frame.px / frame.E;
    double vy = frame.py / frame.E;
    double vz = frame.pz / frame.E;

    daughter1.E = E1_lrf;
    daughter1.px = p1_lrf_x;
    daughter1.py = p1_lrf_y;
    daughter1.pz = p1_lrf_z;
    daughter2.E = E2_lrf;
    daughter2.px = p2_lrf_x;
    daughter2.py = p2_lrf_y;
    daughter2.pz = p2_lrf_z;
    boost_from_rest_frame(frame, daughter1);
    boost_from_rest_frame(frame, daughter2);

    double life_time = 1e10;
    if (M_width > 1e-10) {
        // compute life-time = gamma*1/\Gamma
        double tau0 = frame.E / frame.mass * 1. / (M_width);
        life_time = -tau0 * log(ran.rand_uniform());
        life_time *= AfterburnerUtil::hbarc;  // convert to fm
    }

    daughter1.t = mother.t + life_time;
    daughter1.x = mother.x + vx * life_time;
    daughter1.y = mother.y + vy * life_time;
    daughter1.z = mother.z + vz * life_time;
    daughter2.t = mother.t + life_time;
    daughter2.x = mother.x + vx * life_time;
    daughter2.y = mother.y + vy * life_time;
    daughter2.z = mother.z + vz * life_time;
}

//! This function perform 3 body decays
void particle_decay::perform_three_body_decay(
    particle_info &mother, particle_info &daughter1, particle_info &daughter2,
    particle_info &daughter3) {
    perform_three_body_decay(
        mother, daughter1, daughter2, daughter3, *ran_gen_ptr_);
}

void particle_decay::perform_three_body_decay(
    particle_info &mother, particle_info &daughter1, particle_info &daughter2,
    particle_info &daughter3, RandomUtil::Random &ran) {
    double M_pole = mother.mass;
    double M_width = get_particle_width(mother);
    double m1 = daughter1.mass;
    double m2 = daughter2.mass;
    double m3 = daughter3.mass;
    double M_min = m1 + m2 + m3;
    if (M_pole < M_min) {
        throw std::runtime_error(
            "Listed mass is below the decay threshold for PDG "
            + std::to_string(mother.monval));
    }
    const auto frame = decay_frame(mother, M_min);
    // Use kinetic energies to avoid subtracting nearly equal squares close
    // to threshold. At equality all daughters are at rest in this frame.
    double p1_lrf = 0.0, p2_lrf = 0.0, cos12_lrf = 1.0;
    const double available = frame.mass - M_min;
    if (available > 0.0) {
        do {
            const double kinetic1 = ran.rand_uniform() * available;
            const double kinetic2 = ran.rand_uniform() * available;
            const double kinetic3 = available - kinetic1 - kinetic2;
            if (kinetic3 < 0.0) {
                cos12_lrf = 2.0;
                continue;
            }
            const double p1_sq = kinetic1 * (kinetic1 + 2 * m1);
            const double p2_sq = kinetic2 * (kinetic2 + 2 * m2);
            const double p3_sq = kinetic3 * (kinetic3 + 2 * m3);
            p1_lrf = std::sqrt(p1_sq);
            p2_lrf = std::sqrt(p2_sq);
            const double divisor = 2 * p1_lrf * p2_lrf;
            cos12_lrf = divisor > 0.0 ? (p3_sq - p1_sq - p2_sq) / divisor : 2.0;
        } while (!(cos12_lrf >= -1.0 && cos12_lrf <= 1.0));
    }

    // sample the lifetime
    double life_time = 1e10;
    if (M_width > 1e-10) {
        double tau = frame.E / frame.mass * 1. / M_width;
        life_time = -tau * log(ran.rand_uniform());
        life_time *= AfterburnerUtil::hbarc;  // convert unit to fm
    }
    // compute the decay position
    double decay_time = mother.t + life_time;
    double decay_x = mother.x + frame.px / frame.E * life_time;
    double decay_y = mother.y + frame.py / frame.E * life_time;
    double decay_z = mother.z + frame.pz / frame.E * life_time;

    // compute the momentum of decay daughters
    double tp2_lrf_x = p2_lrf * sqrt(1. - cos12_lrf * cos12_lrf);
    double tp2_lrf_z = p2_lrf * cos12_lrf;
    double tp3_lrf_x = -tp2_lrf_x;
    double tp3_lrf_z = -(p1_lrf + tp2_lrf_z);
    double phi = 2. * M_PI * ran.rand_uniform();
    double ksi = 2. * M_PI * ran.rand_uniform();
    double cos_theta = 2. * ran.rand_uniform() - 1.0;

    double sin_phi = sin(phi);
    double cos_phi = cos(phi);
    double sin_ksi = sin(ksi);
    double cos_ksi = cos(ksi);
    double sin_theta = sqrt(1. - cos_theta * cos_theta);

    double p1_lrf_x = -p1_lrf * sin_theta * cos_ksi;
    double p1_lrf_y = p1_lrf * sin_theta * sin_ksi;
    double p1_lrf_z = p1_lrf * cos_theta;
    const double E1_lrf = sqrt(
        m1 * m1 + p1_lrf_x * p1_lrf_x + p1_lrf_y * p1_lrf_y
        + p1_lrf_z * p1_lrf_z);
    double p2_lrf_x =
        (tp2_lrf_x * (cos_phi * cos_theta * cos_ksi - sin_phi * sin_ksi)
         - tp2_lrf_z * sin_theta * cos_ksi);
    double p2_lrf_y =
        (tp2_lrf_x * (-cos_phi * cos_theta * sin_ksi - sin_phi * cos_ksi)
         + tp2_lrf_z * sin_theta * sin_ksi);
    double p2_lrf_z = tp2_lrf_x * (cos_phi * sin_theta) + tp2_lrf_z * cos_theta;
    const double E2_lrf = sqrt(
        m2 * m2 + p2_lrf_x * p2_lrf_x + p2_lrf_y * p2_lrf_y
        + p2_lrf_z * p2_lrf_z);
    double p3_lrf_x =
        (tp3_lrf_x * (cos_phi * cos_theta * cos_ksi - sin_phi * sin_ksi)
         - tp3_lrf_z * sin_theta * cos_ksi);
    double p3_lrf_y =
        (tp3_lrf_x * (-cos_phi * cos_theta * sin_ksi - sin_phi * cos_ksi)
         + tp3_lrf_z * (sin_theta * sin_ksi));
    double p3_lrf_z = tp3_lrf_x * cos_phi * sin_theta + tp3_lrf_z * cos_theta;
    const double E3_lrf = sqrt(
        m3 * m3 + p3_lrf_x * p3_lrf_x + p3_lrf_y * p3_lrf_y
        + p3_lrf_z * p3_lrf_z);

    daughter1.E = E1_lrf;
    daughter1.px = p1_lrf_x;
    daughter1.py = p1_lrf_y;
    daughter1.pz = p1_lrf_z;
    daughter2.E = E2_lrf;
    daughter2.px = p2_lrf_x;
    daughter2.py = p2_lrf_y;
    daughter2.pz = p2_lrf_z;
    daughter3.E = E3_lrf;
    daughter3.px = p3_lrf_x;
    daughter3.py = p3_lrf_y;
    daughter3.pz = p3_lrf_z;
    for (auto *daughter : {&daughter1, &daughter2, &daughter3}) {
        boost_from_rest_frame(frame, *daughter);
        daughter->t = decay_time;
        daughter->x = decay_x;
        daughter->y = decay_y;
        daughter->z = decay_z;
    }
}

//! This function sample mother particle mass according to breit-wigner
//! distribution
double particle_decay::sample_breit_wigner(
    const double mass, const double width, const double M_min) const {
    return sample_breit_wigner(mass, width, M_min, *ran_gen_ptr_);
}

double particle_decay::sample_breit_wigner(
    double mass, double width, double M_min, RandomUtil::Random &ran) const {
    // this function sample the Breit Wigner distribution for the mass
    // particle mass sampled is >= M_min
    // From Wiki: https://en.wikipedia.org/wiki/Cauchy_distribution
    // The CDF of BreitWigner (or Cauchy) distribution is
    // y = 1/pi arctan((x-x0)/gamma) + 1/2
    // The inverse CDF is
    // x = x0 + gamma*tan(pi*(y - 1/2))

    // compute the minimum probability given by M_min
    double p_min = atan2((M_min - mass), width / 2.) / M_PI + 0.5;
    // generate a random number from p_min to 1
    double y = (1. - p_min) * (ran.rand_uniform()) + p_min;
    // compute the corresponding sampled mass
    double mass_sampled = mass + width / 2. * tan(M_PI * (y - 0.5));
    return (mass_sampled);
}
