// Copyright @ Chun Shen 2016
// to compile the code:
//     g++ convert_to_binary_SMASH.cpp -lz -o convert_to_binary_SMASH.e

#include <iostream>
#include <sstream>
#include <fstream>
#include <cstdlib>
#include "zlib.h"

using namespace std;

int main(int argc, char *argv[]) {
    string input_filename;
    if (argc < 2) {
        cout << "Usage: " << argv[0] << " filename  [format]" << endl;
        cout << "--- format: gz [default], binary" << endl;
        exit(0);
    } else {
        input_filename = argv[1];
    }
    string output_mode = "gz";
    if (argc == 3) {
        output_mode = argv[2];
    }

    // get input filename for output filename
    string temp_filename = input_filename.substr(0, input_filename.find("."));
    stringstream output_filename;
    output_filename << temp_filename << ".gz";
    stringstream output_filename1;
    output_filename1 << temp_filename << ".bin";

    ifstream SMASH_file(input_filename.c_str());
    if (!SMASH_file.is_open()) {
        cout << "can not open " << input_filename << endl;
        cout << "exit" << endl;
        exit(1);
    }

    ofstream outbin;
    if (output_mode == "binary") {
        outbin.open(output_filename1.str().c_str(), ios::out | ios::binary);
    }

    gzFile fp_gz;
    if (output_mode == "gz") {
        fp_gz = gzopen(output_filename.str().c_str(), "wb");
    }

    string temp_string;
    double dummy;
    string str_dummy;
    int pdg, charge, SMASH_proc_Id, pdg_mother1, pdg_mother2;
    double mass, t, x, y, z, E, px, py, pz;
    double t_form, t_fz, x_fz, y_fz, z_fz;

    // start to read in SMASH outputs
    // skip the header
    for (int i = 0; i < 3; i++) {
        getline(SMASH_file, temp_string);
    }

    int nev = 1;
    // get total number of particles
    getline(SMASH_file, temp_string);
    while (!SMASH_file.eof()) {
        stringstream temp1(temp_string);
        int n_particles[1];
        temp1 >> str_dummy >> str_dummy >> str_dummy >> str_dummy
              >> n_particles[0];
        cout << "nev = " << nev << ", npart = " << n_particles[0] << endl;
        if (output_mode == "gz") {
            gzprintf(fp_gz, "%d \n", n_particles[0]);
        } else if (output_mode == "binary") {
            outbin.write((char*) &(n_particles[0]), sizeof(int));
        }

        for (int i = 0; i < n_particles[0]; i++) {
            getline(SMASH_file, temp_string);
            stringstream temp3(temp_string);
            temp3 >> t >> x >> y >> z
                  >> mass >> E >> px >> py >> pz
                  >> pdg >> dummy >> charge >> dummy
                  >> t_form >> dummy >> dummy >> SMASH_proc_Id >> t_fz
                  >> pdg_mother1 >> pdg_mother2;

            // compute the particle's last collision position
            t_fz = std::max(t_fz, t_form);
            double vx = px/E;
            double vy = py/E;
            double vz = pz/E;
            double dt_backpropagate = t - t_fz;
            x_fz = x - vx*dt_backpropagate;
            y_fz = y - vy*dt_backpropagate;
            z_fz = z - vz*dt_backpropagate;

            if (output_mode == "gz") {
                gzprintf(fp_gz, "%d %d %d %d %d ", pdg, charge, SMASH_proc_Id,
                         pdg_mother1, pdg_mother2);
                gzprintf(fp_gz,
                         "%.7e %.7e %.7e %.7e %.7e %.7e %.7e %.7e %.7e\n",
                         mass, t_fz, x_fz, y_fz, z_fz, E, px, py, pz);
            } else if (output_mode == "binary") {
                int info_array [] = {pdg, charge, SMASH_proc_Id,
                                     pdg_mother1, pdg_mother2};
                float particle_array [] = {mass, t_fz, x_fz, y_fz, z_fz,
                                           E, px, py, pz};
                for (int ii = 0; ii < 5; ii++) {
                    outbin.write((char*) &(info_array[ii]), sizeof(int));
                }
                for (int ii = 0; ii < 9; ii++) {
                    outbin.write((char*) &(particle_array[ii]), sizeof(float));
                }
            }
        }
        getline(SMASH_file, temp_string);  // get the last info line

        // get total number of particles
        getline(SMASH_file, temp_string);
        nev++;
    }
    if (output_mode == "gz") {
        gzclose(fp_gz);
    } else if (output_mode == "binary") {
        outbin.close();
    }
    SMASH_file.close();
    return 0;
}
