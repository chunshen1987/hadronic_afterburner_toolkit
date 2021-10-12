// Copyright @ Chun Shen 2016
// to compile the code:
//     g++ convert_to_binary.cpp -lz -o convert_to_binary.e

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

    string temp_filename = input_filename.substr(0, input_filename.find("."));
    stringstream output_filename;
    output_filename << temp_filename << ".gz";
    stringstream output_filename1;
    output_filename1 << temp_filename << ".bin";

    ifstream urqmd_file(input_filename.c_str());
    if (!urqmd_file.is_open()) {
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
    int urqmd_id, urqmd_iso3, urqmd_charge, n_coll;
    int parent_id, parent_proc_type;
    double mass, t, x, y, z, E, px, py, pz;
    getline(urqmd_file, temp_string);
    while (!urqmd_file.eof()) {
        // skip the header
        for (int i = 0; i < 16; i++)
            getline(urqmd_file, temp_string); 

        // get total number of particles
        getline(urqmd_file, temp_string); 
        stringstream temp1(temp_string);
        int n_particles;
        temp1 >> n_particles;

        // get some information
        int n_info[8];
        getline(urqmd_file, temp_string);
        stringstream temp2(temp_string);
        for (int i = 0; i < 8; i++) {
            temp2 >> n_info[i];
        }
        if (output_mode == "gz") {
            gzprintf(fp_gz, "%d \n", n_particles);
            for (int i = 0; i < 8; i++) {
                gzprintf(fp_gz, "%d ", n_info[i]);
            }
            gzprintf(fp_gz, "\n");
        } else if (output_mode == "binary") {
            outbin.write(reinterpret_cast<char *>(&n_particles),
                         sizeof(int));
            for (int i = 0; i < 8; i++) {
                outbin.write(reinterpret_cast<char *>(&n_info[i]),
                             sizeof(int));
            }
        }

        for (int i = 0; i < n_particles; i++) {
            getline(urqmd_file, temp_string);
            stringstream temp3(temp_string);
            temp3 >> dummy >> dummy >> dummy >> dummy
                  >> dummy >> dummy >> dummy >> dummy
                  >> mass >> urqmd_id >> urqmd_iso3 >> urqmd_charge
                  >> parent_id >> n_coll >> parent_proc_type
                  >> t >> x >> y >> z
                  >> E >> px >> py >> pz;

            if (output_mode == "gz") {
                gzprintf(fp_gz, "%d %d %d %d %d %d ",
                         urqmd_id, urqmd_iso3, urqmd_charge,
                         n_coll, parent_id, parent_proc_type);
                gzprintf(fp_gz,
                         "%.7e %.7e %.7e %.7e %.7e %.7e %.7e %.7e %.7e\n",
                         mass, t, x, y, z, E, px, py, pz);
            } else if (output_mode == "binary") {
                int info_array[] = {urqmd_id, urqmd_iso3, urqmd_charge,
                                    n_coll, parent_id, parent_proc_type};
                float particle_array[] = {static_cast<float>(mass),
                                          static_cast<float>(t),
                                          static_cast<float>(x),
                                          static_cast<float>(y),
                                          static_cast<float>(z),
                                          static_cast<float>(E),
                                          static_cast<float>(px),
                                          static_cast<float>(py),
                                          static_cast<float>(pz)};
                for (int ii = 0; ii < 6; ii++) {
                    outbin.write(reinterpret_cast<char *>(&(info_array[ii])),
                                 sizeof(int));
                }
                for (int ii = 0; ii < 9; ii++) {
                    outbin.write(
                        reinterpret_cast<char *>(&(particle_array[ii])),
                        sizeof(float));
                }
            }
        }
        getline(urqmd_file, temp_string);
    }
    if (output_mode == "gz") {
        gzclose(fp_gz);
    } else if (output_mode == "binary") {
        outbin.close();
    }
    urqmd_file.close();
    return 0;
}
