#include <iostream>
#include <sstream>
#include <fstream>
#include <cstring>
#include "zlib.h"
#include <cstdio>
#include <cstdlib>

using namespace std;

string gzreadline(gzFile gzfp) {
    stringstream line;
    char buffer[1];
    int len = gzread(gzfp, buffer, 1);
    while (len == 1 && buffer[0] != '\n') {
        line << buffer[0];
        len = gzread(gzfp, buffer, 1);
    }
    return(line.str());
}

int main(int argc, char *argv[]) {
    gzFile fp_gz = gzopen("test_gz.dat", "rb");
    string temp_string;

    double dummy;
    int urqmd_id, urqmd_iso3, urqmd_charge, parent_id, n_coll, parent_proc_type;
    double mass, t, x, y, z, E, px, py, pz;
    string test_line = gzreadline(fp_gz);
    cout << test_line << endl;
    test_line = gzreadline(fp_gz);
    cout << test_line << endl;
    test_line = gzreadline(fp_gz);
    cout << test_line << endl;
    test_line = gzreadline(fp_gz);
    cout << test_line << endl;

    //char buffer[sizeof(int)];
    //int len = gzread(fp_gz, buffer, sizeof(int));
    //int N = strtol(buffer, NULL, 10);
    //cout << len << " " << N << "  " << buffer << endl;
    //for (int i = 0; i < 8; i++) {
    //    char buffer1[sizeof(int)];
    //    int len2 = gzread(fp_gz, buffer1, sizeof(int));
    //    int N1 = strtol(buffer1, NULL, 10);
    //    cout << len2 << "  " << N1 << "  " << buffer1 << endl;
    //}
    //while (!urqmd_file.eof()) {
    //    // get total number of particles
    //    getline(urqmd_file, temp_string); 
    //    stringstream temp1(temp_string);
    //    int n_particles[1];
    //    temp1 >> n_particles[0];
    //    fwrite(n_particles, sizeof(int), 1, binary_file);
    //    gzwrite(fp_gz, &n_particles, sizeof(int));

    //    // get some information
    //    int n_info[8];
    //    getline(urqmd_file, temp_string);
    //    stringstream temp2(temp_string);
    //    for (int i = 0; i < 8; i++)
    //        temp2 >> n_info[i];
    //    fwrite(n_info, sizeof(int), 8, binary_file);
    //    gzwrite(fp_gz, &n_info, sizeof(int)*8);

    //    for (int i = 0; i < n_particles[0]; i++) {
    //        getline(urqmd_file, temp_string);
    //        stringstream temp3(temp_string);
    //        temp3 >> dummy >> dummy >> dummy >> dummy
    //              >> dummy >> dummy >> dummy >> dummy
    //              >> mass >> urqmd_id >> urqmd_iso3 >> urqmd_charge
    //              >> parent_id >> n_coll >> parent_proc_type
    //              >> t >> x >> y >> z
    //              >> E >> px >> py >> pz;
    //        int array1[6] = {urqmd_id, urqmd_iso3, urqmd_charge,
    //                        n_coll, parent_id, parent_proc_type};
    //        float array2[9] = {mass, t, x, y, z, E, px, py, pz};
    //        fwrite(array1, sizeof(int), 6, binary_file);
    //        gzwrite(fp_gz, &array1, sizeof(int)*6);
    //        fwrite(array2, sizeof(float), 9, binary_file);
    //        gzwrite(fp_gz, &array2, sizeof(float)*9);
    //    }
    //    getline(urqmd_file, temp_string);
    //}
    gzclose(fp_gz);
    return 0;
}
