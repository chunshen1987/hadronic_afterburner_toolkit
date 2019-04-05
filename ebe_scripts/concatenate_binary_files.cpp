// Copyright @ Chun Shen 2016
// to compile the code:
//     g++ concatenate_binary_files.cpp -lz -o concatenate_binary_files.e

#include <iostream>
#include <sstream>
#include <fstream>
#include <cstdlib>
#include "zlib.h"

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
    string input_filename1;
    string input_filename2;
    if (argc < 3) {
        cout << "Usage: " << argv[0] << " filename1 filename2" << endl;
        cout << "This program concatenates filename2 behind the filename1."
             << endl;
        exit(0);
    } else {
        input_filename1 = argv[1];
        input_filename2 = argv[2];
    }
    cout << "Concatenating " << input_filename2 << " into "
         << input_filename1 << " ... " << endl;
    gzFile f1_gz = gzopen(input_filename1.c_str(), "ab");
    gzFile f2_gz = gzopen(input_filename2.c_str(), "rb");
    if (!f1_gz) {
        cout << "can not open " << input_filename1 << endl;
        cout << "exit." << endl;
        exit(1);
    }
    if (!f2_gz) {
        cout << "can not open " << input_filename2 << endl;
        cout << "exit." << endl;
        exit(1);
    }
    string temp_string = gzreadline(f2_gz);
    int n_particle;
    stringstream temp1(temp_string);
    temp1 >> n_particle;
    while (!gzeof(f2_gz)) {
        gzprintf(f1_gz, "%s\n", temp_string.c_str());
        temp_string = gzreadline(f2_gz);
    }
    gzclose(f1_gz);
    gzclose(f2_gz);
    return 0;
}
