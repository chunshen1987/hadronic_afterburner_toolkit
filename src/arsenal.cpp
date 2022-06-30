// Ver 1.6.2
// Zhi Qiu
/*===========================================================================
Change logs: see arsenal.h
============================================================================*/

#include <iostream>
#include <string>
#include <sstream>
#include <stdlib.h>
#include <vector>
#include <iomanip>
#include <cstdarg>

#include "arsenal.h"

namespace AfterburnerUtil {

//**********************************************************************
std::vector<double> stringToDoubles(std::string str) {
    // Return a vector of doubles from the string "str". "str" should
    // be a string containing a line of data.
    // add a blank at the end so the last data will be read
    std::stringstream sst(str+" "); 
    std::vector<double> valueList;
    double val;
    sst >> val;
    while (sst.eof() == false) {
        valueList.push_back(val);
        sst >> val;
    }
    return valueList;
}


//**********************************************************************
double stringToDouble(std::string str) {
    // Return the 1st doubles number read from the string "str". 
    // "str" should be a string containing a line of data.
    // add a blank at the end so the last data will be read
    std::stringstream sst(str+" "); 
    double val;
    sst >> val;
    return val;
}


//**********************************************************************
std::string toLower(std::string str) {
    // Convert all character in string to lower case
    std::string tmp = str;
    for (std::string::iterator it=tmp.begin(); it<=tmp.end(); it++)
        *it = tolower(*it);
    return tmp;
}


//**********************************************************************
std::string trim(std::string str) {
    // Convert all character in string to lower case
    std::string tmp = str;
    long number_of_char = 0;
    for (size_t ii = 0; ii < str.size(); ii++) {
        if (str[ii] != ' ' && str[ii] != '\t') {
            tmp[number_of_char]=str[ii];
            number_of_char++;
        }
    }
    tmp.resize(number_of_char);
    return tmp;
}


void resize2DVector(std::vector<std::vector<double>> &vec2D,
                    const int nx, const int ny, const double val) {
    for (unsigned int i = 0; i < vec2D.size(); i++) {
        vec2D[i].clear();
    }
    vec2D.clear();
    for (int ix = 0; ix < nx; ix++) {
        std::vector<double> tmp(ny, val);
        vec2D.push_back(tmp);
    }
}


void resize3DVector(std::vector<std::vector<std::vector<double>>> &vec3D,
                    const int nx, const int ny, const int nz,
                    const double val) {
    for (unsigned int i = 0; i < vec3D.size(); i++) {
        for (unsigned int j = 0; j < vec3D[i].size(); j++) {
            vec3D[i][j].clear();
        }
        vec3D[i].clear();
    }
    vec3D.clear();
    for (int ix = 0; ix < nx; ix++) {
        std::vector<std::vector<double>> tmp2D;
        for (int iy = 0; iy < ny; iy++) {
            std::vector<double> tmp(nz, val);
            tmp2D.push_back(tmp);
        }
        vec3D.push_back(tmp2D);
    }
}

}
