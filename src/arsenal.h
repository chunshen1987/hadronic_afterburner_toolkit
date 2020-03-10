// Version 1.6.2
// Zhi Qiu

#ifndef arsenal_h
#define arsenal_h

#include <vector>
#include <string>

namespace AfterburnerUtil {
    const double hbarc = 0.19733;

    double stringToDouble(std::string);
    std::vector<double> stringToDoubles(std::string);
    std::string toLower(std::string str);
    std::string trim(std::string str);
}

#endif

