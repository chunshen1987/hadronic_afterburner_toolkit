// Copyright @ Chun Shen 2019

#ifndef SRC_HISTOGRAM_H_
#define SRC_HISTOGRAM_H_

#include <vector>
#include <string>

namespace HistUtil{

class Histogram {
 private:
    double bin_min_;
    double bin_max_;
    double bin_width_;
    int nbin_;
    std::vector<int> hist;
    std::vector<double> bin_center;

 public:
    Histogram(double bin_min=0.0, double bin_max=1.0, int nbin=10);

    void fill(double val);
    void output_histogram(std::string filename) const;

};

}


#endif // SRC_HISTOGRAM_H_
