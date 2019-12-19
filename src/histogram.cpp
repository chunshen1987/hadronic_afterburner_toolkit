// Copyright @ Chun Shen 2019

#include "histogram.h"
#include <fstream>
#include <iomanip>

namespace HistUtil {

Histogram::Histogram(double bin_min, double bin_max, int nbin) {
    hist.clear();
    bin_center.clear();

    bin_min_ = bin_min;
    bin_max_ = bin_max;
    nbin_ = nbin;

    hist.resize(nbin_, 0);
    bin_center.resize(nbin_);

    bin_width_ = (bin_max - bin_min)/nbin;
    for (int i = 0; i < nbin_; i++) {
        bin_center[i] = bin_min_ + (i + 0.5)*bin_width_;
    }
}


void Histogram::fill(double val) {
    if (val < bin_min_ || val >= bin_max_) return;
    int bin_idx = (val - bin_min_)/bin_width_;
    hist[bin_idx]++;
}


void Histogram::output_histogram(std::string filename) const {
    std::ofstream outf(filename.c_str());
    for (int i = 0; i < nbin_; i++) {
        outf << std::scientific << std::setw(10) << std::setprecision(6)
             << bin_center[i] << "  " << hist[i] << std::endl;
    }
    outf.close();
}

}
