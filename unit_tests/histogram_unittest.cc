#include "../src/histogram.h"
#include "gtest/gtest.h"

using HistUtil::Histogram;

namespace {


TEST(histogram, DefaultConstructor) {
    Histogram hist1;
    Histogram hist2(0.0, 1.0, 1);
    Histogram hist3(0.0, 1.0, 100);
    EXPECT_EQ(0, 0);
}


TEST(histogram, fill) {
    Histogram hist1;
    hist1.fill(0.05);
    auto hist_arr = hist1.get_bin_count_arr();
    EXPECT_EQ(hist_arr[0], 1);
    EXPECT_EQ(hist_arr[2], 0);
    EXPECT_EQ(hist1.get_total_count(), 1);
    hist1.fill(-0.05);
    auto hist_0 = hist1.get_bin_count(0);
    EXPECT_EQ(hist_0, 1);
    EXPECT_EQ(hist1.get_total_count(), 1);
    hist1.fill(0.8);
    EXPECT_EQ(hist1.get_total_count(), 2);
}


}
