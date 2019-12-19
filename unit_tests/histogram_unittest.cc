#include "../src/histogram.h"
#include "gtest/gtest.h"

using HistUtil::Histogram;

namespace {

TEST(histogram, DefaultConstructor) {
    Histogram hist(0.0, 1.0, 10);
    EXPECT_EQ(0, 0);
}

}
