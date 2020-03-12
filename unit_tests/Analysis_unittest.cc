
#include "Analysis.h"
#include "gtest/gtest.h"

#include <string>

namespace {

TEST(Analysis, DefaultConstructor) {
    std::string path = "results";
    Analysis test_analysis(path);
    EXPECT_EQ(0, 0);
}

}  // namespace
