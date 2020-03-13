
#include "Analysis.h"
#include "gtest/gtest.h"

#include <string>

namespace {

TEST(Analysis, Constructor) {
    std::string path = "results";
    Analysis test_analysis(path);
    EXPECT_EQ(0, 0);
}


TEST(Analysis, ReadInParameters) {
    std::string path = "results";
    Analysis test_analysis(path);
    test_analysis.UpdateParameterDict("parameters.dat");
    EXPECT_EQ(0, 0);
}


TEST(Analysis, InitializeAnalysis) {
    std::string path = "test_gzip_reader";
    Analysis test_analysis(path);
    test_analysis.UpdateParameterDict("parameters.dat");
    test_analysis.InitializeAnalysis();
    EXPECT_EQ(0, 0);
}


TEST(Analysis, FlowAnalysis) {
    std::string path = "test_gzip_reader";
    Analysis test_analysis(path);
    test_analysis.UpdateParameterDict("parameters.dat");
    test_analysis.InitializeAnalysis();
    test_analysis.FlowAnalysis();
    EXPECT_EQ(0, 0);
}

}  // namespace
