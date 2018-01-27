#include "../src/HBT_correlation.h"
#include "gtest/gtest.h"

namespace {
// Tests the Glauber class

TEST(HBT_correlation, DefaultConstructor) {
    ParameterReader *paraRdr = new ParameterReader;
    paraRdr->readFromFile("parameters.dat");
    string path="test_files";
    particleSamples particle_list(paraRdr, path);
        
    HBT_correlation test(paraRdr, path, &particle_list);
    EXPECT_EQ(0, 0);
}

//TEST(Glauber, determine_collision_point_deceleration) {
//    Parameters *param = new Parameters();
//    param->set_QCD_string_production_mode(1);
//    param->set_evolve_QCD_string_mode(1);
//
//    Random *random = new Random();
//
//    Glauber *glauber = new Glauber(param, random);
//
//    int collided = 0;
//    double t_test, z_test;
//    // free-streaming
//    collided = (
//        glauber->get_collision_point_deceleration_intercept_deceleration(
//                    1.0, 0.0, 1.5, 0.0, 3.0, -1.5, 0.0, &t_test, &z_test));
//    EXPECT_EQ(1, collided);
//    EXPECT_NEAR(1.5, z_test, 1e-8);
//    EXPECT_DOUBLE_EQ(2.6571870894737679, t_test);
//    
//    // left deceleration
//    collided = glauber->get_collision_point_deceleration_intercept_deceleration(
//                        1.0, 0.0, 1.8, 0.5, 3.0, -1.5, 0.0, &t_test, &z_test);
//    EXPECT_EQ(1, collided);
//    EXPECT_DOUBLE_EQ(1.4888921827547192, z_test);
//    EXPECT_DOUBLE_EQ(2.6694589103611595, t_test);
//
//    // right deceleration
//    collided = glauber->get_collision_point_deceleration_intercept_deceleration(
//                        1.0, 0.0, 1.8, 0.0, 3.0, -1.5, -0.5, &t_test, &z_test);
//    EXPECT_EQ(1, collided);
//    EXPECT_DOUBLE_EQ(1.7013652023362309, z_test);
//    EXPECT_DOUBLE_EQ(2.796952257645283, t_test);
//
//    // both deceleration
//    collided = glauber->get_collision_point_deceleration_intercept_deceleration(
//                        1.0, 0.0, 1.6, 0.5, 3.0, -1.6, -0.5, &t_test, &z_test);
//    EXPECT_EQ(1, collided);
//    EXPECT_NEAR(1.5, z_test, 1e-8);
//    EXPECT_DOUBLE_EQ(2.9744075156729188, t_test);
//}
//
//
//TEST(Glauber, get_zlab_from_tlab_deceleration) {
//    Parameters *param = new Parameters();
//    param->set_QCD_string_production_mode(1);
//    param->set_evolve_QCD_string_mode(1);
//
//    Random *random = new Random();
//
//    Glauber *glauber = new Glauber(param, random);
//
//    double t_test, z_test;
//    double y_f;
//    z_test = glauber->get_zlab_from_tlab_deceleration(
//                            1.0, 0.0, 1.6, 0.5, 2.9744075156729188, &y_f);
//    cout << z_test << "  " << y_f << endl;
//    EXPECT_NEAR(1.5, z_test, 1e-8);
//    EXPECT_NEAR(0.3911125432748781, y_f, 1e-8);
//    z_test = glauber->get_zlab_from_tlab_deceleration(
//                            1.0, 3.0, -1.6, -0.5, 2.9744075156729188, &y_f);
//    cout << z_test << "  " << y_f << endl;
//    EXPECT_NEAR(1.5, z_test, 1e-8);
//    EXPECT_NEAR(-0.3911125432748781, y_f, 1e-8);
//}

}  // namespace
