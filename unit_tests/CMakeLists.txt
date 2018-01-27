add_executable (HBT_test.e
  HBT_unittest.cc
  ../src/arsenal.cpp
  ../src/ParameterReader.cpp
  ../src/particleSamples.cpp
  ../src/HBT_correlation.cpp
  ../src/particle_decay.cpp
  )
target_link_libraries (HBT_test.e ${LIBS} ${ZLIB_LIBRARIES} gtest gtest_main)

install(TARGETS HBT_test.e DESTINATION ${CMAKE_HOME_DIRECTORY}/unit_tests/)
