find_package(GTest REQUIRED)

include_directories(${GTEST_INCLUDE_DIRS})

add_executable (HBT_unittest.e HBT_unittest.cc)
target_link_libraries (HBT_unittest.e ${libname} ${ZLIB_LIBRARIES}
                       ${GTEST_BOTH_LIBRARIES})
install(TARGETS HBT_unittest.e DESTINATION ${CMAKE_HOME_DIRECTORY}/unit_tests)

add_executable (particleSamples_unittest.e particleSamples_unittest.cc)
target_link_libraries (particleSamples_unittest.e ${libname} ${ZLIB_LIBRARIES}
                       ${GTEST_BOTH_LIBRARIES})
install(TARGETS particleSamples_unittest.e
        DESTINATION ${CMAKE_HOME_DIRECTORY}/unit_tests)

add_executable (Analysis_unittest.e Analysis_unittest.cc)
target_link_libraries (Analysis_unittest.e ${libname} ${ZLIB_LIBRARIES}
                       ${GTEST_BOTH_LIBRARIES})
install(TARGETS Analysis_unittest.e
        DESTINATION ${CMAKE_HOME_DIRECTORY}/unit_tests)
