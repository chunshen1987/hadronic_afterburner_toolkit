set (SOURCES
  arsenal.cpp
  ParameterReader.cpp
  particleSamples.cpp
  HBT_correlation.cpp
  single_particleSpectra.cpp
  particle_yield_distribution.cpp
  particle_decay.cpp
  BalanceFunction.cpp
  Random.cpp
  pretty_ostream.cpp
  Analysis.cpp
  )

add_library(${libname} SHARED ${SOURCES})
add_executable (${exename} main.cpp ${SOURCES})
if (APPLE)
    set_target_properties (${libname} PROPERTIES COMPILE_FLAGS "-DAPPLE")
    set_target_properties (${exename} PROPERTIES COMPILE_FLAGS "-DAPPLE")
endif (APPLE)
target_link_libraries (${libname} ${ZLIB_LIBRARIES})
target_link_libraries (${exename} ${ZLIB_LIBRARIES})

install(TARGETS ${exename} DESTINATION ${CMAKE_HOME_DIRECTORY})
install(TARGETS ${libname} DESTINATION ${CMAKE_HOME_DIRECTORY})
