add_executable (hadronic_afterburner_tools.e
  main.cpp
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
  )
if (APPLE)
    set_target_properties (hadronic_afterburner_tools.e PROPERTIES COMPILE_FLAGS "-DAPPLE")
endif (APPLE)
target_link_libraries (hadronic_afterburner_tools.e ${LIBS} ${ZLIB_LIBRARIES})

install(TARGETS hadronic_afterburner_tools.e DESTINATION ${CMAKE_HOME_DIRECTORY})
