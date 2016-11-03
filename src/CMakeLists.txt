add_executable (hadronic_afterburner_tools.e
  main.cpp
  arsenal.cpp
  ParameterReader.cpp
  particleSamples.cpp
  HBT_correlation.cpp
  single_particleSpectra.cpp
  particle_yield_distribution.cpp
  )
target_link_libraries (hadronic_afterburner_tools.e ${LIBS})

install(TARGETS hadronic_afterburner_tools.e DESTINATION ${CMAKE_HOME_DIRECTORY})
