add_executable (HBT_afterburner.e
  main.cpp
  arsenal.cpp
  gauss_quadrature.cpp
  ParameterReader.cpp
  particleSamples.cpp
  HBT_correlation.cpp
  single_particleSpectra.cpp
  )
target_link_libraries (HBT_afterburner.e ${LIBS})

install(TARGETS HBT_afterburner.e DESTINATION ${CMAKE_HOME_DIRECTORY})
