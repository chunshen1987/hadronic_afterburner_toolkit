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
target_link_libraries (${libname} ${ZLIB_LIBRARIES})

if (link_with_lib)
    add_executable (${exename} main.cpp)
    target_link_libraries (${exename} ${libname} ${ZLIB_LIBRARIES})
else (link_with_lib)
    add_executable (${exename} main.cpp ${SOURCES})
    target_link_libraries (${exename} ${ZLIB_LIBRARIES})
endif (link_with_lib)

if (APPLE)
    set_target_properties (${libname} PROPERTIES COMPILE_FLAGS "-DAPPLE")
    set_target_properties (${exename} PROPERTIES COMPILE_FLAGS "-DAPPLE")
endif (APPLE)

install(TARGETS ${exename} DESTINATION ${CMAKE_HOME_DIRECTORY})
install(TARGETS ${libname} DESTINATION ${CMAKE_HOME_DIRECTORY})
