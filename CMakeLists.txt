cmake_minimum_required (VERSION 2.6)
project (hadronic_afterburner_toolkit)

set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_HOME_DIRECTORY}")

find_package ( ZLIB REQUIRED )
if ( ZIB_FOUND )
    include_directories( ${ZLIB_INCLUDE_DIRS} )
endif ( ZIB_FOUND)

message("source directory:" ${CMAKE_SOURCE_DIR})
message("binary directory:" ${CMAKE_BINARY_DIR})

if (APPLE)
    set (CMAKE_CXX_FLAGS "-g -Wall -std=c++11")
else ()
    set (CMAKE_CXX_FLAGS "-g -Wall -pthread -std=c++11")
endif ()

add_subdirectory (src)

option(unit_tests "Build all unit tests." OFF)
if (unit_tests)
    add_subdirectory (unit_tests)
endif (unit_tests)
