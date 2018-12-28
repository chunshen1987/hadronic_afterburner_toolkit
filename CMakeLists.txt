cmake_minimum_required (VERSION 3.0)
project (hadronic_afterburner_toolkit LANGUAGES CXX)

set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_HOME_DIRECTORY}")

find_package ( ZLIB REQUIRED )
if ( ZIB_FOUND )
    include_directories( ${ZLIB_INCLUDE_DIRS} )
endif ( ZIB_FOUND)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

message("source directory:" ${CMAKE_SOURCE_DIR})
message("binary directory:" ${CMAKE_BINARY_DIR})

if (APPLE)
    set (CMAKE_CXX_FLAGS "-g -Wall -Wextra")
else ()
    set (CMAKE_CXX_FLAGS "-g -Wall -pthread -Wextra")
endif ()

add_subdirectory (src)

option(unit_tests "Build all unit tests." OFF)
if (unit_tests)
    add_subdirectory (unit_tests)
endif (unit_tests)
