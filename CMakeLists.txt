cmake_minimum_required (VERSION 2.6)
project (hadronic_afterburner_toolkit)

set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_HOME_DIRECTORY}")

find_package ( ZLIB REQUIRED )
if ( ZIB_FOUND )
    include_directories( ${ZLIB_INCLUDE_DIRS} )
endif ( ZIB_FOUND)

set (CMAKE_CXX_FLAGS "-g -Wall")

add_subdirectory (src)
