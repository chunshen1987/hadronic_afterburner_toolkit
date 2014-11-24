cmake_minimum_required (VERSION 2.6)
project (HBTcorrelation_MCafterburner)

set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_HOME_DIRECTORY}")

set (CMAKE_CXX_FLAGS "-O3")

add_subdirectory (src)
