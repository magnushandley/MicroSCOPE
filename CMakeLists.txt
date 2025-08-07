cmake_minimum_required(VERSION 3.16)

# ---------------------------------------------------------------------------
#  Project-wide settings
# ---------------------------------------------------------------------------
project(uboone_bsm VERSION 0.1
                              DESCRIPTION "Generic BSM analysis framework for MicroBooNE"
                              LANGUAGES CXX)

set(CMAKE_CXX_STANDARD         17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Optionally build everything static (useful for grid jobs)
option(STATIC "Link static libraries where possible" OFF)
if(STATIC)
  set(CMAKE_FIND_LIBRARY_SUFFIXES .a ${CMAKE_FIND_LIBRARY_SUFFIXES})
endif()

# ---------------------------------------------------------------------------
#  Locate ROOT
# ---------------------------------------------------------------------------
find_package(ROOT 6.20)

# `ROOT_USE_FILE` pulls in the compile flags for ACLiC-compatible dictionary code.
include(${ROOT_USE_FILE})

# ---------------------------------------------------------------------------
#  Add subdirectories
# ---------------------------------------------------------------------------
add_subdirectory(src)   # library and dictionaries
add_subdirectory(run)   # tiny main() wrappers / macros

# ---------------------------------------------------------------------------
#  Installation (optional)
# ---------------------------------------------------------------------------
install(DIRECTORY include/ DESTINATION include)

# Export a CMake package for downstream users (optional)
export(PACKAGE bsm_uboone_framework)
