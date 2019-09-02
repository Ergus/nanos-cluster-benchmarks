###############################################################################
# Find libErep
#
# This sets the following variables:
# LIBEREP_FOUND        - True if liberep was found.
# LIBEREP_INCLUDES     - Directories containing the liberep include files.
# LIBEREP_LIBRARIES    - Libraries needed to use liberep.
###############################################################################

INCLUDE(FindPackageHandleStandardArgs)

MESSAGE(STATUS "CMAKE_LIBEREP_ROOT: ${CMAKE_LIBEREP_ROOT}")
MESSAGE(STATUS "LIBEREP_ROOT: $ENV{LIBEREP_ROOT}")

# Set extended paths
SET(INC_PATHS
  /usr/include
  /usr/local/include
  /usr/include )

SET(LIB_PATHS
  /usr/lib
  /usr/lib64
  /usr/local/lib
  /usr/lib/x86_64-linux-gnu)

# Find where are the header files.
FIND_PATH(LIBEREP_INCLUDE
  NAMES CommandLine.hpp
  HINTS ${CMAKE_LIBEREP_ROOT} ENV LIBEREP_ROOT
  PATHS ${INC_PATHS} ENV CPLUS_INCLUDE_PATH
  DOC "LIBEREP include path"
  PATH_SUFFIXES erep
  )

FIND_LIBRARY(LIBEREP_LIBRARY
  NAMES erep
  HINTS ${CMAKE_LIBEREP_ROOT} ENV LIBEREP_ROOT
  PATHS ${LIB_PATHS} ENV LD_LIBRARY_PATH
  DOC "LIBEREP library"
  PATH_SUFFIXES lib lib64 lib32)

FIND_PACKAGE_HANDLE_STANDARD_ARGS(LIBEREP
  LIBEREP_INCLUDE LIBEREP_LIBRARY)

MARK_AS_ADVANCED(LIBEREP_INCLUDE LIBEREP_LIBRARY)

