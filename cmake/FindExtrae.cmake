###############################################################################
# Find OpenBlas
#
# This sets the following variables:
# EXTRAE_FOUND        - True if Extrae was found.
# EXTRAE_INCLUDES     - Directories containing the Extrae include files.
# EXTRAE_LIBRARIES    - Libraries needed to use Extrae.
###############################################################################

include (FindPackageHandleStandardArgs)

set (EXTRAE_HOME $ENV{EXTRAE_HOME})
message (STATUS "EXTRAE_HOME: ${EXTRAE_HOME}")

if (EXTRAE_HOME)

  # Find where are the header files.
  find_path (EXTRAE_INCLUDE
    NAMES extrae_user_events.h
    HINTS ENV EXTRAE_HOME
    PATHS ENV C_INCLUDE_PATH
    DOC "Extrae include path"
    PATH_SUFFIXES include
    NO_DEFAULT_PATH
    )

  find_library (EXTRAE_LIBRARY
    NAMES ${Extrae_FIND_COMPONENTS}
    HINTS ENV EXTRAE_HOME
    PATHS ENV LD_LIBRARY_PATH
    REQUIRED
    DOC "Extrae library"
    PATH_SUFFIXES lib lib64 lib32
    NO_DEFAULT_PATH)

  find_package_handle_standard_args (Extrae
    REQUIRED_VARS EXTRAE_LIBRARY EXTRAE_INCLUDE)

endif ()

mark_as_advanced(EXTRAE_INCLUDES EXTRAE_LIBRARIES)

set (EXTRAE_INCLUDES  ${EXTRAE_INCLUDE})
set (EXTRAE_LIBRARIES ${EXTRAE_LIBRARY})
