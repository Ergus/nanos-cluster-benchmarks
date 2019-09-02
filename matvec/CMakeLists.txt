
find_package (LIBEREP REQUIRED)

include_directories (${LIBEREP_INCLUDE})

add_executable (matvec_strong matvec_strong.cpp)
target_link_libraries (matvec_strong ${LIBEREP_LIBRARY})
