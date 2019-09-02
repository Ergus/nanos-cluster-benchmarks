include(FindPackageHandleStandardArgs)


find_program(MCXX mcxx)
find_program(MCC mcc)
find_program(MFC mfc)

if (MCXX AND MCC AND MFC)
  set(CMAKE_CXX_COMPILER ${MCXX} CACHE INTERNAL "" FORCE)
  set(CMAKE_Fortran_COMPILER ${MFC} CACHE INTERNAL "" FORCE)
  set(CMAKE_C_COMPILER ${MCC} CACHE INTERNAL "" FORCE)

  add_definitions(--ompss-2 -k -DNANOS6)
  set(CMAKE_EXE_LINKER_FLAGS --ompss-2)

  message ("-- Found MCXX : ${MCXX}")
  message ("-- Found  MCC : ${MCC}")
  message ("-- Found  MFC : ${MFC}")
endif ()

find_package_handle_standard_args(MCC
  "Mercurium not found in path"
  MCXX MCC MFC
  )

