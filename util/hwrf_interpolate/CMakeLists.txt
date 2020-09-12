#=======================================================================
#$$$ CMAKEFILE DOCUMENTATION BLOCK
# Biju Thomas
# Email: biju.thomas@noaa.gov
#=======================================================================

set(fortran_srcs 
  namelist.f90
  general_tll2xy_mod.f90
  wrf_nmm_io.f90
  wrf_nmm_interpolate.f90
  wrf_interpolate_main.f90)

set(exe_name hwrf_interpolate.exe)
set(exec_dir ${CMAKE_SOURCE_DIR}/exec)
set(mpeu_util_lib ${CMAKE_SOURCE_DIR}/build/src/gsi/CMakeFiles/gsilib_shrd.dir/mpeu_util.F90.o)


if(CMAKE_Fortran_COMPILER_ID MATCHES "^(Intel)$")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -free -O3 -convert big_endian -DLINUX -fp-model precise -assume byterecl ${OpenMP_Fortran_FLAGS}")
endif()

add_executable(${exe_name} ${fortran_srcs})

target_include_directories(
  ${exe_name} PUBLIC
  ${CMAKE_SOURCE_DIR}/lib 
  ${CMAKE_SOURCE_DIR}/include
  ${MPI_Fortran_INCLUDE_DIRS} 
  ${MPI_Fortran_INCLUDE_PATH} 
  ${NETCDF_INCLUDES})

target_link_libraries(
  ${exe_name}
  ${mpeu_util_lib}
  ${GSILIB} 
  ${GSISHAREDLIB} 
  ${NETCDF_LIBRARIES_F90} 
  ${NETCDF_LIBRARIES}  
  ${HDF5_HL_LIBRARIES} 
  ${HDF5_LIBRARIES} 
  ${ZLIB_LIBRARIES})


install(TARGETS ${exe_name} DESTINATION ${exec_dir})
