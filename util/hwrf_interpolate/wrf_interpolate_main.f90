program wrf_interpolate_main

  use kinds
  use mpimod
  use namelist, only : namelistparams,intrp_from_file,intrp_to_file
  use constants, only: init_constants_derived, init_constants
  use wrf_nmm_io
  use general_sub2grid_mod, only: sub2grid_info 
  use wrf_nmm_interpolate, only : interpolate_s2d

  implicit none

  !-------------------------------------------------------------------
  !$$$  main program documentation block
  !
  ! main program: wrf_interpolate
  ! PRGMMR: TONG                           DATE: 2016-02-17
  !
  ! abstract: The wrf interoplation code interpolate wrf fields from
  !   one wrf nmm grid to another wrf nmm grid using modules from
  !   GSI
  !
  !$$$  end documentation bloc
  !-------------------------------------------------------------------

  ! Define variables computed within routine

  real(r_kind)                                       :: exectime_start
  real(r_kind)                                       :: exectime_finish
  type(sub2grid_info) :: grd_src, grd_dst
  type(wrfnmm_variables) :: var_src, var_dst

  !=====================================================================

  ! Initialize MPI session

  call mpi_init(ierror)

  ! Get the number of processes
  call mpi_comm_size(mpi_comm_world,npe,ierror)

  ! Get the number of processes
  call mpi_comm_rank(mpi_comm_world,mype,ierror)

  if(mype==0) then
     call cpu_time(exectime_start)
  end if

  ! Read in namelist
  call namelistparams()

  ! Initialize constants
  call init_constants_derived
  call init_constants(.true.)

  ! Read input grids
  if(mype == 0)then
     call wrf_nmm_interface(intrp_from_file)
     call wrf_nmm_interface(intrp_to_file)
  end if

  call mpi_barrier(mpi_comm_world,ierror)

  call wrf_nmm_read(intrp_from_file,grd_src,var_src)
  call wrf_nmm_read(intrp_to_file,grd_dst,var_dst)

  ! Run interpolation
  call interpolate_s2d(grd_src,var_src,grd_dst,var_dst)

  ! Update target grid
  call wrf_nmm_write(intrp_to_file,grd_dst,var_dst)

  if(mype==0) then
     call cpu_time(exectime_finish)
     write(6,500) exectime_finish - exectime_start
  end if 

  ! Enable the root task to catch up from I/O and calculations

  call mpi_barrier(mpi_comm_world,ierror)

  ! Finalize MPI session

  call mpi_finalize(ierror)

  !=====================================================================

  ! Define format statements

500 format('MAIN: Execution time: ', f13.5, ' seconds.')

  !=====================================================================

end program wrf_interpolate_main
