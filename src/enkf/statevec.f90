module statevec
!$$$  module documentation block
!
! module: statevec            read ensemble members, distribute each
!                             to each task.  Collect updated ensemble
!                             members on root task, write out.
!
! prgmmr: whitaker         org: esrl/psd               date: 2009-02-23
!
! abstract: ensemble IO.
!
! Public Subroutines:
!  read_ensemble: read ensemble members on root, distribute pieces (defined by module loadbal)
!    to each task.
!  write_ensemble: retrieve pieces of updated ensemble from each task on root,
!    write out. Optionally save ensemble mean analysis increment.
!  statevec_cleanup: deallocate allocatable arrays.
!
! Public Variables:
!  nanals: (integer scalar) number of ensemble members (from module params)
!  npts_max: (integer scalar) maximum number of grid points assigned to a task.
!  nlevs: number of analysis vertical levels (from module params).
!  nvars: number of 3d 'non-tracer' variables updated by analysis (from module params).
!  ndim: (nvars + ntrac_update) * nlevs (from module params).
!  anal_chunk(nanals,npts_max,ndim): real array of ensemble perturbations 
!   updated on each task.
!  anal_chunk_prior(nanals,npts_max,ndim): real array of prior ensemble 
!   perturbations.  Before analysis anal_chunk=anal_chunk_prior, after
!   analysis anal_chunk contains posterior perturbations.
!  ensmean_chunk(npts_max,ndim): real array containing pieces of ensemble 
!   mean to be updated on each task.
!  ensmean_chunk_prior(npts_max,ndim): as above, for ensemble mean prior.
!   Before analysis ensmean_chunk=ensmean_chunk_prior, after analysis
!   ensmean_chunk contains posterior ensemble mean.
!   
! Modules Used: mpisetup, params, kinds, loadbal, gridio, gridinfo
!
! program history log:
!   2009-02-23  Initial version.
!   2009-11-28  revamped to improve IO speed
!
! attributes:
!   language: f95
!
!$$$

use mpisetup
use params, only: nlevs,nvars,ndim,&
                  nanals,pseudo_rh,massbal_adjust
use kinds, only: r_kind, i_kind, r_double, r_single
use loadbal, only: npts_max,indxproc,numptsperproc
use gridio, only: readgriddata, writegriddata
use gridinfo, only: lonsgrd, latsgrd, ptop, npts, nvarhumid
use enkf_obsmod, only: nobsgood
implicit none
private
public :: read_ensemble, write_ensemble, statevec_cleanup
real(r_single),public, allocatable, dimension(:,:,:) :: anal_chunk, anal_chunk_prior
real(r_single),public, allocatable, dimension(:,:) :: ensmean_chunk, ensmean_chunk_prior
real(r_single),public, allocatable, dimension(:,:) :: grdin
real(r_double),public, allocatable, dimension(:,:) :: qsat
integer(i_kind), allocatable, dimension(:) :: scounts, displs, rcounts

contains

subroutine read_ensemble()
! read ensemble members on IO tasks,
! distribute pieces (defined by module loadbal) to each task.
! for now, first nanals tasks are IO tasks.
implicit none
real(r_single), allocatable, dimension(:) :: sendbuf,recvbuf
real(r_double) t1,t2
integer(i_kind) nanal,nn,i,n
! npts,nlevs,ntrac arrays
integer(i_kind) ierr, np

! must at least nanals tasks allocated.
if (numproc < nanals) then
  print *,'need at least nanals =',nanals,'MPI tasks, exiting ...'
  call mpi_barrier(mpi_comm_world,ierr)
  call mpi_finalize(ierr)
end if
if (npts < numproc) then
  print *,'cannot allocate more than npts =',npts,'MPI tasks, exiting ...'
  call mpi_barrier(mpi_comm_world,ierr)
  call mpi_finalize(ierr)
end if

allocate(scounts(0:numproc-1))
allocate(displs(0:numproc-1))
allocate(rcounts(0:numproc-1))
! only IO tasks send any data.
! scounts is number of data elements to send to processor np.
! rcounts is number of data elements to recv from processor np.
! displs is displacement into send array for data to go to proc np
do np=0,numproc-1
   displs(np) = np*npts_max*ndim
enddo
if (nproc <= nanals-1) then
   scounts = npts_max*ndim
else
   scounts = 0
endif
! displs is also the displacement into recv array for data to go into anal_chunk on
! task np.
do np=0,numproc-1
   if (np <= nanals-1) then
      rcounts(np) = npts_max*ndim
   else
      rcounts(np) = 0
   end if
enddo

! allocate array to hold pieces of state vector on each proc.
allocate(anal_chunk(nanals,npts_max,ndim))
if (nproc == 0) print *,'anal_chunk size = ',size(anal_chunk)
! send and receive buffers.
allocate(sendbuf(numproc*npts_max*ndim))
allocate(recvbuf(nanals*npts_max*ndim))

! read in whole state vector on i/o procs - keep in memory 
! (needed in write_ensemble)
if (nproc <= nanals-1) then
   allocate(grdin(npts,ndim))
   allocate(qsat(npts,nlevs))
   nanal = nproc + 1
   t1 = mpi_wtime()
   call readgriddata(nanal,grdin,qsat)
   if (nproc == 0) then
     t2 = mpi_wtime()
     print *,'time in readgridata on root',t2-t1,'secs'
     t1 = mpi_wtime()
   end if
   !print *,'min/max ps ens mem',nanal,'=',&
   !         minval(grdin(:,ndim)),maxval(grdin(:,ndim))
   if (pseudo_rh .and. nvarhumid > 0) then
      ! create normalized humidity analysis variable.
      grdin(:,(nvarhumid-1)*nlevs+1:nvarhumid*nlevs) = grdin(:,(nvarhumid-1)*nlevs+1:nvarhumid*nlevs)/qsat
   end if
   ! fill up send buffer.
   do np=1,numproc
    do nn=1,ndim
     do i=1,numptsperproc(np)
      n = ((np-1)*ndim + (nn-1))*npts_max + i
      sendbuf(n) = grdin(indxproc(np,i),nn)
     enddo
    enddo
   enddo
end if
call mpi_alltoallv(sendbuf, scounts, displs, mpi_real4, recvbuf, rcounts, displs,&
                   mpi_real4, mpi_comm_world, ierr)
deallocate(sendbuf)
allocate(anal_chunk_prior(nanals,npts_max,ndim))
allocate(ensmean_chunk(npts_max,ndim))
allocate(ensmean_chunk_prior(npts_max,ndim))
ensmean_chunk = 0.

!==> compute ensemble of first guesses, remove mean from anal.
!$omp parallel do schedule(dynamic,1)  private(nn,i,nanal,n)
do nn=1,ndim 
   do i=1,numptsperproc(nproc+1)

      do nanal=1,nanals
         n = ((nanal-1)*ndim + (nn-1))*npts_max + i
         anal_chunk(nanal,i,nn) = recvbuf(n)
      enddo
      ensmean_chunk(i,nn) = sum(anal_chunk(:,i,nn))/float(nanals)
      ensmean_chunk_prior(i,nn) = ensmean_chunk(i,nn)

! remove mean from ensemble.
      do nanal=1,nanals
         anal_chunk(nanal,i,nn) = anal_chunk(nanal,i,nn)-ensmean_chunk(i,nn)
         anal_chunk_prior(nanal,i,nn)=anal_chunk(nanal,i,nn)
      end do

   end do
end do
deallocate(recvbuf)

if (nproc == 0) then
  t2 = mpi_wtime()
  print *,'time to scatter state on root',t2-t1,'secs'
endif


end subroutine read_ensemble

subroutine write_ensemble()
! retrieve pieces of updated ensemble from each task to IO tasks,
! write out each ensemble member to a separate file.
! for now, first nanals tasks are IO tasks.
implicit none
real(r_single), allocatable, dimension(:) :: sendbuf, recvbuf
real(r_kind), allocatable, dimension(:,:) :: ensmean
real(r_double) t1,t2
integer(i_kind) nanal,i,nvar
integer(i_kind) ierr, np, n, nn

! all tasks send data, but only IO tasks receive any data.
! scounts is number of data elements to send to processor np.
! rcounts is number of data elements to recv from processor np.
! displs is displacement into send array for data to go to proc np
if (nproc <= nanals-1) then
   rcounts = npts_max*ndim
else
   rcounts = 0
endif
do np=0,numproc-1
   displs(np) = np*npts_max*ndim
   if (np <= nanals-1) then
      scounts(np) = npts_max*ndim
   else
      scounts(np) = 0
   end if
enddo
allocate(recvbuf(numproc*npts_max*ndim))
allocate(sendbuf(nanals*npts_max*ndim))
t1 = mpi_wtime()
do nn=1,ndim
 do i=1,numptsperproc(nproc+1)
  do nanal=1,nanals
   n = ((nanal-1)*ndim + (nn-1))*npts_max + i
   ! add ensemble mean back in.
   sendbuf(n) = anal_chunk(nanal,i,nn)+ensmean_chunk(i,nn)
   ! convert to increment (A-F).
   sendbuf(n) = sendbuf(n)-(anal_chunk_prior(nanal,i,nn)+ensmean_chunk_prior(i,nn))
  enddo
 enddo
enddo
call mpi_alltoallv(sendbuf, scounts, displs, mpi_real4, recvbuf, rcounts, displs,&
                   mpi_real4, mpi_comm_world, ierr)
if (nproc <= nanals-1) then
   do np=1,numproc
    do nn=1,ndim
     do i=1,numptsperproc(np)
      n = ((np-1)*ndim + (nn-1))*npts_max + i
      grdin(indxproc(np,i),nn) = recvbuf(n)
     enddo
    enddo
   enddo
   !print *,nproc,'min/max ps',minval(grdin(:,ndim)),maxval(grdin(:,ndim))
end if

if (nproc == 0) then
  t2 = mpi_wtime()
  print *,'time to gather state on root',t2-t1,'secs'
endif

deallocate(sendbuf,recvbuf)
allocate(sendbuf(npts_max*ndim))
allocate(recvbuf(npts*ndim))
if (nproc == 0) t1 = mpi_wtime()
! gather ens. mean anal. increment on root, print out max/mins.
n = 0
do nn=1,ndim
 do i=1,numptsperproc(nproc+1)
   n = n + 1
   ! anal. increment.
   sendbuf(n) = ensmean_chunk(i,nn)-ensmean_chunk_prior(i,nn)
 enddo
enddo
do np=0,numproc-1
   scounts(np) = numptsperproc(np+1)*ndim
   n = 0
   do nn=1,np
      n = n + numptsperproc(nn)*ndim
   enddo
   displs(np) = n
enddo
call mpi_gatherv(sendbuf, numptsperproc(nproc+1)*ndim, mpi_real4, recvbuf, &
      scounts, displs, mpi_real4, 0, mpi_comm_world, ierr)
if (nproc == 0) then
   allocate(ensmean(npts,ndim))
   n = 0
   do np=1,numproc
      do nn=1,ndim
       do i=1,numptsperproc(np)
         n = n + 1
         ensmean(indxproc(np,i),nn) = recvbuf(n)
       enddo
      enddo
   enddo
   if (massbal_adjust) then
     print *,'ens. mean anal. increment min/max ps tend', minval(ensmean(:,ndim-1)),maxval(ensmean(:,ndim-1))
   endif 
   print *,'ens. mean anal. increment min/max ps', minval(ensmean(:,ndim)),maxval(ensmean(:,ndim))
   do nvar=1,nvars
      print *,'ens. mean anal. increment min/max var',nvar,minval(ensmean(:,(nvar-1)*nlevs+1:nvar*nlevs)),maxval(ensmean(:,(nvar-1)*nlevs+1:nvar*nlevs))
   enddo
   deallocate(ensmean)
   t2 = mpi_wtime()
   print *,'time to gather ens mean increment on root',t2-t1,'secs'
end if

deallocate(sendbuf,recvbuf)

! write out state vectors on i/o procs.
if (nproc <= nanals-1) then
   nanal = nproc + 1
   t1 = mpi_wtime()
   if (pseudo_rh .and. nvarhumid > 0) then
      ! re-scale normalized spfh with sat. sphf of first guess
      grdin(:,(nvarhumid-1)*nlevs+1:nvarhumid*nlevs) = grdin(:,(nvarhumid-1)*nlevs+1:nvarhumid*nlevs)*qsat
   end if
   call writegriddata(nanal,grdin)
   if (nproc == 0) then
     t2 = mpi_wtime()
     print *,'time in writegriddata on root',t2-t1,'secs'
   endif 
end if


end subroutine write_ensemble

subroutine statevec_cleanup()
! deallocate module-level allocatable arrays.
if (allocated(anal_chunk)) deallocate(anal_chunk)
if (allocated(anal_chunk_prior)) deallocate(anal_chunk_prior)
if (allocated(ensmean_chunk)) deallocate(ensmean_chunk)
if (allocated(ensmean_chunk_prior)) deallocate(ensmean_chunk_prior)
if (nproc <= nanals-1 .and. allocated(grdin)) deallocate(grdin)
if (nproc <= nanals-1 .and. allocated(qsat)) deallocate(qsat)
deallocate(displs,scounts,rcounts)
end subroutine statevec_cleanup

end module statevec
