module wrf_nmm_interpolate

  !-------------------------------------------------------------------
  !$$$   module documentation block
  !               
  ! module: Interpolate selected model fields from one nmm grid
  !         to another nmm grid, which can have different resolution 
  !         and domain center
  ! 
  ! PRGMMR: Mingjing Tong                          DATE: 2016-02-17
  !
  !$$$ end documentation block
  !-------------------------------------------------------------------

  use kinds
  use namelist
  use wrf_nmm_io

  implicit none

! set default to private
  private
! set subroutines to public
  public :: interpolate_s2d

  !-----------------------------------------------------------------------

contains

  !=======================================================================

  subroutine interpolate_s2d(grd_src,var_src,grd_dst,var_dst)

    !-------------------------------------------------------------------
    !$$$  subprogram documentation block
    ! 
    ! subprogram:   interpolate_s2d 
    ! PRGMMR: Mingjing Tong
    ! abstract: interpolation between nmm domains
    !
    !$$$ end documentation block 
    !-------------------------------------------------------------------

    use wrf_nmm_io, only : wrfnmm_variables
    use general_sub2grid_mod, only: sub2grid_info,general_sub2grid_create_info
    use general_sub2grid_mod, only: general_sub2grid,general_grid2sub
    use general_tll2xy_mod, only: llxy_cons
    use egrid2agrid_mod, only: egrid2agrid_parm,destroy_egrid2agrid
    use constants, only: zero,i_missing
    use mpimod, only: mpi_comm_world,ierror,mype

    implicit none

    type(sub2grid_info), intent(in) :: grd_src
    type(wrfnmm_variables), intent(in) :: var_src
    type(sub2grid_info), intent(inout) :: grd_dst
    type(wrfnmm_variables), intent(inout) :: var_dst

    type(llxy_cons) gt_s,gt_d
    type(egrid2agrid_parm) p_e2a
    integer(i_kind) :: nord_e2a,n3d_half
    integer(i_kind) :: i,j,k,ii,kx,nn,n
    integer(i_kind) :: kpd,ktsk,ksst,ksno,ksice
    integer(i_kind) :: kt,kq,ku,kv,kcwm,kf_ice,kf_rain,kf_rimef,ksmc,kstc
    integer(i_kind) :: npd,ntsk,nsst,nsno,nsice
    integer(i_kind) :: nt,nq,nu,nv,ncwm,nf_ice,nf_rain,nf_rimef,nsmc,nstc
    real(r_kind),allocatable,dimension(:,:,:,:) :: fields_s,fields_d
    real(r_kind),allocatable,dimension(:,:,:,:) :: fields_s2d,work_sub
    real(r_kind),allocatable,dimension(:,:,:,:) :: fields_subs2d
    real(r_single),allocatable,dimension(:,:) :: outwork

    kpd=i_missing ; ktsk=i_missing; ksst=i_missing
    ksno=i_missing ; ksice=i_missing
    ku=i_missing ; kv=i_missing ; kt=i_missing ; kq=i_missing
    kcwm=i_missing ; kf_ice=i_missing ; kf_rain=i_missing
    kf_rimef=i_missing ; ksmc=i_missing ; kstc=i_missing

    nord_e2a=4
    call merge_grid_e_to_grid_a_initialize(var_src%region_lat,     &
         var_src%region_lon,var_dst%region_lat,var_dst%region_lon, &
         grd_src%nlat,grd_src%nlon,grd_dst%nlat,grd_dst%nlon,      &
         nord_e2a,nord_blend,nmix,gt_s,gt_d,p_e2a)

    if(mype == 0)print *,'p_e2a%identity=', p_e2a%identity

    if( .not. p_e2a%identity)then
       if(iblend == 2)then
          call redefineblend(grd_dst%nlat,grd_dst%nlon, &
               var_dst%region_lat,var_dst%region_lon,p_e2a%blend)
       end if

!       if(mype == 0)then
!          allocate(outwork(grd_dst%nlon,grd_dst%nlat))
!          outwork=zero
!          ii=0
!          do j=1,grd_dst%nlon
!             do i=1,grd_dst%nlat
!                ii=ii+1
!                outwork(j,i)=p_e2a%blend(ii)
!             end do
!          end do
!          call outgrads1(outwork,grd_dst%nlon,grd_dst%nlat,'pblend')
!          deallocate(outwork)
!       end if

       do n=1,n3d
          nn=mod(n+1,2)+1
          kx=(n-1)/2
          kx=kx*grd_src%nsig
          select case (trim(var3d(n)))
             case('U')
                ku=kx
                nu=nn
             case('V')
                kv=kx
                nv=nn
             case('T')
                kt=kx
                nt=nn
             case('Q')
                kq=kx
                nq=nn
             case('CWM')
                kcwm=kx
                ncwm=nn
             case('F_ICE')
                kf_ice=kx
                nf_ice=nn
             case('F_RAIN')
                kf_rain=kx
                nf_rain=nn
             case('F_RIMEF')
                kf_rimef=kx
                nf_rimef=nn
             case('SMC')
                ksmc=kx
                nsmc=nn
             case('STC')
                kstc=kx
                nstc=nn
          end select
       end do

       n3d_half=(n3d+1)/2
       kx=n3d_half*grd_src%nsig
       do n=1,n2d
          nn=mod(n+1,2)+1
          kx=kx+n
          select case (trim(var2d(n)))
             case('PD')
                kpd=kx
                npd=nn
             case('TSK')
                ktsk=kx
                ntsk=nn
             case('SST')
                ksst=kx
                nsst=nn
             case('SNO')
                ksno=kx
                nsno=nn
             case('SICE')
                ksice=kx
                nsice=nn
          end select
       end do

       allocate(work_sub(grd_src%inner_vars,grd_src%lat2,grd_src%lon2,grd_src%num_fields))
       allocate(fields_s(grd_src%inner_vars,grd_src%nlat,grd_src%nlon,&
                   grd_src%kbegin_loc:grd_src%kend_alloc))
       work_sub=zero
       do k=1,grd_src%nsig
          do j=1,grd_src%lon2
             do i=1,grd_src%lat2
                work_sub(nu,i,j,ku+k)=var_src%u(i,j,k)       
                work_sub(nv,i,j,kv+k)=var_src%v(i,j,k)
                work_sub(nt,i,j,kt+k)=var_src%t(i,j,k)
                work_sub(nq,i,j,kq+k)=var_src%q(i,j,k)
                if(kcwm > 0)work_sub(ncwm,i,j,kcwm+k)=var_src%cwm(i,j,k)
                if(kf_ice > 0)work_sub(nf_ice,i,j,kf_ice+k)=var_src%f_ice(i,j,k)
                if(kf_rain > 0)work_sub(nf_rain,i,j,kf_rain+k)=var_src%f_rain(i,j,k)
                if(kf_rimef > 0)work_sub(nf_rimef,i,j,kf_rimef+k)=var_src%f_rimef(i,j,k)
                if(ksmc > 0)work_sub(nsmc,i,j,ksmc+k)=var_src%smc(i,j,k)
                if(kstc > 0)work_sub(nstc,i,j,kstc+k)=var_src%stc(i,j,k)
             end do
          end do
       end do
       do j=1,grd_src%lon2
          do i=1,grd_src%lat2
             work_sub(npd,i,j,kpd)=var_src%pd(i,j)
             if(ktsk > 0)work_sub(ntsk,i,j,ktsk)=var_src%tsk(i,j)
             if(ksst > 0)work_sub(nsst,i,j,ksst)=var_src%sst(i,j)
             if(ksno > 0)work_sub(nsno,i,j,ksno)=var_src%sno(i,j)
          end do
       end do
       call general_sub2grid(grd_src,work_sub,fields_s)
       deallocate(work_sub)

       allocate(work_sub(grd_dst%inner_vars,grd_dst%lat2,grd_dst%lon2,grd_dst%num_fields))
       allocate(fields_d(grd_dst%inner_vars,grd_dst%nlat,grd_dst%nlon,&
                   grd_dst%kbegin_loc:grd_dst%kend_alloc))
       work_sub=zero
       do k=1,grd_dst%nsig
          do j=1,grd_dst%lon2
             do i=1,grd_dst%lat2
                work_sub(nu,i,j,ku+k)=var_dst%u(i,j,k)
                work_sub(nv,i,j,kv+k)=var_dst%v(i,j,k)
                work_sub(nt,i,j,kt+k)=var_dst%t(i,j,k)
                work_sub(nq,i,j,kq+k)=var_dst%q(i,j,k)
                if(kcwm > 0)work_sub(ncwm,i,j,kcwm+k)=var_dst%cwm(i,j,k)
                if(kf_ice > 0)work_sub(nf_ice,i,j,kf_ice+k)=var_dst%f_ice(i,j,k)
                if(kf_rain > 0)work_sub(nf_rain,i,j,kf_rain+k)=var_dst%f_rain(i,j,k)
                if(kf_rimef > 0)work_sub(nf_rimef,i,j,kf_rimef+k)=var_dst%f_rimef(i,j,k)
                if(ksmc > 0)work_sub(nsmc,i,j,ksmc+k)=var_dst%smc(i,j,k)
                if(kstc > 0)work_sub(nstc,i,j,kstc+k)=var_dst%stc(i,j,k)
             end do
          end do
       end do
       do j=1,grd_dst%lon2
          do i=1,grd_dst%lat2
             work_sub(npd,i,j,kpd)=var_dst%pd(i,j)
             if(ktsk > 0)work_sub(ntsk,i,j,ktsk)=var_dst%tsk(i,j)
             if(ksst > 0)work_sub(nsst,i,j,ksst)=var_dst%sst(i,j)
             if(ksno > 0)work_sub(nsno,i,j,ksno)=var_dst%sno(i,j)
          end do
       end do
       call general_sub2grid(grd_dst,work_sub,fields_d)
       deallocate(work_sub)

       allocate(fields_s2d(grd_dst%inner_vars,grd_dst%nlat,grd_dst%nlon,&
                   grd_dst%kbegin_loc:grd_dst%kend_alloc))
       fields_s2d=zero
       do k=grd_dst%kbegin_loc,grd_dst%kend_alloc
          if(grd_dst%vector(k))then
             call merge_vgrid_e_to_vgrid_a(fields_s(1,:,:,k), &
                  fields_s(2,:,:,k),fields_d(1,:,:,k),        &
                  fields_d(2,:,:,k),fields_s2d(1,:,:,k),      &
                  fields_s2d(2,:,:,k),gt_s,gt_d,p_e2a)
          else
             do ii=1,grd_dst%inner_vars
                call merge_grid_e_to_grid_a(fields_s(ii,:,:,k), &
                     fields_d(ii,:,:,k),fields_s2d(ii,:,:,k),   &
                     gt_s,gt_d,p_e2a)
             end do
          end if
       end do

       deallocate(fields_s,fields_d)
       call destroy_egrid2agrid(p_e2a)

       allocate(fields_subs2d(grd_dst%inner_vars,grd_dst%lat2,grd_dst%lon2,grd_dst%num_fields))
       call general_grid2sub(grd_dst,fields_s2d,fields_subs2d)
       deallocate(fields_s2d)

       do k=1,grd_dst%nsig
          do j=1,grd_dst%lon2
             do i=1,grd_dst%lat2
                var_dst%u(i,j,k)=fields_subs2d(nu,i,j,ku+k)
                var_dst%v(i,j,k)=fields_subs2d(nv,i,j,kv+k)
                var_dst%t(i,j,k)=fields_subs2d(nt,i,j,kt+k)
                var_dst%q(i,j,k)=fields_subs2d(nq,i,j,kq+k)
                if(kcwm > 0)var_dst%cwm(i,j,k)=fields_subs2d(ncwm,i,j,kcwm+k)
                if(kf_ice > 0)var_dst%f_ice(i,j,k)=fields_subs2d(nf_ice,i,j,kf_ice+k)
                if(kf_rain > 0)var_dst%f_rain(i,j,k)=fields_subs2d(nf_rain,i,j,kf_rain+k)
                if(kf_rimef > 0)var_dst%f_rimef(i,j,k)=fields_subs2d(nf_rimef,i,j,kf_rimef+k)
                if(ksmc > 0)var_dst%smc(i,j,k)=fields_subs2d(nsmc,i,j,ksmc+k)
                if(kstc > 0)var_dst%stc(i,j,k)=fields_subs2d(nstc,i,j,kstc+k)
             end do
          end do
       end do
       do j=1,grd_dst%lon2
          do i=1,grd_dst%lat2
             var_dst%pd(i,j)=fields_subs2d(npd,i,j,kpd)
             if(ktsk > 0)var_dst%tsk(i,j)=fields_subs2d(ntsk,i,j,ktsk)
             if(ksst > 0)var_dst%sst(i,j)=fields_subs2d(nsst,i,j,ksst)
             if(ksno > 0)var_dst%sno(i,j)=fields_subs2d(nsno,i,j,ksno)
             if(ksice > 0)var_dst%sice(i,j)=fields_subs2d(nsice,i,j,ksice)
          end do
       end do
       deallocate(fields_subs2d)
    else
       var_dst%u=var_src%u
       var_dst%v=var_src%v
       var_dst%t=var_src%t
       var_dst%q=var_src%q
       if(kcwm > 0)var_dst%cwm=var_src%cwm
       if(kf_ice > 0)var_dst%f_ice=var_src%f_ice
       if(kf_rain > 0)var_dst%f_rain=var_src%f_rain
       if(kf_rimef > 0)var_dst%f_rimef=var_src%f_rimef
       if(ksmc > 0)var_dst%smc=var_src%smc
       if(kstc > 0)var_dst%stc=var_src%stc
       var_dst%pd=var_src%pd
       if(ktsk > 0)var_dst%tsk=var_src%tsk
       if(ksst > 0)var_dst%sst=var_src%sst
       if(ksno > 0)var_dst%sno=var_src%sno
       if(ksice > 0)var_dst%sice=var_src%sice
    end if

  end subroutine interpolate_s2d

  subroutine redefineblend(nlat,nlon,region_lat,region_lon,pblend)

    !-------------------------------------------------------------------
    !$$$  subprogram documentation block
    !
    ! subprogram:  redefineblend
    ! PRGMMR: Mingjing Tong
    ! abstract: create a circled blending zone around storm center
    !
    !$$$ end documentation block
    !-------------------------------------------------------------------

    use constants, only: zero, one, rad2deg
    use blendmod, only: blend
    use mpimod, only: mype

    implicit none

    integer(i_kind), intent(in) :: nlat
    integer(i_kind), intent(in) :: nlon
    real(r_kind), intent(in) :: region_lat(nlat,nlon)
    real(r_kind), intent(in) :: region_lon(nlat,nlon)
    real(r_kind), intent(out) :: pblend(nlat*nlon)

    integer(i_kind) :: id_storm,iclat,iclon,ipsfc
    integer(i_kind) :: ipcls,irmax,ivobs,ir_vobs
    integer(i_kind) :: ir_v4(4)
    integer(i_kind) :: i,j,ii,mm,n
    integer(i_kind),dimension(0:40):: iblend
    character(1) :: sn,ew,depth
    real(r_kind) :: stormrad,rad1,rad2
    real(r_kind) :: clat,clon
    real(r_kind) :: dist,dist1,dist2,x,y

    open(85,file='storm_radius',form='unformatted')
    read(85)stormrad
    rad2=min(stormrad+5.0_r_kind,radmax)
    rad1=rad2*(one-bwdth)
    if(mype == 0)print *,'stormrad,rad1,rad2=', stormrad,rad1,rad2 
    close(85)

    open(11,file='tcvitals.as',form='formatted')
    read(11,11) id_storm,iclat,sn,iclon,ew,ipsfc,ipcls,  &
                  irmax,ivobs,ir_vobs,(ir_v4(i),i=1,4),depth

 11   format(5x,I2,26x,I3,A1,I5,A1,9x,I4,1x,I4,1x,I4,I3,I4,4I5,1x,A1)
    close(11)
    clat=iclat*0.1_r_kind
    clon=iclon*0.1_r_kind

    if(sn=='S')clat=-clat
    if(ew=='W')clon=-clon

    dist1=gc_dist(clat+rad1,clon,clat,clon)
    dist2=gc_dist(clat+rad2,clon,clat,clon)

    if(mype == 0)print *,'clat,clon=',clat,clon
    if(mype == 0)print *,'dist1,dist2=',dist1,dist2

    pblend=zero

    mm=nord_blend
    call blend(mm,iblend)

    ii=0
    do j=1,nlon
       do i=1,nlat
          ii=ii+1
          dist=gc_dist(region_lat(i,j)*rad2deg,region_lon(i,j)*rad2deg,clat,clon)
          if(dist <= dist1)then
             pblend(ii)=one
          else if(dist < dist2)then
             x=one-(dist-dist1)/(dist2-dist1)
             y=iblend(mm)
             do n=mm-1,0,-1
                y=x*y+iblend(n)
             end do
             y=y*x**(mm+1)
             pblend(ii)=y
          else
             pblend(ii)=zero
          end if
       end do
    end do
 
  end subroutine redefineblend

  real(r_kind) function gc_dist(inlat1,inlon1,inlat2,inlon2)
  !$$$  subprogram documentation block
  !                .      .    .                                       .
  ! subprogram:    gc_dist
  !   prgmmr:                  org:                     date:
  !
  ! abstract: Return the great circle distance (m) between a two pairs of lat/lon points
  !
  ! program history log:
  !   2014-09-09  carley - added subprogram doc block
  !
  !   input argument list:
  !   lat1,lon1,lat2,lon2 in degrees
  !
  !   output argument list:
  !    gc_dist
  !
  ! attributes:
  !   language: f90
  !   machine:
  !
  !$$$ end documentation block

  use constants, only: rearth,deg2rad,one,two
  use kinds, only: r_kind
  implicit none
  real(r_kind),intent (in) :: inlat1,inlon1,inlat2,inlon2
  real(r_kind) :: lat1,lon1,lat2,lon2
  real(r_kind) :: dLat,dLon,a,c

  lat2=inlat2*deg2rad
  lat1=inlat1*deg2rad
  lon1=inlon1*deg2rad
  lon2=inlon2*deg2rad
  dLat = lat2 - lat1
  dLon = lon2 - lon1
  a = sin(dLat / two) * sin(dLat / two) +  cos(lat1) * cos(lat2) * &
      sin(dLon / two) * sin(dLon / two)
  c = two * atan2(sqrt(a), sqrt(one - a))
  gc_dist = rearth * c

  end function gc_dist

subroutine outgrads1(f,nx,ny,label)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    outgrads1
!   prgmmr:
!
! abstract:
!
! program history log:
!   2009-09-18  lueken - added subprogram doc block
!   2012-12-11  parrish - assign np a value.
!
!   input argument list:
!    label
!    nx,ny
!    f
!
!   output argument list:
!
! attributes:
!   language: f90
!   machine:
!
!$$$ end documentation block
  use kinds, only: i_kind,r_single
  implicit none

  character(*)   ,intent(in   ) :: label
  integer(i_kind),intent(in   ) :: nx,ny
  real(r_single) ,intent(in   ) :: f(nx,ny)

  integer(i_kind) i,l,next,last,np,ntime,ioutdat,ioutcor,koutmax
  real(r_single) rlonmap0,undef,dlonmap,pinc,startp,rlatmap0,dlatmap
  character(80) dsdes,dsdat
  character(80) datdes(1000)
  character(1) blank
  data blank/' '/
  data undef/-9.99e33_r_single/

  ioutcor=10
  ioutdat=11
  np=1

  write(dsdes,'(a,".des")')trim(label)
  write(dsdat,'(a,".dat")')trim(label)
  open(unit=ioutcor,file=dsdes,form='formatted')
  open(unit=ioutdat,file=dsdat,form='unformatted')
  ntime=1
  rlonmap0=1._r_single
  dlonmap=1._r_single
  rlatmap0=1._r_single
  dlatmap=1._r_single
  startp=1._r_single
  pinc=1._r_single
  koutmax=1
  do i=1,1000
     write(datdes(i),'(80a1)')(blank,l=1,80)
  end do
  write(datdes(1),'("DSET ",a)')trim(dsdat)
  write(datdes(2),'("options big_endian sequential")')
  write(datdes(3),'("TITLE ",a)')trim(label)
  write(datdes(4),'("UNDEF ",e11.2)')undef
  write(datdes(5),'("XDEF ",i5," LINEAR ",f7.2,f7.2)')nx,rlonmap0,dlonmap
  write(datdes(6),'("YDEF ",i5," LINEAR ",f7.2,f7.2)')ny,rlatmap0,dlatmap
  next=7
  write(datdes(next),'("ZDEF ",i5," LINEAR ",f7.2,f7.2)')np,startp,pinc
  next=next+1
  write(datdes(next),'("TDEF ",i5," LINEAR 0Z23may1992 24hr")')koutmax
  next=next+1
  write(datdes(next),'("VARS 1")')
  next=next+1
  write(datdes(next),'("f   ",i5," 99 f   ")')np
  next=next+1
  write(datdes(next),'("ENDVARS")')
  last=next
  write(ioutcor,'(a80)')(datdes(i),i=1,last)
  close(ioutcor)

  write(ioutdat) f
  close(ioutdat)

end subroutine outgrads1

end module wrf_nmm_interpolate
