module wrf_nmm_io

  !-------------------------------------------------------------------
  !$$$   module documentation block
  !               
  ! module: read and write nmm variables from and to wrf input/output
  !         file
  !   
  ! PRGMMR: Mingjing Tong                          DATE: 2016-02-17
  !
  !$$$ end documentation block
  !-------------------------------------------------------------------

  use kinds
  use netcdf
  use namelist

  implicit none
! set default to private
  private
! set subroutines to public
  public :: wrfnmm_variables
  public :: wrf_nmm_interface, wrf_nmm_read, wrf_nmm_write,grads3d

  type wrfnmm_variables
     real(r_kind)                                                  :: pt
     real(r_kind)                                                  :: pdtop
     real(r_kind),                 dimension(:),       allocatable :: eta1
     real(r_kind),                 dimension(:),       allocatable :: eta2
     real(r_kind),                 dimension(:,:),     allocatable :: region_lon
     real(r_kind),                 dimension(:,:),     allocatable :: region_lat
     real(r_kind),                 dimension(:,:),     allocatable :: pd
     real(r_kind),                 dimension(:,:),     allocatable :: tsk
     real(r_kind),                 dimension(:,:),     allocatable :: sst
     real(r_kind),                 dimension(:,:),     allocatable :: sno
     real(r_kind),                 dimension(:,:),     allocatable :: sice
     real(r_kind),                 dimension(:,:,:),   allocatable :: pint
     real(r_kind),                 dimension(:,:,:),   allocatable :: t
     real(r_kind),                 dimension(:,:,:),   allocatable :: q
     real(r_kind),                 dimension(:,:,:),   allocatable :: u
     real(r_kind),                 dimension(:,:,:),   allocatable :: v
     real(r_kind),                 dimension(:,:,:),   allocatable :: cwm
     real(r_kind),                 dimension(:,:,:),   allocatable :: f_ice
     real(r_kind),                 dimension(:,:,:),   allocatable :: f_rain
     real(r_kind),                 dimension(:,:,:),   allocatable :: f_rimef
     real(r_kind),                 dimension(:,:,:),   allocatable :: smc
     real(r_kind),                 dimension(:,:,:),   allocatable :: stc
  end type wrfnmm_variables

  ! Define global variables

  integer(i_kind) :: ncfileid
  integer(i_kind) :: ncvarid
  integer(i_kind) :: ncdimid
  integer(i_kind) :: ncstatus
  integer(i_kind) :: iunit
  
  data iunit / 15 /

contains


  subroutine wrf_nmm_interface(filename)

    !-------------------------------------------------------------------
    !$$$  subprogram documentation block
    ! 
    ! subprogram:   wrf_nmm_interface
    ! PRGMMR: Mingjing Tong
    ! abstract: read model variables from wrf netcdf file
    !           and write the result to temporary file expected by
    !           wrf_nmm_read
    !
    !$$$ end documentation block 
    !-------------------------------------------------------------------

    implicit none

    character(120), intent(in) :: filename

    integer(i_kind) :: nlon,nlat,nsig,i,j,k,n
    real(r_single) :: pt,pdtop
    real(r_single),allocatable :: eta1(:), eta2(:)
    real(r_single),allocatable :: f2d(:,:),f3d(:,:,:)

    print *,'filename=',trim(filename)
    ncstatus = nf90_open(path=trim(filename),            &
         & mode=nf90_nowrite,ncid=ncfileid)
    ncstatus = nf90_inq_dimid(ncfileid,'west_east',ncdimid)
    ncstatus = nf90_inquire_dimension(ncfileid,ncdimid,                    &
         & len=nlon)
    ncstatus = nf90_inq_dimid(ncfileid,'south_north',ncdimid)
    ncstatus = nf90_inquire_dimension(ncfileid,ncdimid,                    &
         & len=nlat)
    ncstatus = nf90_inq_dimid(ncfileid,'bottom_top',ncdimid)
    ncstatus = nf90_inquire_dimension(ncfileid,ncdimid,                    &
         & len=nsig)

    open(iunit,file=trim(filename)//'_tmp',form='unformatted')
    write(iunit)nlon,nlat,nsig
    ncstatus = nf90_inq_varid(ncfileid,'PT',ncvarid)
    ncstatus = nf90_get_var(ncfileid,ncvarid,pt)
    ncstatus = nf90_inq_varid(ncfileid,'PDTOP',ncvarid)
    ncstatus = nf90_get_var(ncfileid,ncvarid,pdtop)
    write(iunit)pt,pdtop
    allocate(eta1(nsig+1),eta2(nsig+1))
    ncstatus = nf90_inq_varid(ncfileid,'ETA1',ncvarid)
    ncstatus = nf90_get_var(ncfileid,ncvarid,eta1)
    ncstatus = nf90_inq_varid(ncfileid,'ETA2',ncvarid)
    ncstatus = nf90_get_var(ncfileid,ncvarid,eta2)
    write(iunit)eta1,eta2
    deallocate(eta1,eta2)

    allocate(f2d(nlon,nlat))

    ncstatus = nf90_inq_varid(ncfileid,'GLAT',ncvarid)
    ncstatus = nf90_get_var(ncfileid,ncvarid,f2d)
    write(iunit)f2d

    ncstatus = nf90_inq_varid(ncfileid,'GLON',ncvarid)
    ncstatus = nf90_get_var(ncfileid,ncvarid,f2d)
    write(iunit)f2d

    do n=1,n2d
       print *,'2d variable=', trim(var2d(n))
       ncstatus = nf90_inq_varid(ncfileid,trim(var2d(n)),ncvarid)
       ncstatus = nf90_get_var(ncfileid,ncvarid,f2d)
    write(iunit)f2d
    end do
    deallocate(f2d)
    allocate(f3d(nlon,nlat,nsig))
    do n=1,n3d
       print *,'3d variable ', trim(var3d(n))
       ncstatus = nf90_inq_varid(ncfileid,trim(var3d(n)),ncvarid)
       ncstatus = nf90_get_var(ncfileid,ncvarid,f3d)
       do k=1,nsig
          write(iunit)((f3d(i,j,k),i=1,nlon),j=1,nlat) 
       end do
    end do
    deallocate(f3d)
    allocate(f3d(nlon,nlat,nsig+1))
    ncstatus = nf90_inq_varid(ncfileid,'PINT',ncvarid)
    ncstatus = nf90_get_var(ncfileid,ncvarid,f3d)
    do k=1,nsig+1
        write(iunit)((f3d(i,j,k),i=1,nlon),j=1,nlat)   ! T
    end do
    deallocate(f3d)
    ncstatus = nf90_close(ncfileid)
    close(iunit)
    print *,'interface completed'

  end subroutine wrf_nmm_interface

  subroutine wrf_nmm_read(filename,grd,vrbls)

    !-------------------------------------------------------------------
    !$$$  subprogram documentation block
    ! 
    ! subprogram:   wrf_nmm_read
    ! PRGMMR: Mingjing Tong
    ! abstract: read model variables and distribute to subdomains 
    !
    !$$$ end documentation block 
    !-------------------------------------------------------------------

    use general_sub2grid_mod, only: sub2grid_info,general_sub2grid_create_info
    use mpimod, only: ierror,mpi_real4,mpi_comm_world,npe,mype
    use constants, only: i_missing,pi

    implicit none

    character(120), intent(in) :: filename
    type(sub2grid_info), intent(out) :: grd
    type(wrfnmm_variables), intent(out) :: vrbls
    integer(i_kind) :: nlon, nlat, nsig, num_nmm_fields
    integer(i_kind) :: nlat_a, nlon_a
    logical :: regional
    logical,allocatable :: vector(:)
    real(r_single) :: pt,pdtop
    real(r_single),allocatable :: eta1(:),eta2(:)
    real(r_single),allocatable :: f2d(:,:)
    real(r_kind),allocatable :: glat(:,:),glon(:,:)
    real(r_kind),allocatable :: glat_an(:,:),glon_an(:,:)
    real(r_kind),allocatable :: gxtemp(:,:),gytemp(:,:)
    real(r_kind),allocatable :: gxtemp_an(:,:),gytemp_an(:,:)
    real(r_single),allocatable :: tempa(:)
    real(r_single),allocatable :: temp1(:,:)
    real(r_single),allocatable :: all_loc(:,:,:)
    integer(i_kind),allocatable :: igtype(:)
    integer(i_kind) :: irc_s_reg(npe),ird_s_reg(npe)
    integer(i_kind) :: inner_vars, n3d_half, n2d_half, num_fields
    integer(i_kind) num_loc_groups,num_all_pad
    integer(i_kind) i,j,k,n,ifld,icount,icount_prev
    integer(i_kind) i0,j0
    integer(i_kind) :: kpd,ktsk,ksst,ksno,ksice
    integer(i_kind) :: kt,kq,ku,kv,kcwm,kf_ice,kf_rain,kf_rimef,ksmc,kstc
    real(r_kind) rotate3

    kpd=i_missing ; ktsk=i_missing; ksst=i_missing 
    ksno=i_missing ; ksice=i_missing 
    ku=i_missing ; kv=i_missing ; kt=i_missing ; kq=i_missing
    kcwm=i_missing ; kf_ice=i_missing ; kf_rain=i_missing 
    kf_rimef=i_missing ; ksmc=i_missing ; kstc=i_missing

    open(iunit,file=trim(filename)//'_tmp',form='unformatted')

    read(iunit)nlon,nlat,nsig

    num_nmm_fields=max(0,n3d)*nsig+max(0,n2d)
    num_loc_groups=num_nmm_fields/npe
    if(e2a)then
       nlon_a=2*nlon-1
       nlat_a=nlat
    else
       nlon_a=nlon
       nlat_a=nlat
    end if

    if(mype == 0)print *,'nlon,nlat,nlon_a,nlat_a,nsig=',nlon,nlat,nlon_a,nlat_a,nsig

    regional=.true.
    inner_vars=2
    n3d_half=(n3d+1)/2
    n2d_half=(n2d+1)/2
    num_fields=max(0,n3d_half)*nsig+max(0,n2d_half)
    allocate(vector(num_fields))
    vector=.false.
    vector(1:nsig)=.true.
    call general_sub2grid_create_info(grd,inner_vars,nlat_a,nlon_a,nsig, &
                                      num_fields,regional,vector)
    deallocate(vector)

    read(iunit)pt,pdtop
    vrbls%pt=pt
    vrbls%pdtop=pdtop

    allocate(eta1(nsig+1),eta2(nsig+1))
    read(iunit)eta1,eta2

    if(.not. allocated(vrbls%eta1))                                         &
         & allocate(vrbls%eta1(nsig+1))
    if(.not. allocated(vrbls%eta2))                                         &
         & allocate(vrbls%eta2(nsig+1))
    vrbls%eta1=eta1
    vrbls%eta2=eta2
    deallocate(eta1,eta2)

    allocate(f2d(nlon,nlat))
    allocate(glat(nlon,nlat),glon(nlon,nlat))
    read(iunit)f2d
    glat=f2d
    read(iunit)f2d
    glon=f2d

    if(.not. allocated(vrbls%region_lat))                                   &
         & allocate(vrbls%region_lat(nlat_a,nlon_a))
    if(.not. allocated(vrbls%region_lon))                                   &
         & allocate(vrbls%region_lon(nlat_a,nlon_a))

    if(e2a) then
       rotate3=pi/4._r_kind
       allocate(gxtemp(nlon,nlat),gytemp(nlon,nlat))
       i0=nlon/2
       j0=nlat/2
       call ll2rpolar(glat,glon,nlon*nlat, &
                      gxtemp,gytemp,glat(i0,j0),glon(i0,j0),rotate3)
       allocate(gxtemp_an(nlon_a,nlat_a),gytemp_an(nlon_a,nlat_a))
       call fill_nmm_grid2a3(gxtemp,nlon,nlat,gxtemp_an)
       call fill_nmm_grid2a3(gytemp,nlon,nlat,gytemp_an)
       allocate(glat_an(nlon_a,nlat_a),glon_an(nlon_a,nlat_a))
       call rpolar2ll(gxtemp_an,gytemp_an,nlon_a*nlat_a, &
                      glat_an,glon_an,glat(i0,j0),glon(i0,j0),rotate3)
       do k=1,nlon_a
          do i=1,nlat_a
             vrbls%region_lat(i,k)=glat_an(k,i)
             vrbls%region_lon(i,k)=glon_an(k,i)
          end do
       end do

       deallocate(gxtemp,gytemp,gxtemp_an,gytemp_an,glat_an,glon_an)
    else
       do k=1,nlon_a
          do i=1,nlat_a
             vrbls%region_lat(i,k)=glat(k,i)
             vrbls%region_lon(i,k)=glon(k,i)
          end do
       end do
    end if
    deallocate(glon,glat)

    do
       num_all_pad=num_loc_groups*npe
       if(num_all_pad >= num_nmm_fields) exit
       num_loc_groups=num_loc_groups+1
    end do

    allocate(all_loc(grd%lat2,grd%lon2,num_all_pad))
    allocate(igtype(num_nmm_fields))
    igtype=1

    i=1
    do n=1,n2d
       select case (trim(var2d(n)))
          case('PD')
             kpd=i
          case('TSK')
             ktsk=i
          case('SST')
             ksst=i
          case('SNO')
             ksno=i
          case('SICE')
             ksice=i
       end select
       i=i+1
    end do
       
    do n=1,n3d
       select case (trim(var3d(n)))
          case('U')
             ku=i
             igtype(ku)=2
          case('V')
             kv=i
             igtype(kv)=2
          case('T')
             kt=i
          case('Q')
             kq=i
          case('CWM')
             kcwm=i
          case('F_ICE')
             kf_ice=i
          case('F_RAIN')
             kf_rain=i
          case('F_RIMEF')
             kf_rimef=i
          case('SMC')
             ksmc=i
          case('STC')
             kstc=i
       end select
       i=i+grd%nsig
    end do

    do k=1,grd%nsig-1
       i=ku+k                      
       igtype(i)=2
       j=kv+k
       igtype(j)=2
    end do

    do i=1,npe
       irc_s_reg(i)=grd%ijn_s(mype+1)
    end do
    ird_s_reg(1)=0
    do i=1,npe
       if(i /= 1) ird_s_reg(i)=ird_s_reg(i-1)+irc_s_reg(i-1)
    end do

    allocate(temp1(nlon,nlat))
    allocate(tempa(grd%itotsub))
    icount=0
    icount_prev=1
    do ifld=1,num_nmm_fields
       icount=icount+1
       if(mype == mod(icount-1,npe)) then
          read(iunit)((temp1(i,j),i=1,nlon),j=1,nlat)
          if(e2a) call fill_nmm_grid(grd,temp1,nlon,nlat,tempa,abs(igtype(ifld)),1)
       else
          read(iunit)
       end if

       if(mod(icount,npe) == 0.or.icount == num_nmm_fields) then
           call mpi_alltoallv(tempa,grd%ijn_s,grd%displs_s,mpi_real4, &
                all_loc(1,1,icount_prev),irc_s_reg,ird_s_reg,mpi_real4,mpi_comm_world,ierror)
           icount_prev=icount+1
       end if
    end do
    close(iunit)

    if(.not. allocated(vrbls%t) .and. kt > 0)                            &
         & allocate(vrbls%t(grd%lat2,grd%lon2,nsig))
    if(.not. allocated(vrbls%q) .and. kq > 0)                            &
         & allocate(vrbls%q(grd%lat2,grd%lon2,nsig))
    if(.not. allocated(vrbls%u) .and. ku > 0)                            &
         & allocate(vrbls%u(grd%lat2,grd%lon2,nsig))
    if(.not. allocated(vrbls%v) .and. kv > 0)                            &
         & allocate(vrbls%v(grd%lat2,grd%lon2,nsig))
    if(.not. allocated(vrbls%cwm) .and. kcwm > 0)                        &
         & allocate(vrbls%cwm(grd%lat2,grd%lon2,nsig))
    if(.not. allocated(vrbls%f_ice) .and. kf_ice > 0)                    &
         & allocate(vrbls%f_ice(grd%lat2,grd%lon2,nsig))
    if(.not. allocated(vrbls%f_rain) .and. kf_rain > 0)                  &
         & allocate(vrbls%f_rain(grd%lat2,grd%lon2,nsig))
    if(.not. allocated(vrbls%f_rimef) .and. kf_rimef > 0)                &
         & allocate(vrbls%f_rimef(grd%lat2,grd%lon2,nsig))
    if(.not. allocated(vrbls%smc) .and. ksmc > 0)                        &
         & allocate(vrbls%smc(grd%lat2,grd%lon2,nsig))
    if(.not. allocated(vrbls%stc) .and. kstc > 0)                        &
         & allocate(vrbls%stc(grd%lat2,grd%lon2,nsig))

    do k=1,grd%nsig
       do i=1,grd%lon2
          do j=1,grd%lat2
             if(ku>0)vrbls%u(j,i,k) = all_loc(j,i,ku)
             if(kv>0)vrbls%v(j,i,k) = all_loc(j,i,kv)
             if(kt>0)vrbls%t(j,i,k) = all_loc(j,i,kt)
             if(kq>0)vrbls%q(j,i,k) = all_loc(j,i,kq)
             if(kcwm>0)vrbls%cwm(j,i,k) = all_loc(j,i,kcwm)
             if(kf_ice>0)vrbls%f_ice(j,i,k) = all_loc(j,i,kf_ice)
             if(kf_rain>0)vrbls%f_rain(j,i,k) = all_loc(j,i,kf_rain)
             if(kf_rimef>0)vrbls%f_rimef(j,i,k) = all_loc(j,i,kf_rimef)
             if(ksmc>0)vrbls%smc(j,i,k) = all_loc(j,i,ksmc)
             if(kstc>0)vrbls%stc(j,i,k) = all_loc(j,i,kstc)
          end do
       end do
       ku=ku+1
       kv=kv+1
       kt=kt+1
       kq=kq+1
       kcwm=kcwm+1
       kf_ice=kf_ice+1
       kf_rain=kf_rain+1
       kf_rimef=kf_rimef+1
       ksmc=ksmc+1
       kstc=kstc+1
    end do

    if(.not. allocated(vrbls%pd) .and. kpd > 0)                          &
         & allocate(vrbls%pd(grd%lat2,grd%lon2))
    if(.not. allocated(vrbls%tsk) .and. ktsk > 0)                        &
         & allocate(vrbls%tsk(grd%lat2,grd%lon2))
    if(.not. allocated(vrbls%sst) .and. ksst > 0)                        &
         & allocate(vrbls%sst(grd%lat2,grd%lon2))
    if(.not. allocated(vrbls%sno) .and. ksno > 0)                        &
         & allocate(vrbls%sno(grd%lat2,grd%lon2))
    if(.not. allocated(vrbls%sice) .and. ksice > 0)                      &
         & allocate(vrbls%sice(grd%lat2,grd%lon2))

    do i=1,grd%lon2
       do j=1,grd%lat2
          if(kpd>0)vrbls%pd(j,i)=all_loc(j,i,kpd)
          if(ktsk>0)vrbls%tsk(j,i)=all_loc(j,i,ktsk)
          if(ksst>0)vrbls%sst(j,i)=all_loc(j,i,ksst)
          if(ksno>0)vrbls%sno(j,i)=all_loc(j,i,ksno)
          if(ksice>0)vrbls%sice(j,i)=all_loc(j,i,ksice)
       end do
    end do

    deallocate(f2d)
    deallocate(igtype,temp1,tempa,all_loc)

    !call grads3d(grd,vrbls%u,vrbls%v,vrbls%t,vrbls%q,vrbls%pd, &
    !             grd%nsig,mype,trim(filename)//'_tmp')

    if(mype == 0)print *,'wrf_nmm_read completed'

  end subroutine wrf_nmm_read

  subroutine wrf_nmm_write(filename,grd,vrbls)

    !-------------------------------------------------------------------
    !$$$  subprogram documentation block
    !
    ! subprogram:   wrf_nmm_write
    ! PRGMMR: Mingjing Tong
    ! abstract: write updated model variables to wrf netcdf file 
    !
    !$$$ end documentation block
    !-------------------------------------------------------------------

    use general_sub2grid_mod, only: sub2grid_info
    use mpimod, only: ierror,mpi_real4,mpi_comm_world,npe,mype
    use constants, only: zero_single

    implicit none

    character(120), intent(in) :: filename
    type(sub2grid_info), intent(in) :: grd
    type(wrfnmm_variables), intent(in) :: vrbls

    integer(i_kind) :: nlon,nlat,nsig
    integer(i_kind) :: igtype
    real(r_single), allocatable :: f2d(:,:),f3d(:,:,:)
    real(r_single), allocatable :: temp(:,:),temp1(:,:),tempa(:),strp(:)
    real(r_single),allocatable::pint(:,:,:),pd(:,:)
    integer(i_kind) i,j,k,n

    !=====================================================================

    ! Check local variable and proceed accordingly

    !call grads3d(grd,vrbls%u,vrbls%v,vrbls%t,vrbls%q,vrbls%pd, &
    !             grd%nsig,mype,'output')

    ncstatus = nf90_open(path=trim(filename),          &
               mode=nf90_nowrite,ncid=ncfileid)
    ncstatus = nf90_inq_dimid(ncfileid,'west_east',ncdimid)
    ncstatus = nf90_inquire_dimension(ncfileid,ncdimid,                    &
         & len=nlon)
    ncstatus = nf90_inq_dimid(ncfileid,'south_north',ncdimid)
    ncstatus = nf90_inquire_dimension(ncfileid,ncdimid,                    &
         & len=nlat)
    ncstatus = nf90_inq_dimid(ncfileid,'bottom_top',ncdimid)
    ncstatus = nf90_inquire_dimension(ncfileid,ncdimid,                    &
         & len=nsig)
    ncstatus = nf90_close(ncfileid)

    if(e2a)then
       if(grd%nlon /= 2*nlon-1 .or. grd%nlat /= nlat .or. &
          grd%nsig /= nsig)then 
          print *,'grid points are not consistant'
          call stop2(111)
       end if
    else
       if(grd%nlon /= nlon .or. grd%nlat /= nlat .or. &
          grd%nsig /= nsig)then
          print *,'grid points are not consistant'
          call stop2(111)
       end if
    end if

    if(mype == 0)then
       print *,'write to file: ',trim(filename)
       ncstatus = nf90_open(path=trim(filename),          &
                  mode=nf90_write,ncid=ncfileid)
    end if

    allocate(tempa(grd%itotsub))
    allocate(temp1(grd%lat2,grd%lon2))
    allocate(strp(grd%lat1*grd%lon1))
    if(mype == 0) then
       allocate(temp(grd%nlon,grd%nlat))
       allocate(f3d(nlon,nlat,nsig))
       allocate(f2d(nlon,nlat))
    end if

    do n=1,n2d
       temp1=zero_single
       igtype=1
       select case (trim(var2d(n)))
          case('PD')
             temp1=vrbls%pd
          case('TSK')
             temp1=vrbls%tsk
          case('SST')
             temp1=vrbls%sst
          case('SNO')
             temp1=vrbls%sno
          case('SICE')
             temp1=vrbls%sice
       end select

       call strip_grd_single(grd,temp1,strp)
       call mpi_gatherv(strp,grd%ijn(mype+1),mpi_real4, &
            tempa,grd%ijn,grd%displs_g,mpi_real4,0,mpi_comm_world,ierror)
       if(mype == 0) then
          if(trim(var2d(n)) == 'PD')then
             allocate(pint(nlon,nlat,nsig+1))
             allocate(pd(nlon,nlat))
             ncstatus = nf90_inq_varid(ncfileid,'PINT',ncvarid)
             ncstatus = nf90_get_var(ncfileid,ncvarid,pint)
             ncstatus = nf90_inq_varid(ncfileid,'PD',ncvarid)
             ncstatus = nf90_get_var(ncfileid,ncvarid,pd)
          end if
          ncstatus = nf90_inq_varid(ncfileid,trim(var2d(n)),ncvarid)
          ncstatus = nf90_get_var(ncfileid,ncvarid,f2d)
          do i=1,grd%iglobal
             temp(grd%ltosj(i),grd%ltosi(i))=tempa(i)
          end do
          if(e2a) then
             call unfill_nmm_grid(temp,nlon,nlat,f2d,igtype)
          else
             f2d=temp
          end if
          ncstatus = nf90_put_var(ncfileid,ncvarid,f2d)
          if(trim(var2d(n)) == 'PD')then
             do k=1,nsig+1
                do j=1,nlat
                   do i=1,nlon 
                      pint(i,j,k)=pint(i,j,k) &
                          +vrbls%eta2(k)*(f2d(i,j)-pd(i,j))                
                   end do
                end do
             end do
             ncstatus = nf90_inq_varid(ncfileid,'PINT',ncvarid)
             ncstatus = nf90_put_var(ncfileid,ncvarid,pint)
             deallocate(pd,pint)
          end if
       end if
    end do

    do n=1,n3d
       temp1=zero_single
       if(mype == 0) then
          ncstatus = nf90_inq_varid(ncfileid,trim(var3d(n)),ncvarid)
          ncstatus = nf90_get_var(ncfileid,ncvarid,f3d)
       end if
       do k=1,nsig
          select case (trim(var3d(n)))
             case('U')
                temp1=vrbls%u(:,:,k); igtype=2
             case('V')
                temp1=vrbls%v(:,:,k); igtype=2
             case('T')
                temp1=vrbls%t(:,:,k); igtype=1
             case('Q')
                temp1=vrbls%q(:,:,k); igtype=1
             case('CWM')
                temp1=vrbls%cwm(:,:,k); igtype=1
             case('F_ICE')
                temp1=vrbls%f_ice(:,:,k); igtype=1
             case('F_RAIN')
                temp1=vrbls%f_rain(:,:,k); igtype=1
             case('F_RIMEF')
                temp1=vrbls%f_rimef(:,:,k); igtype=1
             case('SMC')
                temp1=vrbls%smc(:,:,k); igtype=1
             case('STC')
                temp1=vrbls%stc(:,:,k); igtype=1
          end select
          call strip_grd_single(grd,temp1,strp)
          call mpi_gatherv(strp,grd%ijn(mype+1),mpi_real4, &
               tempa,grd%ijn,grd%displs_g,mpi_real4,0,mpi_comm_world,ierror)
          if(mype == 0) then
             do i=1,grd%iglobal
                temp(grd%ltosj(i),grd%ltosi(i))=tempa(i)
             end do
             if(e2a) then
                call unfill_nmm_grid(temp,nlon,nlat,f3d(:,:,k),igtype)
             else
                f3d(:,:,k)=temp
             end if
          end if
       end do
       if(mype == 0) then
          ncstatus = nf90_put_var(ncfileid,ncvarid,f3d)
       end if
    end do

    if(mype == 0)then
       deallocate(f2d,f3d,temp)
       ncstatus = nf90_close(ncfileid)
    end if
    deallocate(temp1,tempa,strp)

  end subroutine wrf_nmm_write

  subroutine fill_nmm_grid(grd,gin,nx,ny,gout,igtype,iorder)
  !$$$  subprogram documentation block
  !                .      .    .                                       .
  ! subprogram:    fill_nmm_grid2         fill holes in (wrf) nmm e-grid
  !   prgmmr: parrish          org: np22                date: 2004-06-22
  !
  ! abstract: creates an unstaggered A grid from the staggered E grid used 
  !           by the wrf nmm.  This is done by interpolation to fill the 
  !           holes in the E grid.  This is necessary because the gsi is 
  !           not yet able to work with anything other than unstaggered 
  !           grids.  This solution minimizes additional interpolation error
  !           but doubles the number of grid points.  This routine will be
  !           eliminated when the gsi has the capability to work directly 
  !           with staggered grids.
  !
  ! program history log:
  !   2004-06-22  parrish, document
  !   2013-10-25  todling - reposition ltosi and others to commvars
  !   2016-02-17  tong - modified the use grid info from grd
  !
  !   input argument list:
  !     gin      - input staggered E grid field over entire horizontal domain
  !     nx,ny    - input grid dimensions
  !     igtype   - =1, then (1,1) on staggered grid is at corner of grid 
  !                (mass point for nmm)
  !              - =2, then (1,1) is staggered (wind point for nmm, 
  !                see illustration below)
  !
  !                   igtype=1:
  !
  !
  !       ^   3             x     x     x     x
  !       |
  !       y   2                x     x     x     x
  !
  !           1             x     x     x     x
  !
  !                         1     2     3
  !
  !                           x -->
  !
  !                   igtype=2:
  !
  !
  !
  !       ^   3             x     x     x     x
  !       |
  !       y   2          x     x     x     x
  !
  !           1             x     x     x     x
  !
  !                         1     2     3
  !
  !                           x -->
  !
  !   output argument list
  !     gout     - output filled grid  (reorganized for distibution to local
  !     domains)
  !
  ! attributes:
  !   language: f90
  !   machine:  ibm RS/6000 SP
  !
  !$$$
    use kinds, only: r_single,r_kind,i_kind
    use general_sub2grid_mod, only: sub2grid_info
    use constants, only: half,quarter,zero

    implicit none

    type(sub2grid_info) ,intent(in   ) :: grd
    integer(i_kind),intent(in   ) :: nx,ny,igtype,iorder
    real(r_single) ,intent(in   ) :: gin(nx,ny)
    real(r_single) ,intent(  out) :: gout(grd%itotsub)

    real(r_single) b(2*nx-1,ny)
    integer(i_kind) i,im,ip,j,jm,jp
    real(r_single) fill,test
  
    fill=0.95_r_kind*huge(fill) ; test=0.95_r_kind*fill
    do j=1,ny
       do i=1,2*nx-1
          b(i,j)=fill
       end do
    end do
  
  ! First transfer all staggered points to appropriate
  ! points on filled output grid
    if(igtype==1) then
       do j=1,ny,2
          do i=1,nx
             b(2*i-1,j)=gin(i,j)
          end do
       end do
       do j=2,ny,2
          do i=1,nx-1
             b(2*i,j)=gin(i,j)
          end do
       end do
    else
       do j=1,ny,2
          do i=1,nx-1
             b(2*i,j)=gin(i,j)
          end do
       end do
       do j=2,ny,2
          do i=1,nx
             b(2*i-1,j)=gin(i,j)
          end do
       end do
    end if
    
  
  !  Now fill in holes
  
  ! Top and bottom rows:
    do j=1,ny,ny-1
       do i=1,2*nx-1
          if(b(i,j)>test) then
             ip=i+1 ; if(ip>2*nx-1) ip=i-1
             im=i-1 ; if(im<1) im=i+1
             b(i,j)=half*(b(im,j)+b(ip,j))
          end if
       end do
    end do
  
  
  ! Left and right rows:
    do j=1,ny
       jp=j+1 ; if(jp>ny)   jp=j-1
       jm=j-1 ; if(jm<1) jm=j+1
       do i=1,2*nx-1,2*nx-2
          if(b(i,j)>test) b(i,j)=half*(b(i,jm)+b(i,jp))
       end do
    end do
  
  ! Interior points
    do j=1,ny
       jp=j+1 ; if(jp>ny) jp=j-1
       jm=j-1 ; if(jm<1) jm=j+1
       do i=1,2*nx-1
          if(b(i,j)>test) then
             ip=i+1 ; if(ip>2*nx-1) ip=i-1
             im=i-1 ; if(im<1)      im=i+1
             b(i,j)=quarter*(b(ip,j)+b(im,j)+b(i,jp)+b(i,jm))
          end if
       end do
    end do

! Reorganize for eventual distribution to local domains
    do i=1,grd%itotsub
       gout(i)=zero
    end do
    if(iorder==1)then
       do i=1,grd%itotsub
          gout(i)=b(grd%ltosj_s(i),grd%ltosi_s(i))
       end do
    else
       do i=1,grd%iglobal
          gout(i)=b(grd%ltosj(i),grd%ltosi(i))
       end do
    endif
  
  end subroutine fill_nmm_grid

  subroutine unfill_nmm_grid(gout,nx,ny,gin,igtype)

  !$$$  subprogram documentation block
  !                .      .    .                                       .
  ! subprogram:    unfill_nmm_grid          opposite of fill_nmm_grid2
  !   prgmmr: parrish          org: np22                date: 2004-06-22
  !
  ! abstract: This is almost the reverse of subroutine fill_nmm_grid2.
  !           The input field is an analysis increment on an unstaggered
  !           A grid.  This routine extracts the points which coincide
  !           with the output E grid and adds them to the preexisting
  !           contents of the E grid.  Before this is done, the input
  !           grid must be reordered from the special ordering required
  !           for gathering up a full horizontal grid from the horizontal
  !           subdomains scattered across the processors.
  !
  !           See fill_nmm_grid2.f90 for additional comments.
  !
  ! program history log:
  !   2004-06-22  parrish, document
  !   2013-10-25  todling - reposition ltosi and others to commvars
  !   2016-02-17  tong - modified to directly take gout(2*nx-1,ny)
  !
  !   input argument list:
  !     gout     - input filled grid
  !     gin      - preexisting input values to be added to on staggered E grid
  !     nx,ny    - input grid dimensions
  !     igtype   - =1, then (1,1) on staggered grid is at corner of grid
  !                 (mass point for nmm)
  !              - =2, then (1,1) is staggered (wind point for nmm, see
  !                 illustration below)

  !   output argument list
  !     gin      - output result on staggered E grid
  !
  ! attributes:
  !   language: f90
  !   machine:  ibm RS/6000 SP
  !
  !$$$
    use kinds, only: r_single,i_kind
    implicit none

    integer(i_kind), intent(in   ) :: nx,ny,igtype
    real(r_single) , intent(in   ) :: gout(2*nx-1,ny)
    real(r_single) , intent(inout) :: gin(nx,ny)

    integer(i_kind) i,j

    if(igtype==1) then
       do j=1,ny,2
          do i=1,nx
             gin(i,j)=gout(2*i-1,j)
          end do
       end do
       do j=2,ny,2
          do i=1,nx-1
             gin(i,j)=gout(2*i,j)
          end do
       end do
    else
       do j=1,ny,2
          do i=1,nx-1
             gin(i,j)=gout(2*i,j)
          end do
       end do
       do j=2,ny,2
          do i=1,nx
             gin(i,j)=gout(2*i-1,j)
          end do
       end do
    end if

  end subroutine unfill_nmm_grid

  subroutine fill_nmm_grid2a3(gin,nx,ny,gout)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    fill_nmm_grid2a3
!   prgmmr:
!
! abstract:
!
! program history log:
!   2009-08-04  lueken - added subprogram doc block
!
!   input argument list:
!    nx,ny
!    gin
!
!   output argument list:
!    gout
!
! attributes:
!   language: f90
!   machine:
!
!$$$ end documentation block

  implicit none

  integer(i_kind),intent(in   ) :: nx,ny
  real(r_kind)   ,intent(in   ) :: gin(nx,ny)
  real(r_kind)   ,intent(  out) :: gout(2*nx-1,ny)

  integer(i_kind) i,j
  integer(i_kind) i1a(2*nx-1),i2a(2*nx-1)
  integer(i_kind) i3a(2*nx-1),i4a(2*nx-1)
  real(r_kind) r1a(2*nx-1),r2a(2*nx-1)
  real(r_kind) r3a(2*nx-1),r4a(2*nx-1)
  real(r_kind) x,x1,x2,x3,x4

!  first transfer all staggered points to appropriate
!   points on filled output grid

  do j=1,ny,2
     do i=1,nx
        gout(2*i-1,j)=gin(i,j)
     end do
  end do
  do j=2,ny,2
     do i=1,nx-1
        gout(2*i,j)=gin(i,j)
     end do
  end do

!   compute all interpolation constants for even x points on odd y rows

  i=2
  i1a(i)=i-1 ; i2a(i)=i+1 ; i3a(i)=i+3 ; i4a(i)=i+5
  x=i        ; x1=i1a(i)  ; x2=i2a(i)  ; x3=i3a(i)      ; x4=i4a(i)
  r1a(i)=       (x-x2)*(x-x3)*(x-x4)/(        (x1-x2)*(x1-x3)*(x1-x4))
  r2a(i)=(x-x1)       *(x-x3)*(x-x4)/((x2-x1)        *(x2-x3)*(x2-x4))
  r3a(i)=(x-x1)*(x-x2)       *(x-x4)/((x3-x1)*(x3-x2)        *(x3-x4))
  r4a(i)=(x-x1)*(x-x2)*(x-x3)       /((x4-x1)*(x4-x2)*(x4-x3)        )

  do i=4,2*nx-4,2
     i1a(i)=i-3 ; i2a(i)=i-1 ; i3a(i)=i+1 ; i4a(i)=i+3
     x=i        ; x1=i1a(i)  ; x2=i2a(i)  ; x3=i3a(i)      ; x4=i4a(i)
     r1a(i)=       (x-x2)*(x-x3)*(x-x4)/(        (x1-x2)*(x1-x3)*(x1-x4))
     r2a(i)=(x-x1)       *(x-x3)*(x-x4)/((x2-x1)        *(x2-x3)*(x2-x4))
     r3a(i)=(x-x1)*(x-x2)       *(x-x4)/((x3-x1)*(x3-x2)        *(x3-x4))
     r4a(i)=(x-x1)*(x-x2)*(x-x3)       /((x4-x1)*(x4-x2)*(x4-x3)        )
  end do

  i=2*nx-2
  i1a(i)=i-5 ; i2a(i)=i-3 ; i3a(i)=i-1 ; i4a(i)=i+1
  x=i        ; x1=i1a(i)  ; x2=i2a(i)  ; x3=i3a(i)     ; x4=i4a(i)
  r1a(i)=       (x-x2)*(x-x3)*(x-x4)/(        (x1-x2)*(x1-x3)*(x1-x4))
  r2a(i)=(x-x1)       *(x-x3)*(x-x4)/((x2-x1)        *(x2-x3)*(x2-x4))
  r3a(i)=(x-x1)*(x-x2)       *(x-x4)/((x3-x1)*(x3-x2)        *(x3-x4))
  r4a(i)=(x-x1)*(x-x2)*(x-x3)       /((x4-x1)*(x4-x2)*(x4-x3)        )

!   now get all interpolation constants for odd x points on even y rows

  i=1
  i1a(i)=i+1 ; i2a(i)=i+3 ; i3a(i)=i+5 ; i4a(i)=i+7
  x=i        ; x1=i1a(i)  ; x2=i2a(i)  ; x3=i3a(i)     ; x4=i4a(i)
  r1a(i)=       (x-x2)*(x-x3)*(x-x4)/(        (x1-x2)*(x1-x3)*(x1-x4))
  r2a(i)=(x-x1)       *(x-x3)*(x-x4)/((x2-x1)        *(x2-x3)*(x2-x4))
  r3a(i)=(x-x1)*(x-x2)       *(x-x4)/((x3-x1)*(x3-x2)        *(x3-x4))
  r4a(i)=(x-x1)*(x-x2)*(x-x3)       /((x4-x1)*(x4-x2)*(x4-x3)        )

  i=3
  i1a(i)=i-1 ; i2a(i)=i+1 ; i3a(i)=i+3 ; i4a(i)=i+5
  x=i        ; x1=i1a(i)  ; x2=i2a(i)  ; x3=i3a(i)         ; x4=i4a(i)
  r1a(i)=       (x-x2)*(x-x3)*(x-x4)/(        (x1-x2)*(x1-x3)*(x1-x4))
  r2a(i)=(x-x1)       *(x-x3)*(x-x4)/((x2-x1)        *(x2-x3)*(x2-x4))
  r3a(i)=(x-x1)*(x-x2)       *(x-x4)/((x3-x1)*(x3-x2)        *(x3-x4))
  r4a(i)=(x-x1)*(x-x2)*(x-x3)       /((x4-x1)*(x4-x2)*(x4-x3)        )

  do i=5,2*nx-5,2
     i1a(i)=i-3 ; i2a(i)=i-1 ; i3a(i)=i+1 ; i4a(i)=i+3
     x=i        ; x1=i1a(i)  ; x2=i2a(i)  ; x3=i3a(i)         ; x4=i4a(i)
     r1a(i)=       (x-x2)*(x-x3)*(x-x4)/(        (x1-x2)*(x1-x3)*(x1-x4))
     r2a(i)=(x-x1)       *(x-x3)*(x-x4)/((x2-x1)        *(x2-x3)*(x2-x4))
     r3a(i)=(x-x1)*(x-x2)       *(x-x4)/((x3-x1)*(x3-x2)        *(x3-x4))
     r4a(i)=(x-x1)*(x-x2)*(x-x3)       /((x4-x1)*(x4-x2)*(x4-x3)        )
  end do

  i=2*nx-3
  i1a(i)=i-5 ; i2a(i)=i-3 ; i3a(i)=i-1 ; i4a(i)=i+1
  x=i        ; x1=i1a(i)  ; x2=i2a(i)  ; x3=i3a(i)     ; x4=i4a(i)
  r1a(i)=       (x-x2)*(x-x3)*(x-x4)/(        (x1-x2)*(x1-x3)*(x1-x4))
  r2a(i)=(x-x1)       *(x-x3)*(x-x4)/((x2-x1)        *(x2-x3)*(x2-x4))
  r3a(i)=(x-x1)*(x-x2)       *(x-x4)/((x3-x1)*(x3-x2)        *(x3-x4))
  r4a(i)=(x-x1)*(x-x2)*(x-x3)       /((x4-x1)*(x4-x2)*(x4-x3)        )

  i=2*nx-1
  i1a(i)=i-7 ; i2a(i)=i-5 ; i3a(i)=i-3 ; i4a(i)=i-1
  x=i        ; x1=i1a(i)  ; x2=i2a(i)  ; x3=i3a(i) ; x4=i4a(i)
  r1a(i)=       (x-x2)*(x-x3)*(x-x4)/(        (x1-x2)*(x1-x3)*(x1-x4))
  r2a(i)=(x-x1)       *(x-x3)*(x-x4)/((x2-x1)        *(x2-x3)*(x2-x4))
  r3a(i)=(x-x1)*(x-x2)       *(x-x4)/((x3-x1)*(x3-x2)        *(x3-x4))
  r4a(i)=(x-x1)*(x-x2)*(x-x3)       /((x4-x1)*(x4-x2)*(x4-x3)        )

  do j=1,ny,2
     do i=2,2*nx-2,2
        gout(i,j)=r1a(i)*gout(i1a(i),j)+r2a(i)*gout(i2a(i),j)+ &
                  r3a(i)*gout(i3a(i),j)+r4a(i)*gout(i4a(i),j)
     end do
  end do
  do j=2,ny,2
     do i=1,2*nx-1,2
        gout(i,j)=r1a(i)*gout(i1a(i),j)+r2a(i)*gout(i2a(i),j)+ &
                  r3a(i)*gout(i3a(i),j)+r4a(i)*gout(i4a(i),j)
     end do
  end do

  end subroutine fill_nmm_grid2a3

  subroutine strip_grd_single(grd,field_in,field_out)

! !USES:

    use kinds, only: i_kind,r_kind
    use general_sub2grid_mod, only: sub2grid_info
    implicit none

! !INPUT PARAMETERS:

    type(sub2grid_info)                  ,intent(in   ) :: grd
    real(r_single),dimension(grd%lat2,grd%lon2), intent(in   ) :: field_in !
    ! full domain array containing buffer points

! !OUTPUT PARAMETERS:

    real(r_single),dimension(grd%lat1,grd%lon1), intent(  out) :: field_out !
! subdomain array with buffer points stripped off

! !DESCRIPTION: strip off buffer points froms subdomains for mpi comm
!               purposes
!
! !REVISION HISTORY:
!
!   2004-01-25  kleist
!   2004-05-14  kleist, documentation
!   2004-07-15  todling, protex-compliant prologue
!
! !REMARKS:
!
!   language: f90
!   machine:  ibm rs/6000 sp; sgi origin 2000; compaq/hp
!
! !AUTHOR:
!    kleist           org: np20                date: 2004-01-25
!
!EOP
!-------------------------------------------------------------------------

    integer(i_kind) i,j,jp1

    do j=1,grd%lon1
       jp1 = j+1
       do i=1,grd%lat1
          field_out(i,j)=field_in(i+1,jp1)
       end do
    end do

    return
  end subroutine strip_grd_single

  subroutine ll2rpolar(rlat,rlon,n,x_rpolar,y_rpolar,rlat0,rlon0,rotate3)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    ll2rpolar  convert earth lat-lon to rotated polar stereo
!   prgmmr:
!
! abstract:  Convert earth latitude and longitude to polar stereographic
! coordinates where
!               the reference pole is centered at earth coordinates rlat0,rlon0.
!               The polar stereographic positive x axis is oriented
!               counterclockwise to
!               earth direction south at rlat0,rlon0 by angle rotate3.  The
!               transformed
!               lat-lon coordinate which is the basis of the polar stereographic
!               coordinate consists
!               of a sequence of 3 rotations in 3-d x-y-z space with origin at
!               center of earth, where
!               the x-y plane intersects the equator with the positive x axis
!               intersecting the 0 meridian.
!               and the positive y axis the 90E meridian.  The positive z axis
!               intersects the north pole.
!               1st rotation: counterclockwise in x-y plane by amount rlon0.
!               2nd rotation: counterclockwise in z-x plane by amount pi/2 -
!               rlat0.
!               3rd rotation: counterclockwise in x-y plane by amount rotate3.
!
! program history log:
!   2010-09-09  parrish - initial documentation
!
!   input argument list:
!    rlat,rlon:    input earth latitude and longitude coordinates in radians
!    n:            number of points to compute transform coordinates
!    rlat0,rlon0:  earth coordinates of north pole of new coordinate system
!    rotate3:      angle counterclockwise from earth direction south at
!    rlat0,rlon0 to positive x axis
!                     of output coordinates
!
!   output argument list:
!    x_rpolar,y_rpolar:  x and y polar stereographic coordinates of new
!    coordinate with origin at
!                           rlat0,rlon0.
!
! attributes:
!   language: f90
!   machine:
!
!$$$ end documentation block

   use kinds, only: r_kind,i_kind
   use constants, only: zero,one,two

   implicit none

   integer(i_kind),intent(in)::n
   real(r_kind),intent(in)::rlat(n),rlon(n)
   real(r_kind),intent(in)::rlat0,rlon0,rotate3
   real(r_kind),intent(out)::x_rpolar(n),y_rpolar(n)

!  Declare local variables
   integer(i_kind) i
   real(r_kind) clat0,slat0,clon0,slon0
   real(r_kind) clat(n),slat(n),clon(n),slon(n)
   real(r_kind) x,y,z,xt,yt,zt,x2,y2,z2,rlat2,rlon2,epstest,r_polar

   epstest=epsilon(epstest)
   clat0=cos(rlat0) ; slat0=sin(rlat0)
   clon0=cos(rlon0) ; slon0=sin(rlon0)
   clat =cos(rlat ) ; slat =sin(rlat )
   clon =cos(rlon ) ; slon =sin(rlon )

   do i=1,n
      x=clat(i)*clon(i) ; y=clat(i)*slon(i) ; z=slat(i)
      xt=x*clon0+y*slon0 ; yt=-x*slon0+y*clon0 ; zt=z
      z2=zt*slat0+xt*clat0 ; x2 = -zt*clat0 + xt*slat0 ; y2 = yt
      z2=min(one,max(-one,z2))
      rlat2=asin(z2)
      if(sqrt(y2**2+x2**2) < epstest) then
         rlon2=zero
      else
         rlon2=atan2(y2,x2)
      end if
      r_polar=cos(rlat2)/(one+sin(rlat2))
      x_rpolar(i)=r_polar*cos(rlon2-rotate3)
      y_rpolar(i)=r_polar*sin(rlon2-rotate3)

   end do

  end subroutine ll2rpolar

  subroutine rpolar2ll(x_rpolar,y_rpolar,n,rlat,rlon,rlat0,rlon0,rotate3)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    rpolar2ll  inverse of ll2rpolar
!   prgmmr:
!
! abstract:  Inverse transformation of subroutine ll2rpolar.
!
! program history log:
!   2010-09-09  parrish - initial documentation
!
!   input argument list:
!    x_rpolar,y_rpolar:  x and y polar stereographic coordinates of new
!    coordinate with origin at
!                           rlat0,rlon0.
!    n:            number of points to compute transform coordinates
!    rlat0,rlon0:  earth coordinates of north pole of new coordinate system
!    rotate3:      angle counterclockwise from earth direction south at
!    rlat0,rlon0 to positive x axis
!                     of output coordinates
!
!   output argument list:
!    rlat,rlon:    input earth latitude and longitude coordinates in radians
!
! attributes:
!   language: f90
!   machine:
!
!$$$ end documentation block

   use kinds, only: r_kind,i_kind
   use constants, only: zero,one,two,quarter,pi

   implicit none

   integer(i_kind),intent(in)::n
   real(r_kind),intent(in)::x_rpolar(n),y_rpolar(n)
   real(r_kind),intent(in)::rlat0,rlon0,rotate3
   real(r_kind),intent(out)::rlat(n),rlon(n)

!  Declare local variables
   integer(i_kind) i
   real(r_kind) clat0,slat0,clon0,slon0
   real(r_kind) x,y,z,xt,yt,zt,x2,y2,z2,rlat2,rlon2,epstest,r_polar,pi_quarter
   real(r_kind) slat2,clat2,slon2,clon2

   epstest=epsilon(epstest)
   pi_quarter=quarter*pi
   clat0=cos(rlat0) ; slat0=sin(rlat0)
   clon0=cos(rlon0) ; slon0=sin(rlon0)

   do i=1,n
      r_polar=sqrt(x_rpolar(i)**2+y_rpolar(i)**2)
      rlat2=two*(pi_quarter-atan(r_polar))
      slat2=sin(rlat2) ; clat2=cos(rlat2)
      if(r_polar < epstest) then
         rlon2=rotate3
      else
         rlon2=atan2(y_rpolar(i),x_rpolar(i))+rotate3
      end if
      slon2=sin(rlon2) ; clon2=cos(rlon2)
      x2=clat2*clon2 ; y2=clat2*slon2 ; z2=slat2
      zt=slat0*z2-clat0*x2 ; xt=clat0*z2+slat0*x2 ; yt=y2
      x=xt*clon0-yt*slon0 ; y=xt*slon0+yt*clon0 ; z=zt
      z=min(one,max(-one,z))
      rlat(i)=asin(z)
      if(sqrt(x**2+y**2) < epstest) then
         rlon(i)=zero
      else
         rlon(i)=atan2(y,x)
      end if
   end do

  end subroutine rpolar2ll

subroutine strip_grd_double(grd,field_in,field_out)

! !USES:

    use kinds, only: i_kind,r_kind
    use general_sub2grid_mod, only: sub2grid_info
    implicit none

! !INPUT PARAMETERS:

    type(sub2grid_info)                  ,intent(in   ) :: grd
    real(r_kind),dimension(grd%lat2,grd%lon2), intent(in   ) :: field_in    ! full subdomain
                                                                       !    array containing
                                                                       !    buffer points
! !OUTPUT PARAMETERS:

    real(r_kind),dimension(grd%lat1,grd%lon1), intent(  out) :: field_out  ! subdomain array
                                                                      !   with buffer points
                                                                      !   stripped off

! !DESCRIPTION: strip off buffer points froms subdomains for mpi comm
!               purposes
!
! !REVISION HISTORY:
!
!   2004-01-25  kleist
!   2004-05-14  kleist, documentation
!   2004-07-15  todling, protex-compliant prologue
!
! !REMARKS:
!
!   language: f90
!   machine:  ibm rs/6000 sp; sgi origin 2000; compaq/hp
!
! !AUTHOR:
!    kleist           org: np20                date: 2004-01-25
!
!EOP
!-------------------------------------------------------------------------

    integer(i_kind) i,j,jp1

    do j=1,grd%lon1
       jp1 = j+1
       do i=1,grd%lat1
          field_out(i,j)=field_in(i+1,jp1)
       end do
    end do

    return
end subroutine strip_grd_double

subroutine sub2grid_3a(grd,sub,grid,gridpe,mype)

!     straightforward, but inefficient code to convert a single variable on subdomains to complete
!      slab on one processor.
!  2013-10-24 todling - revisit strip interface

  use kinds, only: r_kind,i_kind
  use constants, only: zero
  use mpimod, only: mpi_comm_world,ierror,mpi_rtype
  use general_sub2grid_mod, only: sub2grid_info
  implicit none

  type(sub2grid_info)                  ,intent(in   ) :: grd
  integer(i_kind), intent(in)::gridpe,mype
  real(r_kind),dimension(grd%lat2,grd%lon2),intent(in):: sub
  real(r_kind),dimension(grd%nlat,grd%nlon),intent(out)::grid

  real(r_kind),dimension(grd%lat1*grd%lon1):: zsm
  real(r_kind),dimension(grd%itotsub):: work1
  integer(i_kind) mm1,i,j,k

  mm1=mype+1

  do j=1,grd%lon1*grd%lat1
    zsm(j)=zero
  end do
  call strip_grd_double(grd,sub,zsm)
  call mpi_gatherv(zsm,grd%ijn(mm1),mpi_rtype, &
                 work1,grd%ijn,grd%displs_g,mpi_rtype, &
                 gridpe,mpi_comm_world,ierror)
  if(mype==gridpe) then
    do k=1,grd%iglobal
      i=grd%ltosi(k) ; j=grd%ltosj(k)
      grid(i,j)=work1(k)
    end do
  end if

end subroutine sub2grid_3a

subroutine grads3d(grd,u,v,tsen,q,pd,nvert,mype,fname)

  use kinds, only: r_kind,i_kind,r_single
  use general_sub2grid_mod, only: sub2grid_info
  implicit none

  type(sub2grid_info)                  ,intent(in   ) :: grd
  integer(i_kind) nvert
  integer(i_kind), intent(in)::mype
  character(*) fname
  real(r_kind),dimension(grd%lat2,grd%lon2,nvert):: u,v,tsen,q
  real(r_kind),dimension(grd%lat2,grd%lon2):: pd

  real(r_kind),dimension(grd%nlat,grd%nlon)::work
  real(r_single) outfield(grd%nlon,grd%nlat)

  character(50) dsname,title,filename
! data dsname/'test.dat'/
  data title/'inmi'/
  character(112) datdes(50000)
  character(1) blank
  data blank/' '/
  data undef/-9.99e33_r_single/

  integer(i_kind) i,k,next,ioutdes,ioutdat
  integer(i_kind) last,j,koutmax
  real(r_single) undef
  real(r_single) startp,pinc

  if(mype==0) then
    startp=1._r_single
    pinc=1._r_single
    ioutdes=98750
    ioutdat=98751
    write(filename,'(a,".des")')trim(fname)
    write(dsname,'(a,".dat")')trim(fname)
    open(unit=ioutdes,file=trim(filename),form='formatted')
    open(unit=ioutdat,file=trim(dsname),form='unformatted')
    rewind ioutdes
    rewind ioutdat
    do i=1,50000
      write(datdes(i),'(112a1)')(blank,k=1,112)
    end do
    write(datdes(1),'("DSET ^",a50)')dsname
    write(datdes(2),'("options big_endian sequential")')
    write(datdes(3),'("TITLE ",a50)')title
    write(datdes(4),'("UNDEF ",e11.2)')undef
    next=5
    write(datdes(next),'("XDEF ",i5," LINEAR ",f7.2,f7.2)')grd%nlon,startp,pinc
    next=next+1
    write(datdes(next),'("YDEF ",i5," LINEAR ",f7.2,f7.2)')grd%nlat,startp,pinc
    next=next+1
    write(datdes(next),'("ZDEF ",i5," LINEAR ",f7.2,f7.2)')nvert,startp,pinc
    next=next+1
    koutmax=1
    write(datdes(next),'("TDEF ",i5," LINEAR 00Z01Jan2000 12hr")')koutmax
    next=next+1
    write(datdes(next),'("VARS 5")')
    next=next+1
    write(datdes(next),'("u   ",i5," 99 u   ")')nvert
    next=next+1
    write(datdes(next),'("v   ",i5," 99 v   ")')nvert
    next=next+1
    write(datdes(next),'("t   ",i5," 99 t   ")')nvert
    next=next+1
    write(datdes(next),'("q   ",i5," 99 q   ")')nvert
    next=next+1
    write(datdes(next),'("pd  ",i5," 99 pd  ")')0
    next=next+1
    write(datdes(next),'("ENDVARS")')
    last=next
    write(ioutdes,'(a112)')(datdes(i),i=1,last)
  endif
  do k=1,nvert
    call sub2grid_3a(grd,u(1,1,k),work,0,mype)
    if(mype==0) then
      do j=1,grd%nlon ; do i=1,grd%nlat
          outfield(j,i)=work(i,j)
      end do ; end do
      write(ioutdat)outfield
    end if
  end do

  do k=1,nvert
    call sub2grid_3a(grd,v(1,1,k),work,0,mype)
    if(mype==0) then
      do j=1,grd%nlon ; do i=1,grd%nlat
          outfield(j,i)=work(i,j)
      end do ; end do
      write(ioutdat)outfield
    end if
  end do

  do k=1,nvert
    call sub2grid_3a(grd,tsen(1,1,k),work,0,mype)
    if(mype==0) then
      do j=1,grd%nlon ; do i=1,grd%nlat
          outfield(j,i)=work(i,j)
      end do ; end do
      write(ioutdat)outfield
    end if
  end do

  do k=1,nvert
    call sub2grid_3a(grd,q(1,1,k),work,0,mype)
    if(mype==0) then
      do j=1,grd%nlon ; do i=1,grd%nlat
          outfield(j,i)=work(i,j)
      end do ; end do
      write(ioutdat)outfield
    end if
  end do

  call sub2grid_3a(grd,pd(1,1),work,0,mype)
  if(mype==0) then
    do j=1,grd%nlon ; do i=1,grd%nlat
        outfield(j,i)=work(i,j)
    end do ; end do
    write(ioutdat)outfield
  end if

  if(mype==0) then
    close(ioutdes)
    close(ioutdat)
  end if
end subroutine grads3d

end module wrf_nmm_io

