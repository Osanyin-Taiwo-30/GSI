module namelist

  !-------------------------------------------------------------------
  !$$$   module documentation block
  !               
  ! module:    define namelist variables and read namelist
  ! PRGMMR: Mingjing Tong                          DATE: 2016-02-15
  !
  ! Variable Definitions:
  ! intrp_from_file: name of source file to be interpolated from
  ! intrp_to_file: name of target file to be interpolated to
  ! e2a: convert E-grid to A-grid before interpolation, if true
  ! n3d: number of 3D variables to be interpolated
  ! n2d: number of 2D variables to be interpolated
  ! var3d: name of the 3D variables
  ! var2d: name of the 2D variables
  ! iblend: blending option 1. around e-grid boundary
  !                         2. circle around storm
  ! nord_blend: order of continuity of blend function 
  !                   (1=continuous 1st derivative, etc)
  ! nmix: width of blend zone on edge of e grid in e grid units. 
  !       (when iblend=1)
  ! radmax: maximum radius from TC center, with in which model fields 
  !         will be interpolated from source domain to target domain
  !         (in degree, iblend=2)
  ! bwdth: blending ozone width factor (the fraction of the radius 
  !        that is used to define the interpolation zone) 
  !$$$ end documentation block
  !-------------------------------------------------------------------

  use kinds
  use mpimod, only: mpi_comm_world,ierror,mype

  implicit none

  ! set default as private
  private

  ! set subroutines as public
  public :: namelistparams
  public :: intrp_from_file,intrp_to_file
  public :: e2a,n3d,n2d,var3d,var2d
  public :: iblend,nord_blend,nmix,radmax,bwdth
  
  ! Define global variables
  character(len=500)           :: intrp_from_file       = 'NOT USED'
  character(len=500)           :: intrp_to_file         = 'NOT USED'
  logical                      :: e2a                   = .true.
  integer(i_kind)              :: n3d                   = 4
  integer(i_kind)              :: n2d                   = 1
  character(8),dimension(10)   :: var3d                 = '        '
  character(8),dimension(5)    :: var2d                 = '        '
  integer(i_kind)              :: iblend                = 1 ! blend option
  integer(i_kind)              :: nord_blend            = 4
  integer(i_kind)              :: nmix                  = 10
  real(r_kind)                 :: radmax                = 10.0
  real(r_kind)                 :: bwdth                 = 0.33
  namelist /setup/ intrp_from_file,intrp_to_file,e2a, &
                   n3d,n2d,var3d,var2d,iblend,nord_blend, &
                   nmix,radmax,bwdth

contains

  subroutine namelistparams()

    !-------------------------------------------------------------------
    !$$$  subprogram documentation block
    ! 
    ! subprogram:   namelistparams 
    ! PRGMMR: Mingjing Tong
    ! abstract: read namelist parameters
    !
    !$$$ end documentation block 
    !-------------------------------------------------------------------


    ! Define variables computed within routine
    logical                                               :: is_it_there

    ! Define counting variables
    integer                                               :: i

    is_it_there = .false.
    inquire(file='hwrf_interpolate.nml',exist = is_it_there)

    if(is_it_there) then

       open(7,file = 'hwrf_interpolate.nml',delim='APOSTROPHE')
       read(7,NML=setup)
       close(7)

       do i=1,n3d
          if(trim(var3d(i))=='U' .and. i /= 1)then
             var3d(i)=var3d(1)
             var3d(1)='U' 
          end if
          if(trim(var3d(i))=='V' .and. i /= 2)then
             var3d(i)=var3d(2)
             var3d(2)='V'
          end if
       end do
    else
       if(mype == 0) write(6,500)
       call stop2(999)
    end if

    if(mype==0) then

       write(6,*) '&setup'
       write(6,*) 'intrp_from_file        = ', trim(intrp_from_file)
       write(6,*) 'intrp_to_file          = ', trim(intrp_to_file)              
       write(6,*) 'e2a                    = ', e2a
       write(6,*) 'n3d                    = ', n3d
       write(6,*) 'n2d                    = ', n2d
       write(6,*) 'var3d                  = ', (var3d(i),i=1,n3d)
       write(6,*) 'var2d                  = ', (var2d(i),i=1,n2d)
       write(6,*) 'iblend                 = ', iblend
       write(6,*) 'nord_blend             = ', nord_blend
       write(6,*) 'nmix                   = ', nmix
       write(6,*) 'radmax                 = ', radmax
       write(6,*) 'bwdth                  = ', bwdth
       write(6,*) '/'

    end if

    call mpi_barrier(mpi_comm_world,ierror)

500 format('NAMELISTPARAMS: hwrf_interpolate.nml not found in the ',       &
         & 'current working directory. ABORTING!!!!')

  end subroutine namelistparams

end module namelist
