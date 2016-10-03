program mabeval
  !-------------------------------------------------------------------
  ! Program to evaluate the molecular aharonov-bohm effect.
  !
  ! cthree-40
  ! Yarkony Group
  ! Department of Chemistry
  ! The Johns Hopkins University
  !
  ! -- Change Log --
  ! 2016-07-21 - Created.
  ! 2016-08-31 - Moved read_colgeom, get_unit to exteral module: ioutil.f90
  !-------------------------------------------------------------------
  use mabutil, only: evaluate_mab
  use ioutil,  only: print_colgeom, read_colgeom, get_unit
  implicit none
  integer :: natoms ! Number of atoms
  character(255) :: gmfile ! Geometry file
  character(255) :: cmdstr ! Command line string
  double precision, dimension(:,:), allocatable :: geom
  double precision :: x0, y0, z0 ! Origin of loop (circle)
  double precision :: rho ! radius of loop (circle)
  double precision :: gphase ! phase of g vector
  integer :: maxint ! integration iterations
  logical :: remtranrot ! remove translations and rotations from Zspace
  double precision :: startangle ! starting angle for loop (degrees)
  integer :: circdir ! circulation direction, >= 0: counterclockwise
  !                                           <  0: clockwise
  integer :: rotations ! number of rotations
  integer :: v1, v2 ! x, y vectors
  integer, dimension(:), allocatable :: v3 !z vectors
  logical :: energy_plot
  integer :: ios
  
  ! Get command line arguments
  call get_command_argument(1, value=cmdstr, status=ios)
  if (ios .eq. 0) then
          read(cmdstr,*) natoms
  else
          stop "*** Error reading argument = 1 ***"
  end if
  call get_command_argument(2, value=gmfile, status=ios)
  if (ios .ne. 0) stop "*** Error reading argument = 2 ***"
  print "(A,i4)", "Number of atoms = ", natoms
  print "(A,A)",  "Geometry file: ", trim(adjustl(gmfile))

  ! Read in geometry
  allocate(geom(3,natoms))
  call read_colgeom(gmfile, geom, natoms)
  call print_colgeom(geom, natoms)

  ! Read input
  allocate(v3(3*natoms - 6))
  call read_mabinput(x0, y0, z0, rho, gphase, v1, v2, v3, energy_plot, &
          (3*natoms - 6), maxint, remtranrot, startangle, circdir,     &
          rotations)
  
  ! Evaluate MAB
  call evaluate_mab(geom, x0, y0, z0, rho, natoms, gphase, v1, v2, v3, &
          maxint, remtranrot, startangle, circdir, rotations)
contains
  subroutine read_mabinput (x0, y0, z0, rho, gphase, vector1, vector2, vector3, &
          energy_plot, na3, maxint, remtranrot, startangle, circdir, rotations)
    implicit none
    ! ..input scalars..
    integer, intent(in) :: na3
    ! ..output scalars..
    integer, intent(out) :: vector1, vector2, maxint, circdir, rotations
    double precision, intent(out) :: startangle
    ! ..output arrays..
    integer, dimension(na3), intent(out) :: vector3
    double precision, intent(out) :: x0, y0, z0, rho, gphase
    ! ..output logicals..
    logical, intent(out) :: energy_plot, remtranrot
    ! ..local scalars..
    integer :: flun
    integer :: ios
    ! ..local arrays..
    character(255) :: flnm
    ! ..namelist ..
    namelist /general/ x0, y0, z0, rho, gphase, vector1, vector2, vector3, &
            maxint, remtranrot, startangle, circdir, rotations
    namelist /plotting/ energy_plot

    x0=0d0
    y0=0d0
    z0=0d0
    gphase=1d0
    rho=1d-2
    vector1 = 1
    vector2 = 2
    vector3 = 0
    maxint = 5000
    startangle = 0d0
    circdir = 100
    rotations = 1
    energy_plot = .false.
    remtranrot = .false.
    flun=get_unit()
    flnm="mab.input"
    open(file=trim(adjustl(flnm)),unit=flun,action="read",position="rewind",&
            status="old",iostat=ios)
    if (ios .eq. 0) then
            read(unit=flun,nml=general,iostat=ios)
            print *, "NAMELIST: &GENERAL READ: ios = ", ios
            read(unit=flun,nml=plotting,iostat=ios)
            print *, "NAMELIST: &PLOTTING READ: ios = ", ios
    end if
    close(flun)
    return
  end subroutine read_mabinput
end program mabeval
