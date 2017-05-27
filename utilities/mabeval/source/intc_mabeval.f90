program intc_mabeval
  use ioutil, only: read_colgeom
  use intc_mabutil, only: evaluate_intc_mabeval, convertgeom_bohr2ang
  implicit none
  integer :: natoms
  double precision, dimension(:,:), allocatable :: geom
  character(255) :: cmdstr, gmfile
  integer :: ios
  integer :: tmp
  integer :: coord1, coord2, maxint, circdir, rotations
  integer :: state1, state2
  ! epsilon = value of loop radius
  double precision :: epsilon
  
  call get_command_argument(1, value=cmdstr, status=ios)
  if (ios .eq. 0) then
          read(cmdstr, *) natoms
  else
          stop "Usage: intc_mabeval.x [atoms]* [geom file]"
  end if

  call get_command_argument(2, value=gmfile, status=ios)
  if (ios .ne. 0) then
          stop "Usage: intc_mabeval.x [atoms] [geom file]*"
  end if
  print "(A,A)",  "MEX Geometry file: ", trim(adjustl(gmfile))
  allocate(geom(3,natoms))
  call read_colgeom(gmfile, geom, natoms)

  call read_input(coord1, coord2, epsilon, maxint, circdir, rotations, &
          state1, state2)
  
  call evaluate_intc_mabeval(natoms, geom, coord1, coord2, epsilon, &
          maxint, circdir, rotations, state1, state2)

contains
  subroutine read_input (coord1, coord2, epsilon, maxint, circdir, rotations,&
          state1, state2)
    use ioutil, only: get_unit
    implicit none
    integer, intent(out) :: coord1, coord2, maxint, circdir, rotations
    integer, intent(out) :: state1, state2
    double precision, intent(out) :: epsilon
    character(255) :: flnm
    integer :: flun, ios
    namelist /general/ coord1, coord2, maxint, circdir, rotations, epsilon, &
            state1, state2

    coord1 = 1
    coord2 = 2
    maxint = 5000
    circdir = 100
    rotations = 1
    epsilon = 0.1d0
    state1=1
    state2=2
    flnm = "intcmab.input"
    flun = get_unit()
    open(file=trim(adjustl(flnm)),unit=flun,action="read",status="old",&
            iostat=ios)
    if (ios .eq. 0) then
            read(unit=flun,nml=general,iostat=ios)
            print *, "NAMELIST: &GENERAL READ: ios = ", ios
    end if
    close(flun)
    return
  end subroutine read_input
end program intc_mabeval
