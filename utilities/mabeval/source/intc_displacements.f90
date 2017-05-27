program intc_displacements
  use ioutil, only: read_colgeom2
  implicit none
  integer :: natoms
  character(255) :: cmdstr
  character(255) :: gmfile
  integer :: ios
  integer :: coord1, coord2, disp
  double precision :: dispsize
  double precision, dimension(:,:), allocatable :: geom
  double precision, dimension(:), allocatable :: masses, nums
  character(2), dimension(:), allocatable :: atoms
  
  call get_command_argument(1, value=cmdstr, status=ios)
  if (ios .eq. 0) then
          read(cmdstr, *) natoms
  else
          stop "Usage: intc_displacements.x [natoms]* [geom file] "//&
                  "[coord1] [coord2] [displacements in each dir]"
  end if
  call get_command_argument(2, value=gmfile, status=ios)
  if (ios .ne. 0) then
          stop  "Usage: intc_displacements.x [natoms] [geom file]* "//&
                  "[coord1] [coord2] [displacements in each dir]"
  end if
  call get_command_argument(3, value=cmdstr, status=ios)
  if (ios .eq. 0) then
          read(cmdstr, *) coord1
  else
          stop "Usage: intc_displacements.x [natoms] [geom file] "//&
                  "[coord1]* [coord2] [displacements in each dir]"
  end if
  call get_command_argument(4, value=cmdstr, status=ios)
  if (ios .eq. 0) then
          read(cmdstr, *) coord2
  else
          stop "Usage: intc_displacements.x [natoms] [geom file] "//&
                  "[coord1] [coord2]* [displacements in each dir]"
  end if
  call get_command_argument(5, value=cmdstr, status=ios)
  if (ios .eq. 0) then
          read(cmdstr, *) disp
  else
          stop "Usage: intc_displacements.x [natoms] [geom file] "//&
                  "[coord1] [coord2] [displacements in each dir]*"
  end if
  call get_command_argument(6, value=cmdstr, status=ios)
  if (ios .eq. 0) then
          read(cmdstr, *) dispsize
  else
          stop "Usage: intc_displacements.x [natoms] [geom file] "//&
                  "[coord1] [coord2] [displacements in each dir] "//&
                  "[disp size]*"
  end if

  print "(A,A)", "Geometry file: ", trim(adjustl(gmfile))
  allocate(geom(3,natoms))
  allocate(atoms(natoms))
  allocate(nums(natoms))
  allocate(masses(natoms))
  call read_colgeom2(gmfile, geom, natoms, masses, atoms, nums)
  print "(3f14.8)", geom
  call make_displacements(geom, natoms, coord1, coord2, disp, dispsize, atoms, &
          masses, nums)

contains
  subroutine make_displacements(gm, na, c1, c2, disp, dispsize, atoms, masses, &
          nums)
    use intc_mabutil
    use ioutil, only: get_unit, print_colgeom
    implicit none
    integer, intent(in) :: na, c1, c2
    integer, intent(inout) :: disp
    double precision, dimension(3,na), intent(in) :: gm
    double precision, dimension(na), intent(in) :: masses, nums
    character(2), dimension(na), intent(in) :: atoms
    double precision, intent(in) :: dispsize

    double precision, dimension(3,na) :: gm_0, gm_ang, gm_bohr
    double precision, dimension(594) :: qq_0, qq, qq_scr
    integer :: intcfl_unit
    integer, dimension(6,1782)  :: ibcode
    integer, dimension(20,594) :: ibcontr
    double precision, dimension(54,594) :: bmat
    double precision, dimension(600,594) :: b

    integer :: i, j, k
    
    ! Convert original geometry to angstrom for cart2int routines
    gm_0 = gm
    call convertgeom_bohr2ang(gm_0, na)
    intcfl_unit = get_unit()
    call open_intcfl_file(intcfl_unit)
    call read_intcfl(intcfl_unit, na, ibcode, ibcontr)
    call generate_bmat(natoms, gm_0, ibcode, ibcontr, bmat, qq_0)

    if (mod(disp,2) .ne. 0) disp = disp+1

    qq_0(c1) = qq_0(c1) - dispsize*(disp/2)
    qq_0(c2) = qq_0(c2) - dispsize*(disp/2)
    do i = 1, disp
            do j = 1, disp
                    qq=0d0
                    qq(c1) = dispsize*i
                    qq(c2) = dispsize*j
                    qq_scr = 0d0
                    call distort((3*na), (3*na-6), na, qq, qq_0, qq_scr, gm_0, &
                            ibcode, ibcontr, bmat, MATH_PI/2, gm_ang)
                    gm_bohr = gm_ang
                    call convertgeom_ang2bohr(gm_bohr, na)
                    print *, ""
                    print *, "New geometry: "
                    do k = 1, na
                            print "(1x,a2,f8.1,3f14.8,f14.8)",atoms(k),nums(k), &
                                    gm_bohr(1,k), gm_bohr(2,k), gm_bohr(3,k),   &
                                    masses(k)
                    end do
            end do
    end do
    return
  end subroutine make_displacements
end program intc_displacements
                    
    
  
  
