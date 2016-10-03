program compute_geomdist
  ! Compute the distance between two Columbus geometries.
  use ioutil
  character(255) :: gstart_file, gfinal_file
  character(255) :: cmdstr
  integer :: natoms
  double precision, dimension(:,:), allocatable :: gstart, gfinal
  double precision, dimension(:,:), allocatable :: gdist
  integer :: ios

  call get_command_argument(1, value=cmdstr, status=ios)
  if (ios .eq. 0) then
          read(cmdstr,*) natoms
  else
          stop "*** Error reading argument = 1 ***"
  end if
  call get_command_argument(2, value=gstart_file, status=ios)
  if (ios .ne. 0) stop "*** Error reading argument =2 ***"
  call get_command_argument(3, value=gfinal_file, status=ios)
  if (ios .ne. 0) stop "*** Error reading argument =3 ***"

  allocate(gstart(3,natoms))
  allocate(gfinal(3,natoms))
  allocate(gdist(3,natoms))
  call read_colgeom(gstart_file, gstart, natoms)
  call read_colgeom(gfinal_file, gfinal, natoms)

  do ios=1, natoms
          gdist(1:3,ios) = gfinal(1:3,ios)-gstart(1:3,ios)
  end do
  call print_colgeom(gdist, natoms)
  deallocate(gstart, gfinal, gdist)
end program compute_geomdist
