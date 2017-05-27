program ghz_position
  use ioutil
  ! Compute (x, y, z) coordinates of input geometry.
  implicit none
  double precision, dimension(:,:), allocatable :: gm, xgm
  double precision, dimension(:), allocatable :: gv, hv
  double precision, dimension(:,:), allocatable :: zv
  integer :: natoms
  integer :: flunit
  character(255) :: gmfile, xgmfile
  character(255) :: cmdstr
  character(255) :: gvfl, hvfl, zvfl
  double precision :: x, y, z
  double precision, external :: ddot
  integer :: ios
  integer :: nzvec
  integer :: i
  
  ! Read geometry from geometry input file. Read G, H, and
  ! Z vectors in to respective arrays. Then take the dot product
  ! of G, H, and Z with the input geometry.

  call get_command_argument(1, value=cmdstr, status=ios)
  if (ios .eq. 0) then
          read(cmdstr,*) natoms
  else
          stop "*** Error reading argument = 1 ***"
  end if
  print "(A,i4)", "Number of atoms = ", natoms
  call get_command_argument(2, value=xgmfile, status=ios)
  if (ios .ne. 0) stop "*** Error reading argument = 2 ***"
  print "(A,A)",  "MEX Geometry file: ", trim(adjustl(xgmfile))
  call get_command_argument(3, value=gmfile, status=ios)
  if (ios .ne. 0) stop "*** Error reading argument = 3 ***"
  print "(A,A)",  "Geometry file: ", trim(adjustl(gmfile))
  call get_command_argument(4, value=cmdstr, status=ios)
  if (ios .eq. 0) then
          read(cmdstr,*) nzvec
  else
          stop "*** Error reading argument = 4 ***"
  end if
  
  allocate(gm(3,natoms))
  allocate(xgm(3,natoms))
  allocate(gv(3*natoms))
  allocate(hv(3*natoms))
  allocate(zv(3*natoms, nzvec))

  call read_colgeom(xgmfile, xgm, natoms)
  call read_colgeom(gmfile, gm, natoms)
  write(gvfl, "(A)") "gvector.dat"
  write(hvfl, "(A)") "hvector.dat"
  write(zvfl, "(A)") "zvector.dat"
  flunit=get_unit()
  call read_vector(gv, natoms, flunit, gvfl, 1)
  if (flunit .lt. 0) stop "*** COULD NOT OPEN G VECTOR FILE ***"
  flunit=get_unit()
  call read_vector(hv, natoms, flunit, hvfl, 1)
  if (flunit .lt. 0) stop "*** COULD NOT OPEN H VECTOR FILE ***"
  flunit=get_unit()
  call read_vector(zv, natoms, flunit, zvfl, nzvec)
  if (flunit .lt. 0) print *, "WARNING: No Z Vector file found!"
  print *, "Norm of Z Vector: ", sqrt(ddot(natoms*3, zv, 1, zv, 1))
  zv = zv / sqrt(ddot(natoms*3, zv, 1, zv, 1))
  gm= gm - xgm
  x = ddot(natoms*3, gv, 1, gm, 1)
  y = ddot(natoms*3, hv, 1, gm, 1)
  do i = 1, nzvec
          print "('Norm g.z_i = ',f12.8)", ddot(natoms*3,zv(1:natoms*3,i),1,gv,1)
          print "('Norm h.z_i = ',f12.8)", ddot(natoms*3,zv(1:natoms*3,i),1,hv,1)
          z = ddot(natoms*3, zv(1:natoms*3,i), 1, gm, 1)
          print "('(',f12.8,',',f12.8,',',f12.8,')')", x, y, z
  end do
  deallocate(gv, hv, zv, gm, xgm)
  
end program ghz_position
