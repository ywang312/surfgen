program intc_mabeval
  use ioutil, only: read_colgeom
  use intc_mabutil, only: evaluate_intc_mabeval, convertgeom_bohr2ang
  implicit none
  integer :: natoms
  double precision, dimension(:,:), allocatable :: geom
  character(255) :: cmdstr, gmfile
  integer :: ios

  call get_command_argument(1, value=cmdstr, status=ios)
  if (ios .eq. 0) then
          read(cmdstr, *) natoms
  else
          stop "Usage: intc_mabeval.x [atoms] [geom file]"
  end if
  call get_command_argument(2, value=gmfile, status=ios)
  if (ios .ne. 0) stop "Usage: intc_mabeval.x [atoms] [geom file]"
  print "(A,A)",  "MEX Geometry file: ", trim(adjustl(gmfile))

  allocate(geom(3,natoms))
  call read_colgeom(gmfile, geom, natoms)
  call convertgeom_bohr2ang(geom, natoms)
  call evaluate_intc_mabeval(natoms, geom)
  
end program intc_mabeval
