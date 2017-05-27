!
! intc_mab_f: evaluate f for two internal coordinates given a geometry
program intc_mab_f
  use ioutil, only: read_colgeom
  use intc_mabutil, only: evaluate_intc_f, convertgeom_bohr2ang
  implicit none
  integer :: natoms
  integer :: st1, st2
  double precision, dimension(:,:), allocatable :: geom
  integer, dimension(2) :: intc
  double precision, dimension(:),allocatable :: fintc
  character(255) :: cmdstr, gmfile
  integer :: ios
    
  call get_command_argument(1, value=cmdstr, status=ios)
  if (ios .eq. 0) then
          read(cmdstr, *) natoms
  else
          stop "Usage: intc_mabeval.x [atoms]* [geom file] [int1] [int2]"
  end if
  call get_command_argument(2, value=gmfile, status=ios)
  if (ios .ne. 0) stop "Usage: intc_mabeval.x [atoms] [geom file]* [int1] [int2]"
  print "(A,A)",  "MEX Geometry file: ", trim(adjustl(gmfile))
  call get_command_argument(3, value=cmdstr, status=ios)
  if (ios .eq. 0) then
          read(cmdstr, *) intc(1)
  else
          stop "Usage: intc_mabeval.x [atoms] [geom file] [int1]* [int2]"
  end if
  call get_command_argument(4, value=cmdstr, status=ios)
  if (ios .eq. 0) then
          read(cmdstr, *) intc(2)
  else
          stop "Usage: intc_mabeval.x [atoms] [geom file] [int1] [int2]*"
  end if

  allocate(geom(3,natoms))
  allocate(fintc(3*natoms-6))
  call read_colgeom(gmfile, geom, natoms)

  call read_inputfile(st1, st2)

  call evaluate_intc_f(natoms, geom, fintc, st1, st2)
  print "(a,f15.8,a,f15.8)", 'Fij_x, Fij_y = ',fintc(intc(1)),',',fintc(intc(2))
  deallocate(geom, fintc)

contains
  subroutine read_inputfile(state1, state2)
    use ioutil, only: get_unit
    implicit none
    integer, intent(out) :: state1, state2
    integer :: flu, ios
    namelist /general/ state1, state2
    state1 = 1
    state2 = 2
    flu = get_unit()
    open(file="intc_mab_f.in",unit=flu,status="unknown",action="read",iostat=ios)
    if (ios .eq. 0) then
            read(unit=flu,nml=general)
            close(flu)
    end if
    return
  end subroutine read_inputfile
    
    
end program intc_mab_f
