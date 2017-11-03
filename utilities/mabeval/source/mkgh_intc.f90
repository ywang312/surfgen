! mkgh_intc: this routine builds g and h and transforms vectors
! into internal coordinates.
program mkgh_intc
  use ioutil, only: read_colgeom
  use intc_mabutil, only: get_ghvec_intc, convertgeom_bohr2ang
  implicit none
  integer :: natoms
  integer :: st1, st2
  double precision, dimension(:,:), allocatable :: geom
  double precision, dimension(:),allocatable :: gintc, hintc
  logical :: mw
  character(255) :: cmdstr, gmfile
  integer :: ios
  integer :: i
  
  call get_command_argument(1, value=cmdstr, status=ios)
  if (ios .eq. 0) then
          read(cmdstr, *) natoms
  else
          stop "Usage: intc_mabeval.x [atoms]* [geom file]"
  end if
  call get_command_argument(2, value=gmfile, status=ios)
  if (ios .ne. 0) stop "Usage: intc_mabeval.x [atoms] [geom file]*"
  print "(A,A)",  "MEX Geometry file: ", trim(adjustl(gmfile))

  allocate(geom(3,natoms))
  call read_colgeom(gmfile, geom, natoms)
  call read_inputfile(st1, st2, mw)

  allocate(gintc(3*natoms-6))
  allocate(hintc(3*natoms-6))
  call get_ghvec_intc(natoms,st1,st2,geom,gintc,hintc,mw)
  print *, "GVECTOR:"
  do i = 1, 3*natoms-6
          print "(f13.8)", gintc(i)
  end do
  print *, "HVECTOR:"
  do i = 1, 3*natoms-6
          print "(f13.8)", hintc(i)
  end do
contains
    subroutine read_inputfile(state1, state2, mweight)
    use ioutil, only: get_unit
    implicit none
    integer, intent(out) :: state1, state2
    logical, intent(out) :: mweight
    integer :: flu, ios
    namelist /general/ state1, state2, mweight
    state1 = 1
    state2 = 2
    mweight = .true.
    flu = get_unit()
    open(file="mkgh_intc.in",unit=flu,status="unknown",action="read",iostat=ios)
    if (ios .eq. 0) then
            read(unit=flu,nml=general)
            close(flu)
    end if
    return
  end subroutine read_inputfile
    
    
end program mkgh_intc
