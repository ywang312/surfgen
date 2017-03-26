program ghz_displacement
  ! Generate displacements along g, h, or z vectors (supplied).
  ! Written by cthree-40
  !
  use ioutil
  implicit none
  double precision, dimension(:,:), allocatable :: xgeom
  double precision, dimension(:),   allocatable :: gv, hv, zv
  double precision, dimension(:),   allocatable :: disp_vec
  integer :: natoms
  integer :: gfunit, hfunit, zfunit
  character(255) :: xgmfile
  character(255) :: cmdstr
  character(255) :: gvfl, hvfl, zvfl
  character(255) :: usage
  double precision :: x, y, z
  integer :: num_disps
  double precision :: total_disp
  double precision, dimension(:), allocatable :: disp_coeff
  double precision, external :: dnrm2
  integer :: ios

  ios = 0
  
  ! Process command line arguments. Correct usage of this program is:
  call process_command_line_args(4, natoms, xgmfile, total_disp, num_disps, ios)
  if (ios .ne. 0) stop "Error!"
  print "('Natoms = ',i5,' Displacement = ',f5.3,' Number of disps = ',i5)", &
          natoms, total_disp, num_disps
  print "('Geometry file = ',a)", trim(adjustl(xgmfile))
  
  ! Read input file.
  call read_input_file(x, y, z)
  print "('Displacing along: ',f5.3,'g + ',f5.3,'h + ',f5.3,'z')", x, y, z
  
  ! Read geometry file.
  allocate(xgeom(3,natoms))
  call read_colgeom(xgmfile, xgeom, natoms)
  print *, "MEX: "
  call print_colgeom2(xgeom, natoms, atom_names, atom_numbers, atom_weights)

  ! Read vector files.
  allocate(gv(3*natoms))
  allocate(hv(3*natoms))
  allocate(zv(3*natoms))
  gvfl = 'gvector.dat'
  hvfl = 'hvector.dat'
  zvfl = 'zvector.dat'
  gfunit = 10
  hfunit = 11
  zfunit = 12
  call read_vector(gv, natoms, gfunit, gvfl, 1)
  if (gfunit .eq. -1) stop "*** Error! G vector not found! ***"
  call read_vector(hv, natoms, hfunit, hvfl, 1)
  if (hfunit .eq. -1) stop "*** Error! H vector not found! ***"
  call read_vector(zv, natoms, zfunit, zvfl, 1)
  if (zfunit .eq. -1) then
          print *, "*** Warning! Z vector not found! ***"
          print *, " ...setting z = 0."
          z = 0
  end if

  ! Make displacement vector
  allocate(disp_vec(3*natoms))
  call make_displacement_vector(natoms, gv, x, hv, y, zv, z, disp_vec)
  print *, "Displacement vector: "
  print "(5x,3f14.8)", disp_vec
  print "(2x,'Norm =',f14.8)", dnrm2(3*natoms,disp_vec,1)

  ! Generate displacement coefficients
  allocate(disp_coeff(num_disps + 1)) ! "+1 for 0th displacement"
  call make_displacement_coeff_list(disp_coeff, total_disp, num_disps)

  ! Make displacements, printing geometries.
  call make_vector_displacements(xgeom, natoms, disp_vec, disp_coeff, &
          num_disps)
  
  ! Deallocate arrays
  deallocate(disp_vec, xgeom, disp_coeff, gv, hv, zv)
contains

  !
  ! make_displacement_coeff_list: generate list of coefficients for
  ! displacements.
  !
  subroutine make_displacement_coeff_list(dc_list, total_disp, ndisps)
    implicit none
    integer, intent(in) :: ndisps
    double precision, intent(in) :: total_disp
    double precision, dimension(ndisps+1), intent(out) :: dc_list
    double precision :: dincr ! displacement increment
    integer :: i
    
    dincr = total_disp / dble(ndisps)
    do i = 1, ndisps + 1
            dc_list(i) = dble(i - 1) * dincr
    end do
    return
  end subroutine make_displacement_coeff_list
  
  !
  ! make_displacement_vector: make vector for displacements given g, h, and z
  ! vectors and the coeffiecients for each to form the vector we are displacing
  ! along.
  !
  subroutine make_displacement_vector(natoms, gv, gcoeff, hv, hcoeff, zv, &
          zcoeff, disp_vec)
    implicit none
    integer, intent(in) :: natoms
    double precision, intent(in) :: gcoeff, hcoeff, zcoeff
    double precision, dimension(3*natoms), intent(in) :: gv, hv, zv
    double precision, dimension(3*natoms), intent(out):: disp_vec
    integer :: i
    double precision :: dvnorm
    double precision, external :: dnrm2
    
    do i = 1, 3*natoms
            disp_vec(i) = gv(i)*gcoeff + hv(i)*hcoeff + zv(i)*zcoeff
    end do
    dvnorm = dnrm2(3*natoms, disp_vec, 1)
    do i = 1, 3*natoms
            disp_vec(i) = disp_vec(i)/dvnorm
    end do
    return
  end subroutine make_displacement_vector

  !
  ! make_vector_displacements: make displacements from an original geometry
  ! along some vector. Each new geometry is printed to a new file with
  ! the number of the displacement appended.
  !
  subroutine make_vector_displacements(geom0, natoms, dvec, dcoef, ndisps)
    implicit none
    integer, intent(in) :: natoms, ndisps
    double precision, dimension(ndisps+1), intent(in) :: dcoef
    double precision, dimension(3,natoms), intent(in) :: dvec, geom0
    double precision, dimension(3,natoms) :: geomx
    integer :: i
    character(255) :: flname, flindex
    integer :: flun, ios
    
    do i = 1, ndisps+1
            write(flindex,"(i5)") (i-1)
            flname = "geometry."//trim(adjustl(flindex))
            flun = get_unit()
            open(file=trim(adjustl(flname)), unit=flun, status='unknown', &
                    iostat = ios, action='write', position = 'rewind')
            if (ios .ne. 0) then
                    print *, "Could not write to file ", trim(adjustl(flname))
                    stop "Exiting..."
            end if

            geomx = geom0 + dvec*dcoef(i)
            call print_colgeom2_fileout(geomx, natoms, atom_names, &
                    atom_numbers, atom_weights, flun)
            close(flun)
    end do
    return
  end subroutine make_vector_displacements
    
  !
  ! process_command_line_args: process the command line arguments. Return -1 if
  ! incorrect number of arguments. +n for error reading nth argument.
  !
  subroutine process_command_line_args(nargs, natoms, xgmfile, tdisp, ndisp, &
          error)
    implicit none
    integer, intent(in) :: nargs
    integer, intent(out) :: natoms, ndisp, error
    double precision, intent(out) :: tdisp
    character(255), intent(out) :: xgmfile
    character(255) :: usage, cmdstr
    integer :: ios

    usage="ghz_displacements.x"
    usage=trim(adjustl(usage))//" [natoms] [xgeom file] [total disp]"
    usage=trim(adjustl(usage))//" [num disps]"

    call get_command_argument(1, value=cmdstr, status=ios)
    if (ios .eq. 0) then
            read(cmdstr, *) natoms
    else
            print *, trim(adjustl(usage))
            error = 1
            return
    end if

    call get_command_argument(2, value=xgmfile, status=ios)
    if (ios .ne. 0) then
            print *, trim(adjustl(usage))
            error = 2
            return
    end if

    call get_command_argument(3, value=cmdstr, status=ios)
    if (ios .eq. 0) then
            read(cmdstr, *) tdisp
    else
            print *, trim(adjustl(usage))
            error = 3
            return
    end if

    call get_command_argument(4, value=cmdstr, status=ios)
    if (ios .eq. 0) then
            read(cmdstr, *) ndisp
    else
            print *, trim(adjustl(usage))
            error = 4
            return
    end if
    error = 0
    return
  end subroutine process_command_line_args

  !
  ! read_input_file: read input file, returning x, y, and z values.
  !
  subroutine read_input_file(x, y, z)
    implicit none
    double precision, intent(out) :: x, y, z
    integer :: un, ios
    character(255) :: flname
    namelist /general/ x, y, z

    flname = "ghz_disp.in"
    un = get_unit()
    open(unit=un,file=trim(adjustl(flname)),status='old',iostat=ios)
    if (ios .ne. 0) stop "*** Could not open input file: ghz_disp.in ***"
    read(unit=un, nml=general)
    close(unit=un)
    return
  end subroutine read_input_file
    
end program ghz_displacement
