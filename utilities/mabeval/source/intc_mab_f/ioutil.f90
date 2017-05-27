module ioutil
  !-------------------------------------------------------------------
  ! Input/Output utilities
  !
  ! cthree-40
  ! Yarkony Group
  ! Department of Chemistry
  ! The Johns Hopkins University
  !
  ! -- Change Log --
  ! 2016-07-21 Created.
  ! 2016-08-31 Added: get_unit, read_colgeom, read_vector
  !-------------------------------------------------------------------
  implicit none
  character(len=3), dimension(:), allocatable :: atom_names
  double precision, dimension(:), allocatable :: atom_weights, atom_numbers
contains

  !-------------------------------------------------------------------
  ! get_unit: get available file unit
  function get_unit () result(u)
    implicit none
    integer :: i
    integer :: u
    logical :: uexist, uopen
    u = 0
    do i = 15, 9999
            inquire(unit=i,exist=uexist,opened=uopen)
            if (uexist .and. .not. uopen) then
                    u=i
                    exit
            end if
    end do
    if (u .eq. 0) stop "No unit available."
  end function get_unit

  !-------------------------------------------------------------------
  ! print_colgeom: print columbus geometry
  subroutine print_colgeom (g, na)
    implicit none

    ! ..input scalars..
    integer, intent(in) :: na

    ! ..input arrays..
    double precision, dimension(3,na), intent(in) :: g

    ! ..local scalars..
    integer :: i

    do i = 1, na
            print "(3f15.8)", g(1:3,i)
    end do
    return
  end subroutine print_colgeom

  !-------------------------------------------------------------------
  ! print_colgeom2: print columbus geometry with atom info
  subroutine print_colgeom2 (g, na, name, num, wt)
    implicit none

    ! ..input scalars..
    integer, intent(in) :: na

    ! ..input arrays..
    double precision, dimension(3,na), intent(in) :: g
    double precision, dimension(na),   intent(in) :: num
    double precision, dimension(na),   intent(in) :: wt
    character(len=3), dimension(na),   intent(in) :: name
    
    ! ..local scalars..
    integer :: i

    do i = 1, na
            print "(1x,1a,f7.1,4f14.8)", name(i), num(i), g(1,i), g(2,i), &
                    g(3,i), wt(i)
    end do
    return
  end subroutine print_colgeom2

  !-------------------------------------------------------------------
  ! read_colgeom: read geometry in columbus format
  subroutine read_colgeom (flnm, gm, na)
    implicit none
    character(255), intent(in) :: flnm
    integer, intent(in) :: na
    double precision, dimension(3,na), intent(out) :: gm
    integer :: i, flun, ios

    if (.not. allocated(atom_weights)) allocate(atom_weights(na))
    if (.not. allocated(atom_numbers)) allocate(atom_numbers(na))
    if (.not. allocated(atom_names)) allocate(atom_names(na))
    ! Open file
    flun = get_unit()
    open(unit=flun,file=flnm,status="old",iostat=ios)
    if (ios .ne. 0) stop "*** Can't open input geometry file. ***"
    ! Read geometry
    do i = 1, na
            read(flun,fmt="(1x,a,f7.1,f14.8,f14.8,f14.8,f14.8)") &
                    atom_names(i),atom_numbers(i),gm(1,i),gm(2,i),gm(3,i), &
                    atom_weights(i)
    end do
    close(flun)
    return
  end subroutine read_colgeom

  !-------------------------------------------------------------------
  ! read_colgeom2: read geometry in columbus format
  subroutine read_colgeom2 (flnm, gm, na, masses, atoms, nums)
    implicit none
    character(255), intent(in) :: flnm
    integer, intent(in) :: na
    double precision, dimension(3,na), intent(out) :: gm
    character(2), dimension(na), intent(out) :: atoms
    double precision, dimension(na), intent(out) :: masses, nums
    integer :: i, flun, ios
    
    ! Open file
    flun = get_unit()
    open(unit=flun,file=flnm,status="old",iostat=ios)
    if (ios .ne. 0) stop "*** Can't open input geometry file. ***"
    ! Read geometry, atomic number, and mass
    do i = 1, na
            read(flun,fmt=*) atoms(i), nums(i), &
                    gm(1,i), gm(2,i), gm(3,i), masses(i)
    end do
    close(flun)
    return
  end subroutine read_colgeom2

  !---------------------------------------------------------
  ! read_vector: read vector from file
  subroutine read_vector(v, na, u, nm, nv)
    implicit none
    integer, intent(in) :: na, nv
    integer, intent(inout) :: u
    character(255), intent(in) :: nm
    double precision, dimension(3,na,nv), intent(inout) :: v
    integer :: i, j

    open(unit=u,file=trim(adjustl(nm)),action="read",position="rewind",&
            status="old",iostat=i)
    if (i .ne. 0) then
            u=-1
            return
    end if
    read(unit=u,fmt=*) v
    close(unit = u)
    return
  end subroutine read_vector
    
end module ioutil
