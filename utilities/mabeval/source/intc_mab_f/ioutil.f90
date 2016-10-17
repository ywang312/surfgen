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
  ! read_colgeom: read geometry in columbus format
  subroutine read_colgeom (flnm, gm, na)
    implicit none
    character(255), intent(in) :: flnm
    integer, intent(in) :: na
    double precision, dimension(3,na), intent(out) :: gm
    integer :: i, flun, ios
    
    ! Open file
    flun = get_unit()
    open(unit=flun,file=flnm,status="old",iostat=ios)
    if (ios .ne. 0) stop "*** Can't open input geometry file. ***"
    ! Read geometry
    do i = 1, na
            read(flun,fmt="(11x,3f14.8,14x)") gm(1:3,i)
    end do
    close(flun)
    return
  end subroutine read_colgeom

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
