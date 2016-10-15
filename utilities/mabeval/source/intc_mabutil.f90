module intc_mabutil
  implicit none
  integer, dimension(594) :: intctype
  character(4), dimension(8) :: tipus
  double precision, parameter :: MATH_PI = 3.14159265358979323846264338328d0
  double precision, parameter :: BOHR2ANG = 0.529177249
contains
  !===================================================================
  ! convertgeom_bohr2ang: convert a geometry from bohr 2 angstrom
  subroutine convertgeom_bohr2ang (geom, natoms)
    implicit none
    integer, intent(in) :: natoms
    double precision, dimension(3,natoms),intent(inout) :: geom
    integer :: i, j
    do i = 1, natoms
            do j = 1, 3
                    geom(j,i) = geom(j,i)*BOHR2ANG
            end do
    end do
    return
  end subroutine convertgeom_bohr2ang
  !-------------------------------------------------------------------

  !===================================================================
  ! evaluate_intc_f: evaluate fij in internal coordinates
  subroutine evaluate_intc_f (natoms, geom, fijintc, st1, st2)
    use ioutil, only: get_unit
    implicit none
    integer, intent(in) :: natoms, st1, st2
    double precision, dimension(3,natoms), intent(inout) :: geom
    double precision, dimension(3*natoms-6), intent(out) :: fijintc
    integer :: intcfl_unit
    integer, dimension(6,1782)  :: ibcode
    integer, dimension(20,594) :: ibcontr
    double precision, dimension(54,594) :: bmat
    double precision, dimension(600,594) :: b
    double precision, dimension(594) :: qq
    character(8), dimension(600) :: labr, labc
    integer :: natm, nst
    double precision, dimension(:,:,:), allocatable :: cg, dcg
    double precision, dimension(:,:), allocatable :: h
    double precision, dimension(:), allocatable :: fijcart
    double precision, dimension(:), allocatable :: e
    integer :: na3, nintc
    integer :: mxf
    double precision :: maxforce

    na3=3*natoms
    nintc=na3-6
    
    call initPotential()
    call getInfo(natm, nst)
    allocate(cg(na3, nst, nst))
    allocate(dcg(na3, nst, nst))
    allocate(e(nst))
    allocate(h(nst,nst))
    allocate(fijcart(na3))
    call EvaluateSurfgen(geom, e, cg, h, dcg)
    fijcart(1:na3) = (-1d0)*cg(1:na3, st1, st2)/BOHR2ANG

    intcfl_unit = get_unit()
    call convertgeom_bohr2ang(geom, natoms)
    call open_intcfl_file(intcfl_unit)
    call read_intcfl(intcfl_unit, natoms, ibcode, ibcontr)
    call generate_bmat(natoms, geom, ibcode, ibcontr, bmat, qq)
    call make_complete_bmat(natoms, bmat, labr, labc, ibcode, ibcontr, b)
    call prblkc('B matrix', b, 600, na3, nintc, labr, labc, 1, 6) 
    call intforc(nintc,na3,bmat,ibcontr,fijcart,fijintc,qq,mxf,maxforce)
    return
    
  end subroutine evaluate_intc_f
  !-------------------------------------------------------------------
  
  !===================================================================
  ! evaluate_intc_mabeval: evaluate mab in internal coordinates
  subroutine evaluate_intc_mabeval (natoms, geom)
    use ioutil, only: get_unit
    implicit none
    integer, intent(in) :: natoms
    double precision, dimension(3,natoms), intent(in) :: geom
    integer :: intcfl_unit
    integer, dimension(6,1782)  :: ibcode
    integer, dimension(20,594) :: ibcontr
    double precision, dimension(54,594) :: bmat
    double precision, dimension(594) :: qq
    
    intcfl_unit = get_unit()
    call open_intcfl_file(intcfl_unit)
    call read_intcfl(intcfl_unit, natoms, ibcode, ibcontr)
    call generate_bmat(natoms, geom, ibcode, ibcontr, bmat, qq)

    return
    
  end subroutine evaluate_intc_mabeval
  !-------------------------------------------------------------------

  !===================================================================
  ! generate_bmat: generate bmatrix and internal coordinate values
  ! for a geometry
  subroutine generate_bmat(natoms, geom, ibcode, ibcontr, bmat, qq)
    implicit none
    integer, intent(in) :: natoms
    double precision, dimension(3,natoms) :: geom
    integer, dimension(6,1782), intent(in) :: ibcode
    integer, dimension(6,594), intent(in) :: ibcontr
    double precision, dimension(54,594), intent(inout) :: bmat
    double precision, dimension(594), intent(inout) :: qq
    integer :: nintc
    nintc=3*natoms-6
    call machbnew(natoms, geom, nintc, .false., MATH_PI/2, ibcode, ibcontr, &
            bmat, qq)

    return
  end subroutine generate_bmat
  !-------------------------------------------------------------------

  !===================================================================
  ! make_complete_bmat: copy bmat(54,594) into b(3*na,3*na-6)
  subroutine make_complete_bmat (natoms, bmat, labr, labc, ibcode, ibcontr, b)
    implicit none
    integer, intent(in) :: natoms
    double precision, dimension(54,594), intent(in) :: bmat
    integer, dimension(6,1782), intent(in) :: ibcode
    integer, dimension(20,594), intent(in) :: ibcontr
    character(8), dimension(600), intent(out) :: labr, labc
    double precision, dimension(600,594), intent(out) :: b
    integer :: i, k
    integer :: na, ia, iat3, k3
    do i = 1, 3*natoms
            write(labr(i),'(a4,i3)')'cart',i
    end do
    do i=1,3*natoms-6
            write(labc(i),'(a4,i3)')' int',i
    end do
    
    b=0d0
    do i=1, 3*natoms-6
            na=ibcontr(2,i)
            k3=0
            do k=1, na
                    k3=k3+3
                    ia=ibcontr(k+2,i)
                    iat3=ia*3
                    b(iat3-2,i)=bmat(k3-2,i)
                    b(iat3-1,i)=bmat(k3-1,i)
                    b(iat3,i)=bmat(k3,i)
            end do
    end do
    return
  end subroutine make_complete_bmat
  !-------------------------------------------------------------------
  
  !===================================================================
  ! open_intcfl_file: open and intcfl file
  subroutine open_intcfl_file (flunit)
    implicit none
    integer, intent(in) :: flunit
    integer :: ios
    open(file="intcfl",unit=flunit,status="old",action="read",iostat=ios)
    if (ios .ne. 0) stop "*** Could not open intcfl file. ***"
    return
  end subroutine open_intcfl_file
  !-------------------------------------------------------------------

  !===================================================================
  ! print_bmatrix: print bmatrix
  subroutine print_bmatrix (na, bmat)
    implicit none
    integer, intent(in) :: na
    double precision, dimension(na*3,na*3-6), intent(in) :: bmat
    integer :: i, j
    do i=1, 3*na-6
            print *, "Int = ", i
            do j = 1, 3*na
                    print *, bmat(j,i)
            end do
    end do
    return
  end subroutine print_bmatrix
  !===================================================================
  ! read_intcfl: read an intcfl file.
  subroutine read_intcfl (flunit, natoms, ibcode, ibcontr)
    implicit none
    ! .. input scalars ..
    ! flunit = file unit of intcfl file
    ! natoms = number of atoms
    ! nintc = number of internal coordinates
    integer, intent(in) :: flunit
    integer, intent(in) :: natoms
    ! .. output arrays ..
    integer, dimension(6, 1782), intent(out) :: ibcode
    integer, dimension(20, 594), intent(out) :: ibcontr
    ! .. local scalars..
    integer :: nintc

    nintc = 3*natoms-6
    ! read first line (title).
    read(flunit,*)
    call bread(natoms, flunit, 6, nintc, nintc, ibcode, ibcontr, .false.)
    
  end subroutine read_intcfl
end module intc_mabutil
