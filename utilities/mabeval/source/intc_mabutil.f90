module intc_mabutil
  implicit none
  integer, dimension(594) :: intctype
  character(4), dimension(8) :: tipus
  double precision, parameter :: MATH_PI = 3.14159265358979323846264338328d0
  double precision, parameter :: BOHR2ANG = 0.529177249
  double precision, parameter :: AU2CM1 = 219474.63
contains

  !===================================================================
  ! compute_ftheta: compute f(theta) for circulation
  function compute_ftheta (qq_0, gm_ang_0, c1, c2, eps, theta, qq, geom_bohr,&
          geom_ang, prevfij, nstates, natoms, bmat, cg, dcg, h, e, st1, st2 ,   &
          ibcode, ibcontr) &
          result(fth)
    implicit none
    integer, intent(in) :: c1, c2
    integer, intent(in) :: st1, st2
    integer, intent(in) :: nstates, natoms
    double precision, intent(in) :: eps, theta
    double precision, dimension(3,natoms), intent(in) :: gm_ang_0
    double precision, dimension(594), intent(in) :: qq_0
    integer, dimension(6,1782), intent(in) :: ibcode
    integer, dimension(20,594), intent(in) :: ibcontr
    
    double precision, dimension(594), intent(inout) :: qq
    double precision, dimension(54,594), intent(inout) :: bmat
    double precision, dimension(3,natoms), intent(inout) :: geom_bohr, geom_ang
    double precision, dimension(3*natoms,nstates,nstates), intent(inout) :: cg
    double precision, dimension(3*natoms,nstates,nstates), intent(inout) :: dcg
    double precision, dimension(nstates,nstates), intent(inout) :: h
    double precision, dimension(nstates), intent(inout) :: e
    double precision, dimension(3*natoms), intent(inout) :: prevfij

    double precision, dimension(594) :: qq_scr
    double precision, dimension(54,594) :: bmat_scr
    double precision, dimension(3*natoms) :: fijcart, nfij
    double precision, dimension(3*natoms-6) :: fijintc

    integer :: mxf
    double precision :: maxforce
    
    double precision, external :: ddot
    double precision :: dotp

    double precision :: fth
    integer :: na3, nintc

    na3 = 3*natoms
    nintc = 3*natoms - 6
    ! Set internals to origin. Add displacement in internals.
    qq = 0d0
    qq(c1) = eps * dcos(theta)
    qq(c2) = eps * dsin(theta)
    qq_scr = 0d0

    bmat_scr = bmat
    ! Distort geometry. Convert geometry to bohr.
    call distort(na3, nintc, natoms, qq, qq_0, qq_scr, gm_ang_0, ibcode, &
            ibcontr, bmat_scr, MATH_PI/2, geom_ang)
    geom_bohr = geom_ang
    call convertgeom_ang2bohr(geom_bohr, natoms)

    ! Evaluate surface at new geometry (bohrs). Convert fij to Angstrom for
    ! cart2int usage. Align fij phase. Generate bmat. Trasform fij.
    call EvaluateSurfgen(geom_bohr, e, cg, h, dcg)
    fijcart(1:na3) = (-1d0)*cg(1:na3, st1, st2)/BOHR2ANG
    !nfij = (-1d0)*fijcart
    !if (ddot(na3,nfij,1,prevfij,1) .gt. ddot(na3,fijcart,1,prevfij,1)) then
    !        fijcart = nfij
    !end if
    dotp = ddot(na3, fijcart, 1, prevfij, 1)
    fijcart = fijcart * dotp / abs(dotp)
    !prevfij = fijcart
    call generate_bmat(natoms, geom_ang, ibcode, ibcontr, bmat_scr, qq_scr)
    call intforc(nintc, na3, bmat_scr, ibcontr, fijcart, fijintc, qq_scr, &
            mxf, maxforce)

    fth = (-1d0)*fijintc(c1)*dsin(theta) + fijintc(c2)*dcos(theta)
    fth = fth * eps
    if (abs(fth) .gt. 1d-12) prevfij = fijcart
    return
  end function compute_ftheta
    
  !===================================================================
  ! convertgeom_ang2bohr: convert a geometry from angstrom 2 bohr
  subroutine convertgeom_ang2bohr (geom, natoms)
    implicit none
    integer, intent(in) :: natoms
    double precision, dimension(3,natoms),intent(inout) :: geom
    integer :: i, j
    do i = 1, natoms
            do j = 1, 3
                    geom(j,i) = geom(j,i)/BOHR2ANG
            end do
    end do
    return
  end subroutine convertgeom_ang2bohr
  !-------------------------------------------------------------------

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

    if (st1 .gt. nst .or. st2 .gt. nst) then
            print "('Error! NSTATES=',i3,' Check input.')", nst
            return
    end if

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
   ! call make_complete_bmat(natoms, bmat, labr, labc, ibcode, ibcontr, b)
   ! call prblkc('B matrix', b, 600, na3, nintc, labr, labc, 1, 6) 
    call intforc(nintc,na3,bmat,ibcontr,fijcart,fijintc,qq,mxf,maxforce)
    return
    
  end subroutine evaluate_intc_f
  !-------------------------------------------------------------------
  
  !===================================================================
  ! evaluate_intc_mabeval: evaluate mab in internal coordinates
  subroutine evaluate_intc_mabeval (natoms, geom, c1, c2, eps, maxi, &
          circdir, rots, st1, st2)
    use ioutil, only: get_unit
    implicit none
    integer, intent(in) :: natoms
    integer, intent(in) :: c1, c2
    integer, intent(in) :: st1, st2
    integer, intent(in) :: maxi, circdir, rots
    double precision, intent(in) :: eps ! radius of loop
    double precision, dimension(3,natoms), intent(inout) :: geom
    double precision, dimension(3,natoms) :: geom_ang
    integer :: intcfl_unit
    integer, dimension(6,1782)  :: ibcode
    integer, dimension(20,594) :: ibcontr
    double precision, dimension(54,594) :: bmat
    double precision, dimension(594) :: qq
    double precision :: c1_0, c2_0 ! Origin coordinates (c1, c2)
    integer :: natm, nst

    ! Initialize surface
    call initPotential()
    call getInfo(natm, nst)

    ! Read in intcfl and generate initial internal coordinates.
    intcfl_unit = get_unit()
    geom_ang = geom
    call convertgeom_bohr2ang(geom_ang, natoms)
    call open_intcfl_file(intcfl_unit)
    call read_intcfl(intcfl_unit, natoms, ibcode, ibcontr)
    bmat = 0d0
    qq = 0d0
    call generate_bmat(natoms, geom_ang, ibcode, ibcontr, bmat, qq)
    c1_0 = qq(c1)
    c2_0 = qq(c2)
    print "('Origin: (',f6.3,',',f6.3,')')", c1_0, c2_0
    call intc_mabanalysis(natoms, nst, ibcode, ibcontr, bmat, qq, &
            c1_0, c2_0, eps, c1, c2, maxi, circdir, rots,   &
            geom_ang, st1, st2)
    return
    
  end subroutine evaluate_intc_mabeval
  !-------------------------------------------------------------------

  !-------------------------------------------------------------------
  ! find_degeneracy: locate degeracy in energy array and return states
  subroutine find_degeneracy(e, nst, st1, st2, err)
    implicit none
    integer, intent(in) :: nst
    double precision, dimension(nst), intent(in) :: e
    integer, intent(out) :: st1, st2, err
    integer :: i

    err = 0
    do i = 2, nst
            if (abs(e(i-1) - e(i)) < 1d-5) then
                    err = 1
                    st1 = i-1
                    st2 = i
                    return
            end if
    end do
    return
  end subroutine find_degeneracy

  
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
  ! intc_mabanalysis: perform mab analysis
  subroutine intc_mabanalysis (natm, nst, ibcode, ibcontr, bmat, qq_0, &
          c1_0, c2_0, epsilon, c1, c2, maxi, circdir, rots, geom_ang_0,&
          st1, st2)
    implicit none
    integer, intent(in) :: natm, nst, st1, st2
    integer, intent(in) :: c1, c2
    integer, intent(in) :: maxi, circdir, rots
    integer, dimension(6,1782), intent(in) :: ibcode
    integer, dimension(20,594), intent(in) :: ibcontr
    double precision, dimension(54,594), intent(inout) :: bmat
    double precision, dimension(594), intent(in) :: qq_0
    double precision, intent(in) :: c1_0, c2_0
    double precision, intent(in) :: epsilon
    double precision, dimension(3,natm), intent(in) :: geom_ang_0

    double precision :: intgrl
    double precision :: fija, fijb, fijx
    double precision :: ta, tb, tx
    double precision :: theta ! circulation angle
    double precision :: dth ! d(theta)
    double precision, dimension(3,natm) :: geom_ang, geom_bohr
    double precision, dimension(594) :: geom_qq

    double precision, dimension(3*natm) :: prevfij

    double precision, dimension(3*natm, nst, nst) :: cg, dcg
    double precision, dimension(nst, nst) :: h
    double precision, dimension(nst) :: e

    integer :: i
    integer :: na3

    na3 = 3*natm
    
    ! Set up circulation.
    intgrl = 0d0
    dth = 2*MATH_PI/maxi
    if (circdir .lt. 0) dth = dth * (-1d0)
    ta = 0d0
    prevfij=1d0
    geom_bohr = 0d0
    geom_qq = 0d0
    
    ! Circulation
    fija = compute_ftheta(qq_0, geom_ang_0, c1, c2, epsilon, ta, geom_qq, &
            geom_bohr, geom_ang, prevfij, nst, natm, bmat, cg, dcg, h, e, &
            st1, st2, ibcode, ibcontr)
    print *, "fija(0) = ", fija
    do i = 1, maxi*rots
            print *, "---------------------------------------------"
            print *, "Iteration = ", i
            tb = ta + dth
            tx = (ta + tb) / 2d0
            fijx = compute_ftheta(qq_0, geom_ang_0, c1, c2, epsilon, tx,   &
                    geom_qq, geom_bohr, geom_ang, prevfij, nst, natm, bmat,&
                    cg, dcg, h, e, st1, st2, ibcode, ibcontr)
            print *, "fij(x) = ", fijx
            fijb = compute_ftheta(qq_0, geom_ang_0, c1, c2, epsilon, tb,   &
                    geom_qq, geom_bohr, geom_ang, prevfij, nst, natm, bmat,&
                    cg, dcg, h, e, st1, st2, ibcode, ibcontr)
            intgrl = intgrl + simpsons_rule(fija, fijx, fijb, ta, tx, tb)
            ta = tb
            fija = fijb
    end do
    print *, "integral = ", intgrl
    return
  end subroutine intc_mabanalysis
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
  !-------------------------------------------------------------------
  
  !===================================================================
  ! simpsons_rule: compute SF(x)dx using simpsons rule.
  !  SF(x)dx = [(b-a)/6]*[f(a) + 4f((a+b)/2) + f(b)]
  function simpsons_rule(fa, fx, fb, a, x, b) result(val)
    implicit none

    ! ..input scalars..
    double precision, intent(in) :: fa, fx, fb
    double precision, intent(in) :: a, x, b

    double precision :: val

    val = fa + 4*fx + fb
    val = val * (b - a)
    val = val / 6
    return
  end function simpsons_rule
  !-------------------------------------------------------------------
  
end module intc_mabutil
