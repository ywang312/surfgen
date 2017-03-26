module mabutil
  !-------------------------------------------------------------------
  ! Utilities for evalating the molecular aharanov-bohm effect.
  !
  ! cthree-40
  ! Yarkony Group
  ! Department of Chemistry
  ! The Johns Hopkins University
  !
  ! -- Change Log --
  ! 2016-07-21 - Created.
  !-------------------------------------------------------------------
  use ioutil, only: atom_weights, atom_numbers, atom_names
  implicit none
  double precision :: MATH_PI = 4d0*datan(1d0)
  double precision :: AU2CM1 = 219474.63
contains
  !-------------------------------------------------------------------
  ! addzvec: add a vector to a geometry
  subroutine addzvec (zcoeff, zvector, lv, geometry)
    implicit none

    ! ..input scalars..
    integer, intent(in) :: lv
    double precision, intent(in) :: zcoeff

    ! ..input arrays..
    double precision, dimension(lv), intent(in) :: zvector

    ! ..input/output arrays..
    double precision, dimension(lv), intent(inout) :: geometry

    if (zcoeff**2 .le. 1d-16) return
    geometry = geometry + zcoeff * zvector
    return
  end subroutine addzvec

  !-------------------------------------------------------------------
  ! align_phase: align the phase of one vector to another.
  subroutine align_phase (v1, v2, lv, phase)
    implicit none

    ! ..input scalars..
    integer, intent(in) :: lv

    ! ..input arrays..
    double precision, dimension(lv), intent(in) :: v1

    ! ..input/output scalars..
    double precision, intent(inout) :: phase

    ! ..input/output arrays..
    double precision, dimension(lv), intent(inout) :: v2

    ! ..local scalars..
    double precision :: dp1, dp2
    double precision, parameter :: p1=1d0, p2=-1d0

    ! ..local arrays..
    double precision, dimension(lv) :: tmp
    ! ..external functions..
    double precision, external :: ddot

    tmp = v2 * p1
    dp1 = ddot(lv, v1, 1, tmp, 1)
    tmp = v2 * p2
    dp2 = ddot(lv, v1, 1, tmp, 1)
    if (dp1 .ge. dp2) then
            phase = p1
    else
            phase = p2
    end if
    v2 = v2 * phase
    return
  end subroutine align_phase

  !-------------------------------------------------------------------
  ! align_thetafij_phase: align the phase of the coupling dependent on theta.
  subroutine align_thetafij_phase (ath, bth, phase)
    implicit none
    double precision, intent(in) :: ath
    double precision, intent(inout) :: bth, phase
    double precision :: df1, df2
    double precision :: b1, b2

    b1 = bth * 1d0
    b2 = bth * (-1d0)
    df1 = abs(ath - b1)
    df2 = abs(ath - b2)
    if (df1 .le. df2) then
            phase = 1d0
    else
            phase = -1d0
    end if
    bth = bth * phase
    return
  end subroutine align_thetafij_phase

  !-------------------------------------------------------------------
  ! build_gvector: builds g vector
  ! g_ij = (g_jj - g_ii) / 2
  subroutine build_gvector(cg, na, nst, st1, st2, gv)
    implicit none
    ! ..input scalars..
    ! na = number of atoms
    ! st1 = state 1
    ! st2 = state 2
    ! nst = number of states
    integer, intent(in) :: na, st1, st2, nst

    ! ..input arrays..
    ! cg = gradients
    double precision, dimension(3*na,nst,nst), intent(in) :: cg

    ! ..output arrays..
    ! gv = g vector
    double precision, dimension(3*na), intent(out) :: gv

    gv(1:3*na) = cg(1:3*na,st2,st2)
    gv(1:3*na) = gv(1:3*na) - cg(1:3*na,st1,st1)
    gv(1:3*na) = gv(1:3*na) / 2
    return
  end subroutine build_gvector

  !-------------------------------------------------------------------
  ! build_hvector: builds g vector
  ! h_ij = g_ij
  subroutine build_hvector(cg, na, nst, st1, st2, hv, e)
    implicit none
    ! ..input scalars..
    ! na = number of atoms
    ! st1 = state 1
    ! st2 = state 2
    ! nst = number of states
    integer, intent(in) :: na, st1, st2, nst

    ! ..input arrays..
    ! cg = gradients
    ! e  = energies
    double precision, dimension(3*na,nst,nst), intent(in) :: cg
    double precision, dimension(st2), intent(in) :: e
    
    ! ..output arrays..
    ! hv = h vector
    double precision, dimension(3*na), intent(out) :: hv

    hv(1:3*na) = cg(1:3*na,st1,st2) * abs(e(st2) - e(st1))
    return
  end subroutine build_hvector

  !-------------------------------------------------------------------
  ! build_gvector: builds s vector
  ! g_ij = (g_jj + g_ii) / 2
  subroutine build_svector(cg, na, nst, st1, st2, sv)
    implicit none
    ! ..input scalars..
    ! na = number of atoms
    ! st1 = state 1
    ! st2 = state 2
    ! nst = number of states
    integer, intent(in) :: na, st1, st2, nst

    ! ..input arrays..
    ! cg = gradients
    double precision, dimension(3*na,nst,nst), intent(in) :: cg

    ! ..output arrays..
    ! sv = s vector
    double precision, dimension(3*na), intent(out) :: sv

    sv(1:3*na) = cg(1:3*na,st2,st2)
    sv(1:3*na) = sv(1:3*na) + cg(1:3*na,st1,st1)
    sv(1:3*na) = sv(1:3*na) / 2
    return
  end subroutine build_svector
  
  !-------------------------------------------------------------------
  ! build_tranrot: build translations and rotations vector
  subroutine build_tranrot(xgeom, natoms, trot)
    implicit none

    ! ..input scalars..
    integer, intent(in) :: natoms

    ! ..input arrays..
    double precision, dimension(natoms*3), intent(in) :: xgeom

    ! ..input/output arrays..
    double precision, dimension(natoms*3, 6) :: trot

    ! ..local scalars..
    integer :: i, j

    ! ..local arrays..
    double precision, dimension(natoms*3):: cgeom
    double precision, dimension(3):: center
    ! ..external functions..
    double precision, external :: dnrm2

    ! Build rotation and translation vectors: trot. Orthogonalize vector to
    ! translation and rotation vectors. Renormalize vector.
    trot = 0d0
    center = 0d0
    do i = 1, natoms
            center = center + xgeom(i*3-2:i*3)
    end do
    center=center/natoms
    do i = 1, natoms
            cgeom(i*3-2:i*3)=xgeom(i*3-2:i*3)-center
    end do
    do i = 1, 3
            do j = 0, natoms - 1
                    trot(j*3 + i, i) = 1d0
            end do
            trot(:,i) = trot(:,i)/dnrm2(3*natoms, trot(:,i), 1)
    end do
    do i = 1, natoms
            trot(i*3-1,4) =cgeom(i*3)
            trot(i*3,4) = -cgeom(i*3-1)
    end do
    do i = 1, natoms
            trot(i*3-2,5)=-cgeom(i*3)
            trot(i*3,5) =  cgeom(i*3-2)
    end do
    do i = 1, natoms
            trot(i*3-2,6) = -cgeom(i*3-1)
            trot(i*3-1,6) = cgeom(i*3-2)
    end do
    do i = 4, 6
            trot(:,i) = trot(:,i)/dnrm2(3*natoms, trot(:,i), 1)
    end do
    ! Orthogonalize
    do i = 2, 6
            call orthogonalize_vector(trot(1,1), (3*natoms), i-1, trot(1,i))
            call normalize_vector(trot(1,i), (3*natoms))
    end do
    return
  end subroutine build_tranrot
  
  !-------------------------------------------------------------------
  ! build_zspace: build 3*NATOMS-6 = N orthogonal vectors.
  subroutine build_zspace(natoms, xgeom, lv, gv, hv, zspace, remtr)
    use ioutil, only: get_unit
    implicit none

    ! ..input scalars..
    ! natoms = number of atoms
    ! lv = natoms * 3
    integer, intent(in) :: natoms, lv

    ! ..input logicals..
    logical, intent(in) :: remtr
    
    ! ..input arrays..
    double precision, dimension(lv), intent(in) :: xgeom, gv, hv

    ! ..input/output arrays..
    double precision, dimension(lv, lv-6), intent(inout) :: zspace

    ! ..local scalars..
    integer :: mxg, mxh
    integer :: i, j
    integer :: ios, nv, zflunit
    
    ! ..local arrays..
    double precision, dimension(lv, 6) :: trotv

    double precision, external :: ddot
    
    ! Initialize the zspace vectors. g and h are vectors 1 and 2, respectively.
    ! The 3*NATOMS-8 vectors are unit vectors. The two largest components of
    ! g and h must be found to minimize overlap with these vectors.
    zspace = 0d0
    zspace(1:lv,1) = gv(1:lv)
    zspace(1:lv,2) = hv(1:lv)
    mxg = find_maximum_abs(zspace(1,1), lv)
    mxh = find_maximum_abs(zspace(1,2), lv)
    print *, "Max G element: ", mxg, zspace(mxg, 1)
    print *, "Max H element: ", mxh, zspace(mxh, 2)
    i=3
    j=1
    do while (i .le. lv-6)
            if (j .ne. mxg .and. j .ne. mxh) then
                    zspace(j,i) = 1d0
                    i=i+1
            end if
            j=j+1
    end do
    zspace(lv-6:lv,lv-6)=1d0
    call build_tranrot(xgeom, natoms, trotv)
    zflunit = get_unit()
    open(unit=zflunit,file="zvector.dat",status="old",action="read",iostat=ios)
    if (ios .eq. 0) then
            loop1: do i=3, natoms*3-6
                    loop2: do j = 1, natoms
                            read(unit=zflunit,fmt="(3f15.10)",iostat=ios) &
                                    zspace(j*3-2:j*3,i)
                            if (ios .lt. 0) exit loop1 
                    end do loop2
            end do loop1
            close(unit=zflunit)
    end if
    print *, "Read in ", max(0,i-3), " z vectors."
    do j = 3, i-1
            print "(3f15.10)", zspace(:,j)
    end do
    i=3
    do while (i .le. lv-6)
            call orthogonalize_vector(zspace(1,1), lv, i-1, zspace(1,i))
            if (.not. remtr) then
                    call orthogonalize_vector(trotv(1,1), lv, 5, zspace(1,i))
            end if
            ! Check norm
            if (sqrt(ddot(lv, zspace(1,i), 1, zspace(1,i), 1)) .lt. 1d-8) then
                    print *, "WARNING: Small norm: ", i, ": ",&
                            sqrt(ddot(lv, zspace(1,i), 1, zspace(1,i), 1))
            end if
            call normalize_vector(zspace(1,i), lv)
            i=i+1
    end do
    return
  end subroutine build_zspace

  !-------------------------------------------------------------------
  ! build_zvector: build z vector and its orthogonal compliment for
  ! dispacement not in gh plane
  subroutine build_zvector (zcomp, natoms, zspace, zv, ozv)
    implicit none
    ! ..input scalars..
    integer, intent(in) :: natoms
    ! ..input arrays..
    integer, dimension(3*natoms-6), intent(in) :: zcomp
    double precision, dimension(3*natoms,(3*natoms-6)), intent(in) :: zspace
    ! ..output arrays..
    double precision, dimension(3*natoms), intent(out) :: zv, ozv
    ! ..local scalars..
    integer :: i
    double precision :: norm
    ! ..external functions..
    double precision, external :: dnrm2
    zv = 0d0
    ozv = 0d0
    do i = 1, 3*natoms-6
            if (zcomp(i) .ne. 0) then
                    zv = zv + zspace(:,zcomp(i))
            else
                    ozv = ozv + zspace(:,zcomp(i))
            end if
    end do
    norm = dnrm2(3*natoms, zv, 1)
    zv = zv / norm
    norm = dnrm2(3*natoms, ozv, 1)
    ozv = ozv / norm
    return
  end subroutine build_zvector

  !-------------------------------------------------------------------
  ! check_atom_dist: check distance between two atoms
  function check_atom_dist(atom1, atom2, geom, na3) result(dist)
    implicit none

    ! ..input scalars..
    integer, intent(in) :: atom1, atom2, na3
    ! ..input arrays..
    double precision, dimension(na3), intent(in) :: geom
    ! ..external functions..
    double precision, external :: ddot
    ! ..local arrays..
    double precision, dimension(3) :: tmp

    double precision :: dist

    tmp = geom(atom1*3-2:atom1*3)-geom(atom2*3-2:atom2*3)
    dist = sqrt(ddot(3,tmp,1,tmp,1))
  end function check_atom_dist

  !-------------------------------------------------------------------
  function compute_ftheta(xgeom,gv,hv,zv,lv,rval,atheta,nst,st1,st2,x0,y0,&
          z0, cg, dcg, h, e, frho, fzv, fozv, ozv, ofij,init, plvl) result (val)
    use ioutil, only: print_colgeom2
    implicit none

    ! ..input scalars..
    ! lv = 3*natoms
    ! nst = number of states
    ! st1 = lower state
    ! st2 = upper state
    ! plvl = print level
    integer, intent(in) :: lv, nst, st1, st2, plvl
    double precision, intent(in) :: rval, atheta, x0, y0, z0
    logical, intent(in) :: init
    
    ! ..input arrays..
    double precision, dimension(lv) :: xgeom, gv, hv, zv, ozv, ofij
    double precision, dimension(lv,nst,nst), intent(out) :: cg, dcg
    double precision, dimension(nst,nst), intent(out) :: h
    double precision, dimension(nst), intent(out) :: e

    ! ..output values..
    double precision, intent(out) :: frho, fozv, fzv
    
    ! ..local scalars..
    double precision :: x, y, fx, fy, fth, dpfo
    double precision :: xt, yt, zt
    
    double precision :: val
    double precision, external :: ddot, dnrm2
    
    ! ..local arrays..
    double precision, dimension(lv) :: geom, fij, nfij

    cg = 0d0
    dcg = 0d0
    h = 0d0
    e = 0d0
    frho = 0d0
    
    ! Get x and y of new coordinates
    call polar2xy(rval, atheta, x, y)
    x = x + x0
    y = y + y0

    call gh2geom(x, y, gv, hv, lv, geom)
    call addzvec(z0, zv, lv, geom)
    ! test circulation
    xt = ddot(lv, gv, 1, geom, 1)
    yt = ddot(lv, hv, 1, geom, 1)
    zt = ddot(lv, zv, 1, geom, 1)
    print "('Location = (',f8.5,',',f8.5,',',f8.5,')')", xt, yt, zt
    if ((abs(xt - x) .gt. 1d-5) .or. (abs(yt - y) .gt. 1d-5)) then
            print "('*** WARNING: (x,y) differs from (g.G,h.G):')"
            print "(5x,'(xt, x, yt, t) =',4f15.8)", xt, x, yt, y
    end if
    if (abs(zt - z0) .gt. 1d-5) then
            print "('*** WARNING: z differs from z.G:')"
            print "(5x,'(zt, z0) = ',4f15.8)", zt, z0
    end if
    geom = geom + xgeom

    if (check_atom_dist(7, 13, geom, lv) .ge. 3.5) then
            print *, "WARNING: O-H Distance > 3.5"
    end if
    
    call EvaluateSurfgen(geom, e, cg, h, dcg)
    if (plvl .gt. 4) then
            print *, " -- current geometry --"
            call print_colgeom2(geom,(lv/3),atom_names, atom_numbers, &
                    atom_weights)
            print *, " ----------------------"
    end if
    fij = cg(1:lv, st1, st2)
    nfij = (-1d0)*fij
    if (ddot(lv, nfij, 1, ofij, 1) .gt. ddot(lv, fij, 1, ofij, 1)) then
            fij = nfij
    end if
!    dpfo = ddot(lv, fij, 1, ofij, 1)
!    fij = fij * dpfo / abs(dpfo)
        
    call get_ghcoord(fij, gv, hv, lv, fx, fy)
    call xy2polar_fij(fx, fy, rval, atheta, frho, fth)

    ! get z and oz components of fij
    fzv = ddot(lv, fij, 1, zv, 1)
    fozv = ddot(lv, fij, 1, ozv, 1)
    val = fth
    if (abs(fth) .gt. 1d-12) ofij = fij
  end function compute_ftheta
  !-------------------------------------------------------------------
  ! compute_intgrl: compute integral via trapezoidal rule.
  !  Sf(x)dx = (F(b) + f(a) / 2)*(b - a)
  function compute_intgrl(xgeom,gv,hv,lv,rval,ath,bth,nst,st1,st2,a,b,x0,y0,&
          fijphase, fthphase, fthpred) result(val)
    implicit none

    ! ..input scalars..
    ! lv = 3*natoms
    ! nst = number of states
    ! st1 = lower state
    ! st2 = upper state
    ! a = distance along x-axis (ellipse)
    ! b = distance along y-axis (ellipse)
    ! x0 = x coorindate of origin (ellipse)
    ! y0 = y coordinate of orogin (ellipse)
    ! rval = a
    ! ath = Theta_a
    ! bth = Theta_b
    integer, intent(in) :: lv, nst, st1, st2
    double precision, intent(in) :: rval, ath, bth, a, b, x0, y0
    
    ! ..input arrays..
    double precision, dimension(lv), intent(in) :: xgeom, gv, hv

    ! ..input/output scalars..
    ! fijphase = phase of previous dTh fij of point a
    ! fthphase = phase of previous fij_theta of point a
    ! fthpred = prediction of next fij_theta of point b + 1
    double precision, intent(inout) :: fijphase, fthphase
    double precision, intent(inout) :: fthpred
    
    ! ..output scalars..
    double precision :: val

    ! ..local scalars..
    ! adrdth = dr/dTheta at a
    ! bdrdth = dr/dTheta at b
    double precision :: x, y
    double precision :: ax, ay, bx, by
    double precision :: fax, fay, fbx, fby
    double precision :: fath, fbth, farho, fbrho
    double precision :: arho, brho, theta
    double precision :: atheta, btheta
    double precision :: adrdth, bdrdth
    double precision :: rho
    
    ! ..local arrays..
    double precision, dimension(lv) :: ageom, bgeom
    double precision, dimension(lv) :: afij, bfij
    double precision, dimension(lv, nst, nst) :: dcg, cg
    double precision, dimension(nst, nst) :: h
    double precision, dimension(nst) :: e

    rho = rval
    
    ! Get x and y of new coordinates
    call polar2xy(rho, ath, ax, ay)
    ax = ax + x0
    ay = ay + y0
    call polar2xy(rho, bth, bx, by)
    bx = bx + x0
    by = by + y0


    
    ! Get new geometries
    call gh2geom(ax, ay, gv, hv, lv, ageom)
    ageom = ageom + xgeom
    call gh2geom(bx, by, gv, hv, lv, bgeom)
    bgeom = bgeom + xgeom

    ! Evaluate surfgen at points a and b
    call EvaluateSurfgen(ageom, e, cg, h, dcg)
    afij = cg(1:lv,st1,st2)
    call EvaluateSurfgen(bgeom, e, cg, h, dcg)
    bfij = cg(1:lv,st1,st2)

    ! Get x and y and then rho and theta for points a and b.
    call get_ghcoord(afij, gv, hv, lv, fax, fay)
    call get_ghcoord(bfij, gv, hv, lv, fbx, fby)
    call xy2polar_fij(fax, fay, rho, ath, farho, fath)
    call xy2polar_fij(fbx, fby, rho, bth, fbrho, fbth)

    print *, "Nascent fath = ", fath
    fath = fath * fthphase
    if (abs(fbth - fthpred) .gt. abs(fbth*(-1d0) - fthpred)) then
            print *, "abs(fbth - fthpred) = ", abs(fbth - fthpred)
            print *, "abs(fbth*-1 - fthpred) = ", abs(fbth*(-1d0) - fthpred)
            fthphase = -1d0
            fbth = fbth * fthphase
    else
            fthphase = 1d0
            fbth = fbth * fthphase
    end if
    fthpred = linear_extrap(fath, fbth, (bth - ath))
    
    call print_intgrl_info(rho, ath, bth, ax, ay, bx, by, farho, fath,&
            fbrho, fbth, fthpred, e, h, st1, st2, nst, bfij, lv)

    ! Compute integral. We only need angular components.
    val = (fbth + fath) / 2
    val = val * (bth - ath)

    return
  end function compute_intgrl  

  !-------------------------------------------------------------------
  ! continuity_check: enforce continuity of fij(theta) and fij(rho)
  subroutine continuity_check(fij_pred, ft, fr, fz, foz)
    implicit none
    double precision, dimension(4), intent(in) :: fij_pred
    double precision, intent(inout) :: ft, fr, fz, foz
    double precision, dimension(4) :: diffs
    integer :: max_el
    print *, "Predictions = ", fij_pred
    print *, "Actual      = ", ft, fr, fz, foz
    diffs(1) = abs(fij_pred(1) - ft) - abs(fij_pred(1) - (-1d0)*ft)
    diffs(2) = abs(fij_pred(2) - fr) - abs(fij_pred(2) - (-1d0)*fr)
    diffs(3) = abs(fij_pred(3) - fz) - abs(fij_pred(3) - (-1d0)*fz)
    diffs(4) = abs(fij_pred(4) - foz)- abs(fij_pred(4) - (-1d0)*foz)
    max_el = find_maximum_abs(diffs, 4)
    if (diffs(max_el) .gt. 0) then
            print "(a,i3)", "Switching element: ", max_el
            ft = (-1d0) * ft
            fr = (-1d0) * fr
            fz = (-1d0) * fz
            foz= (-1d0) * foz
    end if
    return
  end subroutine continuity_check
  !-------------------------------------------------------------------
  ! evaluate_mab: evaluate molecular aharanov-bohm effect.
  subroutine evaluate_mab(xgeom, x0, y0, z0, rho, natoms, gphase, v1, v2, v3,&
          maxint, remtranrot, startangle, circdir, rotations, plvl)
    implicit none
    ! ..input scalars..
    ! natoms = number of atoms
    ! x0, y0 = coordinates in IAC of loop center
    ! rho    = radius of loop
    ! gphase = phase of g vector
    ! zvec = vector in branching space to displace origin by
    ! zdisp = value of z displacement
    ! startangle = starting position of loop
    ! circdir = circulation direction. >=0 : counterclockwise
    !                                   <0 : clockwise
    ! rotations = number of rotations
    ! plvl = print level
    integer, intent(in) :: natoms, v1, v2, circdir, rotations
    double precision, intent(in) :: x0, y0, z0, rho, gphase, startangle
    integer, intent(in) :: maxint
    integer, intent(in) :: plvl

    ! ..input logicals..
    logical, intent(in) :: remtranrot
    
    ! ..input arrays..
    ! xgeom = reference geometry
    integer, dimension(3*natoms-6), intent(in) :: v3
    double precision, dimension(3*natoms), intent(in) :: xgeom

    ! ..local scalars..
    ! natm = number of atoms from getInfo()
    ! nst  = number of states
    ! st1  = lowest state of intersection
    ! st2  = highest state of intersection
    ! ios  = input/output status and error flag
    integer :: natm, nst
    integer :: st1, st2
    integer :: ios
    
    ! ..local arrays..
    ! cg = adiabatic gradients
    ! dcg= diabatic gradients
    ! h  = diabatic matrix
    ! e  = energy values
    ! gv = gvector at xgeom
    ! hv = hvector at xgeom
    ! zv = z vector
    ! ozv= orthogonal compliment to z vector in zspace
    double precision, dimension(:,:,:), allocatable :: cg, dcg
    double precision, dimension(:,:),   allocatable :: h
    double precision, dimension(:),     allocatable :: e
    double precision, dimension(:),     allocatable :: gv, hv, sv
    double precision, dimension(:,:),   allocatable :: zspace
    double precision, dimension(:),     allocatable :: zv
    double precision, dimension(:),   allocatable :: ozv
    
    ! ..external functions..
    double precision, external :: ddot
    
    ! Initialize potential and evaluate surface at xgeom. Check to ensure
    ! two states are degenerate.
    call initPotential()
    call getInfo(natm, nst)
    print *, "natm = ", natm
    print *, "nst  = ", nst
    if (natm .ne. natoms) then
            print *, "*** Error: natm .ne. natoms! ***"
            print *, "natm = ", natm, " natoms = ", natoms
            return
    end if
    allocate(cg(3*natoms, nst, nst))
    allocate(dcg(3*natoms, nst, nst))
    allocate(e(nst))
    allocate(h(nst, nst))
    call EvaluateSurfgen(xgeom, e, cg, h, dcg)
    call find_degeneracy(e, nst, st1, st2, ios)
    print "('Adiabatic Energies: ', 5f18.2)", e(1:nst)*AU2CM1
    if (ios .eq. 0) then
            print *, "*** No degeneracy found! ***"
            print *, "Continue? (1=yes,0=no)"
            read *, ios
            if (ios .eq. 0) stop "Stopping."
            print *, "Enter states to evaluate MAB between."
            print *, "st1: "
            read *, st1
            print *, "st2: "
            read *, st2
    end if

    ! Build, orthogonalize and normalize g and h vectors.
    allocate(gv(3*natoms))
    allocate(hv(3*natoms))
    allocate(sv(3*natoms))
    call build_gvector(cg, natoms, nst, st1, st2, gv)
    call build_hvector(cg, natoms, nst, st1, st2, hv, e)
    call build_svector(cg, natoms, nst, st1, st2, sv)
    call orthog_ghvecs(gv, hv, natoms)
    call normalize_vector(gv, natoms*3)
    call normalize_vector(hv, natoms*3)
    print "('Sx =',f8.5)", ddot(natoms*3, sv, 1, gv, 1)
    print "('Sy =',f8.5)", ddot(natoms*3, sv, 1, hv, 1)
    gv = gv * gphase
    print *, "GVECTOR:"
    print "(3f15.10)", gv
    print *, "HVECTOR:"
    print "(3f15.10)", hv

    ! Build 3*NATOMS-6 = N orthogonal vectors (G and H are vectors 1 and 2)
    allocate(zspace(natoms*3,(natoms*3)-6))
    allocate(ozv(natoms*3))
    call build_zspace(natoms, xgeom, natoms*3, gv, hv, zspace, remtranrot)
    if (plvl .gt. 2) then
            print *, "ZSPACE:"
            print "(3f15.10)", zspace
    end if
    allocate(zv(natoms*3))
    call build_zvector(v3, natoms, zspace, zv, ozv)
    if (plvl .gt. 0) then
            print *, "ZVECTOR:"
            print "(3f15.10)", zv
    end if
    
    ! Perform MAB analysis.
    call mab_analysis(xgeom, x0, y0, z0, rho, zspace(1,v1), zspace(1,v2), &
            zv, natoms*3, st1, st2, nst, maxint, ozv, startangle, &
            circdir, rotations, plvl)
    
    return
  end subroutine evaluate_mab

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

  !-------------------------------------------------------------------
  ! find_maximum_abs: return index of element with largest absolute value
  function find_maximum_abs(vec, lv) result(index)
    implicit none
    integer, intent(in) :: lv
    double precision, dimension(lv), intent(in) :: vec
    integer :: index
    integer :: i

    index = 1
    do i = 1, lv
            if (abs(vec(i)) .gt. abs(vec(index))) then
                    index = i
            end if
    end do
    return
  end function find_maximum_abs

  !-------------------------------------------------------------------
  ! get_ellipse_rho: compute the rho of an ellipse given an angle
  ! 1 = (x - x0)^2/a^2 + (y - y0)^2/b^2
  !   = (r*cosT - x0)^2/a^2 + (r*sinT - y0)^2/b^2
  ! Self-consistently solve the follwing:
  !  r^2 = (a^2  * b^2) / (b^2(cosT - x0/r)^2 + a^2(sinT - y0/r)^2)
  subroutine get_ellipse_rho(rho, angle, x0, y0, a, b)
    implicit none

    ! ..input scalars..
    double precision, intent(in) :: angle, x0, y0, a, b

    ! ..output scalars..
    double precision, intent(out) :: rho

    ! ..local scalars..
    double precision :: rhs
    integer :: i
    double precision, parameter :: rtol = 1d-8
    
    
    return
  end subroutine get_ellipse_rho
  
  !-------------------------------------------------------------------
  ! get_ghcoord: get the g and h coordinates of a geometry
  ! x = g, y = h
  subroutine get_ghcoord(geom, gv, hv, lv, x, y)
    implicit none

    ! ..input scalars..
    integer, intent(in) :: lv

    ! ..input arrays..
    ! geom = geometry 
    double precision, dimension(lv), intent(in) :: geom, gv, hv

    ! ..output scalars..
    double precision, intent(out) :: x, y

    ! ..external functions..
    double precision, external :: ddot
    
    x = ddot(lv, geom, 1, gv, 1)
    y = ddot(lv, geom, 1, hv, 1)
    return
  end subroutine get_ghcoord

  !-------------------------------------------------------------------
  ! gh2geom: get geometry from g and h coordinates
  subroutine gh2geom(x, y, gv, hv, lv, geom)
    implicit none

    ! ..input scalars..
    integer, intent(in) :: lv
    double precision, intent(in) :: x, y

    ! ..input arrays..
    double precision, dimension(lv), intent(in) :: gv, hv

    ! ..output arrays..
    double precision, dimension(lv), intent(out) :: geom

    geom = gv*x + hv*y
    return
  end subroutine gh2geom

  !-------------------------------------------------------------------
  ! linear_extrap: linearly extrapolate between two points to predict
  ! the next value at dx
  function linear_extrap(f1, f2, dx) result(val)
    implicit none

    ! ..input scalars..
    double precision, intent(in) :: f1, f2, dx

    ! ..output scalars..
    double precision :: val

    val = f2 + (f2 - f1)
    return
  end function linear_extrap

  !-------------------------------------------------------------------
  ! mab_analysis: evaluates line integral S(fij)dR around conical intersection
  subroutine mab_analysis(xgeom, x0, y0, z0, rho, gv, hv, zv, lv, st1, st2, nst,&
          maxint, ozv, sa, cdir, nrot, plvl)
    implicit none

    ! ..input scalars..
    ! lv = 3 * natoms
    ! st1 = lower state
    ! st2 = upper state
    ! nst = total states
    ! x0, y0 = origin of loop in IAC
    ! rho = rdius of loop in IAC
    ! sa = starting angle
    ! cdir = circulation direction
    ! nrot = number of rotations
    integer, intent(in) :: lv, st1, st2, nst, cdir, nrot
    double precision, intent(in) :: x0, y0, z0, rho, sa
    integer, intent(in) :: maxint, plvl

    ! ..input arrays..
    ! gv = normalized g vector
    ! hv = normalized h vector
    ! zv = displacement vector
    ! ozv= orthononal compliment in Zspace of displacement vector
    ! xgeom = intersection geometry
    double precision, dimension(lv), intent(in) :: gv, hv, zv, ozv
    double precision, dimension(lv), intent(in) :: xgeom
    
    ! ..local scalars..
    ! rval = value of rho
    ! theta= rotation angle
    ! dth = dTheta
    ! x = geom.g
    ! y = geom.h
    ! fth = fij in terms of theta
    ! intgrl = integral value
    ! ointgrl = previous integral value
    ! inttol = integral tolerance
    ! maxiter = integration loop iterations
    ! btheta = theta of new geometry
    ! atheta = theta of old geometry
    ! fijphase = coupling phase
    ! fthphase = theta component of coupling phase
    ! fthpred = prediction of next fth value
    double precision :: rval
    double precision :: dth
    double precision :: intgrl, ointgrl
    double precision :: atheta, xtheta, btheta
    double precision, parameter :: inttol = 1d-12
    integer, parameter :: maxiter = 100
    double precision :: fijphase, fthphase, fthpred, frhopred
    double precision :: fija, fijx, fijb, frhox, frhoa, frhob
    double precision :: fza, fzx, fzb, foza, fozx, fozb
    ! ..local arrays..
    ! geom0 = origin geometry
    ! fijx = f_theta for values of simpson's rule formula
    ! escr = energy scratch array
    ! cgscr = gradient scratch array
    ! dcgscr = diabatic gradient scratch array
    ! hscr = diabatic hamiltonian scratch array
    double precision, dimension(lv) :: geom0
    double precision, dimension(lv) :: prev_fij
    double precision, dimension(lv, nst, nst) :: cgscr, dcgscr
    double precision, dimension(nst, nst) :: hscr
    double precision, dimension(nst) :: escr
    double precision, dimension(4) :: fij_phase, fij_pred
    integer :: i, j

    double precision, external :: ddot,dnrm2

    rval = rho

    ! Get origin (of loop) geometry
    call gh2geom(x0, y0, gv, hv, lv, geom0)
    call addzvec(z0, zv, lv, geom0)
    geom0 = geom0 + xgeom
    print "(A,f8.6,A,f8.6,A,f8.6,A)", "Origin: (",x0,",",y0,",",z0,")"
    print "('Number of rotations: ',i5)", nrot
    ! Integration loop. Compute each value for fTheta and Theta. Check
    ! continuity of fth values.
    intgrl = 0d0
    dth = 2*MATH_PI/maxint
    if (cdir .lt. 0) dth = dth*(-1d0)
    atheta = sa*2*MATH_PI/360d0
    prev_fij = hv
    fija = compute_ftheta(xgeom, gv, hv, zv, lv, rval, atheta, nst, &
            st1, st2, x0, y0, z0, cgscr, dcgscr, hscr, escr, frhoa, &
            fza, foza, ozv, prev_fij,.false., plvl)
    print *, "fija(0) = ", fija
    do i = 1, maxint*nrot
            print *, "============================================"
            print *, "Int = ", i
            btheta = atheta + dth
            xtheta = (atheta + btheta) / 2d0
            fijx = compute_ftheta(xgeom, gv, hv, zv, lv, rval, xtheta, nst, &
                    st1, st2, x0, y0, z0, cgscr, dcgscr, hscr, escr, frhox, &
                    fzx, fozx, ozv, prev_fij,.false., plvl)
            print *, "fijx = ", fijx
            print *, "PREV FIJ: "
            print "(3f15.8)", prev_fij
            print *, "RAW FIJ: "
            print "(3f15.8)", cgscr(:,st1,st2)
            fijb = compute_ftheta(xgeom, gv, hv, zv, lv, rval, btheta, nst, &
                    st1, st2, x0, y0, z0, cgscr, dcgscr, hscr, escr, frhob, &
                    fzb, fozb, ozv, prev_fij,.false., plvl)
            intgrl = intgrl + &
                    simpsons_rule(fija, fijx, fijb, atheta, xtheta, btheta)
            atheta = btheta
            fija = fijb
            frhoa= frhob

            print *, "atheta = ", atheta
            print *, "fij.g = ", ddot(lv, prev_fij, 1, gv, 1)
            print *, "fij.h = ", ddot(lv, prev_fij, 1, hv, 1)
            print *, "|fij| = ", dnrm2(lv, prev_fij, 1)
            if (st2 .ne. nst) then
                    print *, "|fij(3,4)| =", dnrm2(lv, cgscr(1:lv,st2,st2+1), 1)
                    print *, "fij(3,4).g =", ddot(lv,cgscr(1:lv,st2,st2+1),1,&
                            gv,1)
                    print *, "fij(3,4).h =", ddot(lv,cgscr(1:lv,st2,st2+1),1,&
                            hv,1)
            end if
            print *, "prev_fij.cgscr =", ddot(lv, cgscr(1:lv,st1,st2),1,prev_fij,1)
            print *, ""
            print *, "PREV FIJ: "
            print "(3f15.8)", prev_fij
            print *, "RAW FIJ: "
            print "(3f15.8)", cgscr(:,st1,st2)
            print *, "GVECTOR: "
            print "(3f15.8)", gv
            print *, ""
            print *, "fija = ", fija
            print *, "frho = ", frhob
            print *, "intgrl = ", intgrl
            print "(a,5(f12.8,','))", "energies = ", escr(1:nst)
            print *, "============================================"

    end do
    
    print *, "Integral = ", intgrl
    
    
    return
  end subroutine mab_analysis
  !-------------------------------------------------------------------
  ! normalize_vector: normalize a vector
  subroutine normalize_vector(v, lv)
    implicit none
    ! ..input scalars..
    integer, intent(in) :: lv

    ! ..input/output arrays..
    double precision, dimension(lv), intent(inout) :: v

    ! ..local scalars..
    double precision :: vnrm
    
    ! ..external functions..
    double precision, external :: dnrm2

    vnrm = dnrm2(lv, v, 1)
    v = v / vnrm
    return
  end subroutine normalize_vector

  !-------------------------------------------------------------------
  ! orthog_ghvecs: orthogonalize g and h vectors
  subroutine orthog_ghvecs(g, h, na)
    implicit none
    ! ..input scalars..
    ! na = number of atoms
    integer, intent(in) :: na

    ! ..input/output arrays..
    double precision, dimension(3*na), intent(inout) :: g, h

    ! ..local scalars..
    double precision :: gh, gnrm, hnrm, tin, theta

    ! ..local arrays..
    double precision, dimension(3*na) :: scrg, scrh

    ! ..external functions..
    double precision, external :: ddot, dnrm2

    ! Get theta, and rotate vectors. Check dot product
    scrg = g
    scrh = h
    gh = ddot(3*na, g, 1, h, 1)
    print "(A,e15.8)", "Before rotation: G.H = ", gh
    gnrm = dnrm2(3*na, g, 1)
    hnrm = dnrm2(3*na, h, 1)
    tin = 2 * gh /(hnrm**2 - gnrm**2)
    theta = datan(tin)/4d0
    g = dcos(theta*2d0)*scrg - dsin(theta*2d0)*scrh
    h = dsin(theta*2d0)*scrg + dcos(theta*2d0)*scrh
    gh = ddot(3*na, g, 1, h, 1)
    print "(A,f8.5)",  "Theta = ", theta
    print "(A,e15.8)", "After rotation:  G.H = ", gh
    print "(A,e15.8)", "Norm of G: ", dnrm2(3*na, g, 1)
    print "(A,e15.8)", "Norm of H: ", dnrm2(3*na, h, 1)
    return
  end subroutine orthog_ghvecs

  !-------------------------------------------------------------------
  ! orthogonalize_vector: orthogonalize a vector to a space of normalized
  ! vectors.
  subroutine orthogonalize_vector(space, lv, dim, vec)
    implicit none

    ! ..input scalars..
    integer, intent(in) :: lv, dim

    ! ..input arrays..
    double precision, dimension(lv, dim), intent(in) :: space

    ! ..input/output arrays..
    double precision, dimension(lv), intent(inout) :: vec

    ! ..local scalars..
    integer :: i, j
    double precision :: ovrlp

    ! ..external functions..
    double precision, external :: ddot
    
    do i = 1, dim
            ovrlp = ddot(lv, space(1,i), 1, vec, 1)
            do j = 1, lv
                    vec(j) = vec(j) - ovrlp*space(j,i)
            end do
    end do
    return
  end subroutine orthogonalize_vector
  
  !-------------------------------------------------------------------
  ! print_intgrl_info: print integration informatin
  subroutine print_intgrl_info(rho, ath, bth, ax, ay, bx, by, fijrho, fijth, &
          fijrhob, fijthb, fpred, e, hd, st1, st2, nst, fij, nat3)
    implicit none

    ! ..input scalars..
    double precision, intent(in) :: rho, ath, bth, ax, ay, bx, by, fijrho
    double precision, intent(in) :: fpred, fijth, fijrhob, fijthb
    integer, intent(in) :: nst, st1, st2, nat3
    double precision, dimension(nst), intent(in) :: e
    double precision, dimension(nst,nst), intent(in) :: hd
    double precision, dimension(nat3), intent(in) :: fij

    double precision, dimension(nst) :: de
    double precision, external :: dnrm2
    integer :: i

    do i = 1, nst
            de(i) = hd(i,i)
    end do
    
    print "(15('-'))"
    print 10, "A", ax, ay, rho, ath
    print 10, "B", bx, by, rho, bth
    print "('Fij_T(a)= ', f18.12,', Fij_R(a)= ', f18.12,',')", fijth, fijrho
    print "('Fij_T(b)= ', f18.12,', Fij_R(b)= ', f18.12,',')", fijthb, fijrhob
    print "('Prediction = ', f15.8)", fpred
    print "('Adiabatic Fij = ',f20.12)", dnrm2(nat3, fij, 1)
    print "('Diabatic  Energies = ', 5(f15.2,','))", de(:)*AU2CM1
    print "('Adiabatic Energies = ', 5(f15.2,','))", e(:)*AU2CM1
    print "(15('-'))"

    return
10  format(1x,'Point ',A2,'(',f15.8,',',f15.8,') r =',f15.8,' O=',f15.8)
  end subroutine print_intgrl_info
  
  !-------------------------------------------------------------------
  ! polar2xyz: convert polar coordiantes to xy coordinates
  subroutine polar2xy(rho, theta, x, y)
    implicit none

    ! ..input scalars..
    double precision, intent(in) :: rho, theta

    ! ..output scalars..
    double precision, intent(out) :: x, y

    x = rho * dcos(theta)
    y = rho * dsin(theta)
    return
  end subroutine polar2xy

  !-------------------------------------------------------------------
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

  end function simpsons_rule
  !-------------------------------------------------------------------
  ! xy2polar: convert xy coordinates to polar coordiantes
  subroutine xy2polar(x, y, rho, theta)
    implicit none

    ! ..input scalars..
    double precision, intent(in) :: x, y

    ! ..output scalars..
    double precision, intent(out) :: rho, theta

    rho = sqrt(x**2 + y**2)
    theta = datan(y/x)
    return
  end subroutine xy2polar

  !-------------------------------------------------------------------
  ! xy2polar_fij: convert fij(x,y) to fij(r,theta)
  subroutine xy2polar_fij(x, y, r, th, rho, theta)
    implicit none

    ! ..input scalars..
    ! x = fij_x
    ! y = fij_y
    ! r = rho of geometry (in g,h plane)
    ! th = theta of geometry (in g,h plane)
    double precision, intent(in) :: x, y, r, th

    ! ..output scalars..
    ! fij_rho
    ! fij_theta
    double precision, intent(out) :: rho, theta
    
    rho = dcos(th)*x + dsin(th)*y
    theta = r*(dcos(th)*y - dsin(th)*x)
    return
  end subroutine xy2polar_fij
end module mabutil
