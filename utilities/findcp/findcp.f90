! findcp
!-------------------
! Critical point searching program for potential generated by surfgen
!-------------------
! This is part of the surfgen program package.  To compile this program,
! first generate surfgen library with `make libs` in surfgen main directory,
! then run the `install.sh` script in the program's directory.
!-------------------
! (C) 2013 Yarkony Group, Johns Hopkins University
!-------------------
! Feb 2013        Xiaolei Zhu        Created for NH3 and phenol surfaces
! Sep 2013        Xiaolei Zhu        Updated to run for general sytems
! Nov 2013        cthree-40          Updated for use of input file.
! Nov 2013        cthree-40          Updated for generation of molden.freq file
!-------------------
! opttools contains subroutines used to analyse the geometry, to calculate
! numerical hessian matrix, to search for critical points using Newton-Raphson
! procedure, and to calculate harmonic frequencies using hessian matrix. 
module opttools 
 implicit none
 contains

!---calculate harmonic frequencies from hessian matrix
subroutine getFreq(natoms,masses,hess,w,cg,anm,pl,mol)
  implicit none
  integer,intent(in)          :: natoms, pl
  double precision,intent(in) :: masses(natoms),hess(3*natoms,3*natoms)
  double precision,intent(in) :: cg(3*natoms)
  logical,intent(in)          :: mol
  character*3, dimension(natoms), intent(in) :: anm
  double precision,intent(out):: w(3*natoms)

  double precision,  parameter  :: amu2au=1.822888484514D3,au2cm1=219474.6305d0
  integer  :: i,j
  double precision  :: sqrm,  hmw(3*natoms,3*natoms), tmp(1)
  double precision,dimension(:),allocatable :: WORK
  integer,dimension(:),allocatable :: IWORK
  integer           :: LIWORK, LWORK, itmp(1),INFO
  ! convert hessian into mass weighed coordinates
  hmw = hess/amu2au
  do i=1,natoms
    sqrm = sqrt(masses(i))
    do j=i*3-2,i*3
      hmw(j,:)=hmw(j,:)/sqrm
      hmw(:,j)=hmw(:,j)/sqrm
    end do
  end do 
  !calculate eigenvalues of hmw
  call DSYEVD('V','U',3*natoms,hmw,3*natoms,w,tmp,-1,itmp,-1,INFO)
  if(info/=0)print *,"DSYEVD allocation investigation failed.  info=",info
  LWORK = int(tmp(1))
  LIWORK= itmp(1)
  allocate(WORK(LWORK))
  allocate(IWORK(LIWORK))
  
  ! if print level is greater than 0 we want to print the modes.
  if (pl == 0) then
        call DSYEVD('N','U',3*natoms,hmw,3*natoms,w,WORK,LWORK,IWORK,LIWORK,INFO)
  else
        call DSYEVD('V','U',3*natoms,hmw,3*natoms,w,WORK,LWORK,IWORK,LIWORK,INFO)
  endif
  if(info/=0)print *,"DSYEVD failed.  info=",info
  
  
  
  do i=1,3*natoms
    if(w(i)<0)then
      w(i) = -sqrt(-w(i))*au2cm1
    else 
      w(i) = sqrt(w(i))*au2cm1
    end if
  end do
  
  ! print modes
  if (pl > 0) then
        print *," Modes:"
        do i=1,3*natoms
                print "(2x,'Mode ',i5,5x,'Frequency: ',f12.2)", i, w(i)
                do j=1, natoms
                        print "(2x,3f15.8)", &
                        hmw((j-1)*3 + 1,i), &
                        hmw((j-1)*3 + 2,i), & 
                        hmw((j-1)*3 + 3,i)
                end do
        end do
  end if

  ! print molden output?
  if (mol) then
          call gen_molden_file(cg, hmw, w, natoms, anm)
  end if
  return
end subroutine getFreq

! gen_molden_file: generate molden frequency file
subroutine gen_molden_file(cg, evec, eval, natoms, anames)
        implicit none
        integer, intent(in) :: natoms
        character*3,dimension(natoms) :: anames
        real*8, dimension(3*natoms), intent(in) :: cg
        real*8, dimension(3*natoms,3*natoms), intent(in) :: evec
        real*8, dimension(3*natoms), intent(in) :: eval
        
        integer :: mu=11
        character*25 :: mn

        integer :: ios, i, j

        write(mn,"(a)") "molden.freq"
        mn=trim(mn)

        open(file=mn,unit=mu,action='write',status='unknown',iostat=ios)
        if (ios .ne. 0) then
                print "(A)", "Could not open molden file. No file generated."
                return
        end if 

        ! print header
        write(mu,"(1x,A)") "-- > start of molden input"
        write(mu,"(1x,A)") "[Molden Format]"
        ! print frequencies
        write(mu,"(1x,A)") "[FREQ]"
        do i=1,3*natoms
                write(mu,"(f10.2)") eval(i)
        end do
        ! print geometry
        write(mu,"(1x,A)") "[FR-COORD]"
        do i=1, natoms
                write(mu,"(1x,a3,3f13.6)") trim(anames(i)), &
                        cg((i-1)*3+1),cg((i-1)*3+2),cg((i-1)*3+3)
        end do
        ! print modes
        write(mu,"(1x,A)") "[FR-NORM-COORD]"
        do i=1, 3*natoms
                write(mu,"(1x,'vibration',i24)") i
                do j=1,natoms
                        write(mu,"(3f13.5)") evec((j-1)*3+1,i), &
                                evec((j-1)*3+2,i), &
                                evec((j-1)*3+3,i)
                end do
        end do
        write(mu,"(1x,A)") "--> end of molden input"
end subroutine gen_molden_file

!---calculate hessian at a certain geometry
subroutine calcHess(natoms,cgeom,nstate,istate,stepsize,hessian,centerd,skip)
  implicit none
  integer, intent(in)          :: natoms, nstate,istate
  logical, intent(in),optional :: skip(natoms)
  logical, intent(in)          :: centerd !whether to do centered difference or only backward difference
  double precision,intent(in)  :: stepsize
  double precision,intent(in)  :: cgeom(3*natoms)
  double precision,intent(out) :: hessian(3*natoms,3*natoms)
  double precision   ::  mdif
  
  integer   ::   i,  j
  logical   ::   skipdisp(natoms*3)
  double precision  ::  dispgeom(3*natoms), dgrd(3*natoms),gref(3*natoms)
  real*8    ::  h(nstate,nstate),cg(3*natoms,nstate,nstate),dcg(3*natoms,nstate,nstate),e(nstate)
  skipdisp=.false.
  if(present(skip))then
    do i=1,natoms
      if(skip(i))skipdisp(i*3-2:i*3)=.true.
    end do
  end if
  ! to perform backward difference, gradient at references is needed
  if(.not.centerd)then
    call EvaluateSurfgen(cgeom,e,cg,h,dcg)
    gref = cg(:,istate,istate)
  end if
  hessian = 0d0
  do i=1,3*natoms
    if(skipdisp(i))cycle
    dispgeom=cgeom
    dispgeom(i)=dispgeom(i) - stepsize
    call EvaluateSurfgen(dispgeom,e,cg,h,dcg)
    dgrd  =-cg(:,istate,istate)
    if(centerd)then ! centered difference
      dispgeom=cgeom
      dispgeom(i)=dispgeom(i) + stepsize
      call EvaluateSurfgen(dispgeom,e,cg,h,dcg)
      dgrd = dgrd+cg(:,istate,istate)
      hessian(i,:)= dgrd/2/stepsize
    else !backward difference
      dgrd = dgrd+gref
      hessian(i,:)= dgrd/stepsize
    end if
  end do!o=1,3*natoms
  do i=1,3*natoms
    if(skipdisp(i))hessian(:,i)=0d0
  end do
  mdif = maxval(abs(hessian-transpose(hessian)))
  if(mdif>1d-5)print *,"maximum hermiticity breaking : ",mdif
  hessian = (hessian+transpose(hessian))/2
end subroutine calcHess


!----search for minimum on adiabatic surfaces 
subroutine findmin(natoms,nstate,cgeom,isurf,maxiter,shift,Etol,Stol,gscale,hess_disp,MAXD, &
   masses, sadd_srch, converge_test)
  implicit none
  integer, intent(in)                                 ::  natoms,isurf,maxiter,nstate
  double precision,dimension(3*natoms),intent(inout)  ::  cgeom
  integer, intent(inout)                              :: converge_test
  double precision,intent(in)                         ::  shift,Etol,Stol, gscale, hess_disp
  double precision, intent(in)                        :: MAXD !Max displacement size
  character(len=1),intent(in)                         :: sadd_srch
  
  real*8    ::  h(nstate,nstate),cg(3*natoms,nstate,nstate),dcg(3*natoms,nstate,nstate),e(nstate)
  double precision,dimension(3*natoms)  ::  grad, b1, b2, w
  double precision,dimension(3*natoms,3*natoms)  :: hess,hinv
  double precision,dimension(natoms)             :: masses
  double precision,dimension(:),allocatable :: WORK
  integer,dimension(:),allocatable :: IWORK
  integer           :: LIWORK, LWORK, itmp(1),INFO  
  integer  ::  iter  , i, eig_start
  double precision            :: nrmG, nrmD, tmp(1)
  double precision, external  :: dnrm2
  double precision,  parameter  :: amu2au=1.822888484514D3,au2cm1=219474.6305d0
  double precision,dimension(:),allocatable  :: rem_modes
  
  character*3,dimension(natoms) :: ctmp

  ! Derived type for normal mode removal in eigenvalue decomposition procedure.
  type mode_removed
      double precision  :: freq
      integer           :: removed 
  end type mode_removed
  ! use this type for normal mode removal list
  type(mode_removed),dimension(:),allocatable :: rem_mode
  
  ! Set convergence flag
  converge_test=0

  ! initialize work spaces
  call DSYEVD('V','U',3*natoms,hess,3*natoms,w,tmp,-1,itmp,-1,INFO)
  if(info/=0)print *,"DSYEVD allocation failed.  info=",info
  LWORK = int(tmp(1))
  LIWORK= itmp(1)
  allocate(WORK(LWORK))
  allocate(IWORK(LIWORK))
  
  ! Allocate mode matrix
  ! This array will allow users to see what modes are being removed in the eigenvalue decomposition
  !   procedure.
  allocate(rem_modes(3*natoms))
  allocate(rem_mode(3*natoms))

  print "(A,I4,A,I4,A)","Searching for minimum on surface ",isurf," in ",maxiter," iterations."
  print "(A)","  Convergence Tolerances"
  print "(A,E10.2,A,E10.2)","  Energy Gradient: ",Etol,"   Displacement:",Stol
  do iter=1,maxiter
     call EvaluateSurfgen(cgeom,e,cg,h,dcg)
     grad = cg(:,isurf,isurf)
     nrmG=dnrm2(3*natoms,grad,1)
     call calcHess(natoms,cgeom,nstate,isurf,hess_disp,hess,.true.)
     if(iter.eq.1 .or. iter.eq.maxiter)then
        call writeHess(hess,3*natoms)
        call getFreq(natoms,masses,hess,w,cgeom,ctmp,0,.false.)
        do i=1,3*natoms
          print "(I5,F12.2)",i,w(i)
        end do
        do i=1,nstate
            print "(2x,A,I2,A,F12.8)", " Energy of state ", i,"= ", e(i)
        end do
     end if
     hinv = hess
     ! invert the hessian
     ! Call DSYEVD
     call DSYEVD('V','U',3*natoms,hess,3*natoms,w,WORK,LWORK,IWORK,LIWORK,INFO)
     ! hess now contains the orthonormal eigenvectors of hess(old)
     ! w contains the eigenvalues from lowest to highest
     if(info/=0)print *,"DSYEVD failed.  info=",info
     ! [Old Hessian] = [hess][w][hess]^T
     ! Thus, 
     !    x = [Old Hessian]^-1[b] = ([hess]w)^-1[hess]^T[b]
     !
     ! [b1]= [hess]^T[g]       =>
     ! Call DGEMV. perform [hess]^T[grad]=b1
     ! First, scale gradient
     grad=grad*gscale
     call DGEMV('T',3*natoms,3*natoms,1d0,hess,3*natoms,grad,1,0d0,b1,1)
     ! b1' = w^-1*b1
     ! Check if saddle point search
     if( sadd_srch .EQ. 'Y' ) then
       eig_start=2      ! First eigenvalue should be large and negative, so we skip it
       b1(1)=b1(1)/w(1)
     else
       eig_start=1      ! Otherwise, we continue as normal
     end if
     !
     rem_modes=0
     do i=eig_start,3*natoms
      if(abs(w(i))<shift)then
         b1(i)=b1(i)
         rem_mode(i)%freq=w(i)
         rem_mode(i)%removed=1
      else!
         b1(i)=b1(i)/w(i)
      end if!
     end do!
     ! Print out to what modes are being removed
     do i=1, 3*natoms
       if ( rem_mode(i)%removed .eq. 1 ) then      ! If mode was removed, print info
         write(*,1001) i, rem_mode(i)%freq
       end if
     end do
     ! b2 = hess.b1'
     call DGEMV('N',3*natoms,3*natoms,1d0,hess,3*natoms,b1,1,0d0,b2,1)
     ! cgeom' = cgeom - H^-1.g
     nrmD=dnrm2(3*natoms,b2,1)
     if(nrmD>maxD)then
       b2=b2/nrmD*maxD
       nrmD=maxD
     end if
     print 1000,iter,e(isurf)*au2cm1,nrmG,nrmD
     cgeom = cgeom - b2
     if(nrmG<Etol.and.nrmD<Stol)then
       print *,"Optimization converged"
       print "(A,10F20.4)","Final energy of all states : ",e*au2cm1
       converge_test=1
       return
     end if
  end do 
1000 format("Iteration ",I4,": E=",F20.4,", |Grd|=",E12.5,", |Disp|=",E12.5)
1001 format("Modes Removed: ",I4,F14.8)
end subroutine findmin
!==================================================================================
!>analysegeom2
! This subroutine prints the last geometry out to the file new.geom
!----------------------------------------------------------------------------------
subroutine analysegeom2(natoms,geom,aname,anum,masses,gflname)

  implicit none
  integer, intent(in)          ::  natoms
  character*3,intent(in)       ::  aname(natoms)
  character(len=300),intent(IN)      ::  gflname
  double precision,intent(in)  ::  anum(natoms),masses(natoms)
  double precision,intent(in)  ::  geom(3,natoms)
  double precision, parameter  ::  bohr2ang=0.529177249d0
  integer   ::  i, ios
  
  ! Open new file
  open(unit=9,file=gflname,status="unknown",position="rewind",iostat=ios)
  if (ios .ne. 0 ) then ! If file cannot be opened
      print "(1x,A)", "Could not open geom file. (analysegeom2)"
  end if
  do i=1,natoms
     write(unit=9, fmt="(x,a3,1x,f4.1,3F14.8,F14.8)") aname(i),anum(i),geom(:,i),masses(i)
  end do
  close(unit=9)
  return
end subroutine analysegeom2
!===================================================================================
end module opttools 

! main program
program findcp

 use opttools 
 implicit none

  integer :: natm , nst, ngeoms
  integer :: i,j, isurf, ios
  double precision,allocatable  :: anum(:),masses(:),  hess(:,:), w(:)
  character*3,allocatable       :: aname(:)
  double precision,allocatable  :: cgeom(:)
  logical,allocatable           :: skip(:)
  double precision,external     :: dnrm2
  character*300                 ::  geomfile,str,inputfile
  double precision,allocatable  :: dchess(:)
  integer                       :: un_infile
  integer                       :: converge_test
  
  integer               :: niter, printlvl
  double precision      :: egrad_tol, shift, disp_tol, grad_scale, hess_disp, maxdisp
  character*1           :: sadd_search
  logical               :: check_inputfl, molden, ant_output
  character*300         :: new_geomfl, old_geomfl
! Namelist input for inputfile
  namelist /cpsearch/   niter, egrad_tol, shift, disp_tol, grad_scale, hess_disp, &
      maxdisp, sadd_search, printlvl, molden, ant_output
! If not read in, here are the default values:
  niter=100                         ! Max iterations
  egrad_tol=1d-9                    ! Energy gradient tolerance for convergence
  shift=1d-5                        ! Value of shift
  disp_tol=1d-5                     ! Step size tolerance for convergence
  un_infile=40                      ! Unit number of input file
  grad_scale=1.0                    ! Scaling of gradient
  hess_disp=1d-5                    ! Displacement for hessian calculation
  maxdisp =1d-1                     ! Size of maximum displacement
  sadd_search='N'                   ! Saddle point specific searching not implemented yet
  old_geomfl="old.geom"             ! Old geometry file
  new_geomfl="new.geom"             ! New geometry file
  printlvl=0                        ! Print level
  molden=.false.                    ! Generate molden output
  ant_output=.false.                ! Print final geometry in ANT format

  print *," ***************************************** "
  print *," *    findcp.x                           * "
  print *," ***************************************** "
  print *," Critical point search on Surfgen surface" 
  call initPotential()
  call getInfo(natm,nst)

! allocate arrays

  print "(A,I6)","Number of Atoms:  ",natm
  print "(A,I6)","Number of States: ",nst 

  allocate(masses(natm))
  allocate(anum(natm))
  allocate(aname(natm))
  allocate(hess(natm*3,natm*3))
  allocate(w(3*natm))
  allocate(cgeom(3*natm))
  allocate(skip(natm))
  allocate(dchess(3*natm))    !Dimension of diagonal correction to hessian
  dchess=0d0
  skip = .false.

! process arguments
! synposis:    findcp.x geom states
! Default values:
! geom        geom
! states      1
  call get_command_argument(number=1,value=geomfile,status=ios)
  if(ios<-1)  &     ! input argument larger than container
      stop "Filename truncated.  Please use less than 300 characters"
  if(ios>0) then
    print *,"No filename suppied.  Using default."
    write(geomfile,"(A)"), "geom"
    isurf = 1
  else
    call get_command_argument(number=2,value=str,status=ios)
    if(ios/=0)then
      print *,"Cannot get surface number from command line options.  Using default."
      isurf = 1
    else
      read(str,*,iostat=ios)isurf
      if(ios/=0)then
        print *,"Cannot get surface number from command line options.  Using default."
        isurf = 1
      end if
    end if
  end if

! Reading namelist input
  call get_command_argument(number=3,value=inputfile,status=ios)
  if(ios/=0)then      ! If no filename is provided, use default
    print *, "Cannot find filename from command line options.  Using default."
    inputfile="findcp.in"
  endif
! Open input file
  inquire(file=inputfile,exist=check_inputfl)
  if (check_inputfl) then
      open(unit=un_infile,file=inputfile,iostat=ios)
            read(unit=un_infile,nml=cpsearch)
      close(unit=un_infile)
  else
        print *, "No input file found. Continuing with default values."
  end if

  print *,"Reading input from input file "//trim(geomfile)
  call readColGeom(geomfile,1,natm,aname,anum,cgeom,masses)

  print "(/,A)","-------------- Initial Geometry ------------"
  call analysegeom(natm,cgeom,aname,anum,masses,2d0,.true.)
  print "(/,A)"," Printing original geometry "
  call analysegeom2(natm,cgeom,aname,anum,masses,old_geomfl)
  print "(/,A)","----------- Geometry Optimizations ---------"
  converge_test=0
  call findmin(natm,nst,cgeom,isurf,niter,shift,egrad_tol ,disp_tol,grad_scale,hess_disp,&
    maxdisp, masses, sadd_search, converge_test)
  print "(/,A)","--------------- Final Geometry -------------"
  call analysegeom(natm,cgeom,aname,anum,masses,2d0,.true.)
  print "(/,A)","------------ Harmonic Frequencies ----------"
  call calcHess(natm,cgeom,nst,isurf,hess_disp,hess,.true.,skip)
  call writeHess(hess,3*natm)
  call getFreq(natm,masses,hess,w,cgeom,aname,printlvl,molden)
  do i=1,3*natm
    print "(I5,F12.2)",i,w(i)
  end do
  call compute_zeropt(natm, w, 3*natm)
 
  ! Print final geometry
  call analysegeom2(natm,cgeom,aname,anum,masses,new_geomfl)

  if (ant_output) then
      call print_ant_output(natm, cgeom, aname, anum, masses)
  end if
  
! deallocate arrays
  deallocate(masses)
  deallocate(anum)
  deallocate(aname)
  deallocate(hess)
  deallocate(w)
  deallocate(cgeom)
  deallocate(skip)

end program findcp

!*
! compute_zeropt: compute and print zero point energy.
!*
subroutine compute_zeropt(na, freqs, nvibs)
  implicit none
  integer, intent(in) :: na, nvibs
  double precision, dimension(nvibs), intent(in) :: freqs
  integer :: i
  double precision :: zp
  zp = 0d0
  do i = 1, nvibs
          if (freqs(i) .gt. 0d0) then
                  zp = zp + freqs(i)
          end if
  end do
  zp = zp / 2d0
  print "(A,f10.2)", "Zero point energy = ", zp
  return
end subroutine compute_zeropt

!*
! print_ant_output: print geometry in ANT initial geometry format
!*
subroutine print_ant_output(na, g, aname, anum, masses)
    implicit none
    integer, intent(in) :: na
    real*8, dimension(3, na),   intent(in) :: g
    real*8, dimension(na),      intent(in) :: anum, masses
    character(3), dimension(na),intent(in) :: aname
    real*8, parameter :: BOHR2ANG = 0.529177
    integer :: i
    print *, ""
    print *, "Final geometry (ANT format):"
    do i = 1, na
        print 100, trim(aname(i)), masses(i), g(1,i) * BOHR2ANG, &
          g(2,i) * BOHR2ANG, g(3,i) * BOHR2ANG
    end do
    return
100 format(a2,4f14.8)
end subroutine print_ant_output

!------ write hessian file to disk, columbus format
subroutine writeHess(hess,nvibs)
  implicit none
  integer, intent(in)         :: nvibs
  double precision,intent(in) :: hess(nvibs,nvibs)
  integer,parameter     :: hessfl= 32175
  integer  :: ios,i

  open (unit=hessfl,file="hessian",status='replace',position='rewind',access='sequential',&
            action='write',iostat=ios)
  if(ios/=0)then
    print *,"Cannot open hessian file for write."
    return
  end if
  do i=1,nvibs
     write(hessfl,"(8F13.6)")hess(i,:)
  end do
  close(hessfl)
end subroutine
