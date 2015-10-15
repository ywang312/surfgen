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
!-------------------
! opttools contains subroutines used to analyse the geometry, to calculate
! numerical hessian matrix, to search for critical points using Newton-Raphson
! procedure, and to calculate harmonic frequencies using hessian matrix. 
module opttools 
 implicit none
 contains

!---calculate harmonic frequencies from hessian matrix
subroutine getFreq(natoms,masses,hess,w)
  implicit none
  integer,intent(in)          :: natoms
  double precision,intent(in) :: masses(natoms),hess(3*natoms,3*natoms)
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
  call DSYEVD('N','U',3*natoms,hmw,3*natoms,w,tmp,-1,itmp,-1,INFO)
  if(info/=0)print *,"DSYEVD allocation investigation failed.  info=",info
  LWORK = int(tmp(1))
  LIWORK= itmp(1)
  allocate(WORK(LWORK))
  allocate(IWORK(LIWORK))
  
  call DSYEVD('N','U',3*natoms,hmw,3*natoms,w,WORK,LWORK,IWORK,LIWORK,INFO)
  if(info/=0)print *,"DSYEVD failed.  info=",info
  do i=1,3*natoms
    if(w(i)<0)then
      w(i) = -sqrt(-w(i))*au2cm1
    else 
      w(i) = sqrt(w(i))*au2cm1
    end if
  end do
end subroutine getFreq


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
subroutine findmin(natoms,nstate,cgeom,isurf,maxiter,shift,Etol,Stol)
  implicit none
  integer, intent(in)                                 ::  natoms,isurf,maxiter,nstate
  double precision,dimension(3*natoms),intent(inout)  ::  cgeom
  double precision,intent(in)                         ::  shift,Etol,Stol

  real*8    ::  h(nstate,nstate),cg(3*natoms,nstate,nstate),dcg(3*natoms,nstate,nstate),e(nstate)
  double precision,dimension(3*natoms)  ::  grad, b1, b2, w
  double precision,dimension(3*natoms,3*natoms)  :: hess,hinv
  double precision,dimension(:),allocatable :: WORK
  integer,dimension(:),allocatable :: IWORK
  integer           :: LIWORK, LWORK, itmp(1),INFO  
  integer  ::  iter  , i
  double precision            :: nrmG, nrmD, tmp(1)
  double precision, external  :: dnrm2
  double precision,  parameter  :: amu2au=1.822888484514D3,au2cm1=219474.6305d0
  double precision, parameter   :: MAXD = 1D-1

  ! initialize work spaces
  call DSYEVD('V','U',3*natoms,hess,3*natoms,w,tmp,-1,itmp,-1,INFO)
  if(info/=0)print *,"DSYEVD allocation failed.  info=",info
  LWORK = int(tmp(1))
  LIWORK= itmp(1)
  allocate(WORK(LWORK))
  allocate(IWORK(LIWORK))

  print "(A,I4,A,I4,A)","Searching for minimum on surface ",isurf," in ",maxiter," iterations."
  print "(A)","  Convergence Tolerances"
  print "(A,E10.2,A,E10.2)","  Energy Gradient: ",Etol,"   Displacement:",Stol
  do iter=1,maxiter
     call EvaluateSurfgen(cgeom,e,cg,h,dcg)
     grad = cg(:,isurf,isurf)
     nrmG=dnrm2(3*natoms,grad,1)
     call calcHess(natoms,cgeom,nstate,isurf,1D-4,hess,.true.)
     hinv = hess
     ! invert the hessian
     call DSYEVD('V','U',3*natoms,hess,3*natoms,w,WORK,LWORK,IWORK,LIWORK,INFO)
     if(info/=0)print *,"DSYEVD failed.  info=",info
     ! hess. w. hess^T = hess_old
     ! so x= hess_old^-1.b = hess. w^-1. hess^T. b
     ! b1= hess^T.g       =>
     call DGEMV('T',3*natoms,3*natoms,1d0,hess,3*natoms,grad,1,0d0,b1,1)
     ! b1' = w^-1*b1
     do i=1,3*natoms
      if(abs(w(i))<shift)then
         b1(i)=dble(0)
      else!
         b1(i)=b1(i)/w(i)
      end if!
     end do!
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
       return
     end if
  end do 
1000 format("Iteration ",I4,": E=",F20.4,", |Grd|=",E12.5,", |Disp|=",E12.5)
end subroutine findmin
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
  character*300                 ::  geomfile,str


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

  print *,"Reading input from input file "//trim(geomfile)
  call readColGeom(geomfile,1,natm,aname,anum,cgeom,masses)

  print "(/,A)","-------------- Initial Geometry ------------"
  call analysegeom(natm,cgeom,aname,anum,masses,1d8,.true.)
  print "(/,A)","----------- Geometry Optimizations ---------"
  call findmin(natm,nst,cgeom,isurf,100,1d-3,1d-8 ,1d-7)
  print "(/,A)","--------------- Final Geometry -------------"
  call analysegeom(natm,cgeom,aname,anum,masses,1d8,.true.)
  print "(/,A)","------------ Harmonic Frequencies ----------"
  call calcHess(natm,cgeom,nst,isurf,1D-4,hess,.true.,skip)
  call writeHess(hess,3*natm)
  call getFreq(natm,masses,hess,w)
  do i=1,3*natm
    print "(I5,F12.2)",i,w(i)
  end do

! deallocate arrays
  deallocate(masses)
  deallocate(anum)
  deallocate(aname)
  deallocate(hess)
  deallocate(w)
  deallocate(cgeom)
  deallocate(skip)

end


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
