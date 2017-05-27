c
c***************************************************************
c
      subroutine readin
c
       implicit none
c  ##  parameter & common block section
c
c   ma is the max. number of atoms
      integer ma
      parameter (ma=200)
c   lq is the max. number of cartesisans
      integer lq
      parameter (lq=3*ma)
c   mq is the max. number of internal coordinates,
      integer mq, maxcard
      parameter (mq=lq-6)
c
      character*10 calctype
      integer          intgradout,intforceout,rgfcoord,fixc(mq),
     &                 geomch
      common /control/ intgradout,intforceout,rgfcoord,fixc,geomch,
     &                 calctype
c
c  # formats
      character*30    geomfm, gradfm, intforcfm, intgeomfm, intgradfm
      common /format/ geomfm, gradfm, intforcfm, intgeomfm, intgradfm
c
c  # unit numbers
      integer       inp,       nlist,    ngeom,    ngrad,      nintc,
     &              nintgeom,  nintforc, ngeomnew, nintgeomch, nintgrad,
     &              nbmat,     ngconv
      common /tape/ inp,       nlist,    ngeom,    ngrad,      nintc,
     &              nintgeom,  nintforc, ngeomnew, nintgeomch, nintgrad,
     &              nbmat,     ngconv
c
      integer           linonly,    linallow,    maxiter, forceinau
      real*8            maxstep
      common/transform/ linonly,    linallow,    maxstep,  maxiter,
     &                  forceinau
c
c  # file names
      character*40     c2ils,       c2iinp,    geomfl,    gradfl,
     &                 intcfl,      intgeomfl, intforcfl, geomnewfl,
     &                 intgeomchfl, intgradfl, bmatfl,    gconvfl
      common /flnames/ c2ils,       c2iinp,    geomfl,    gradfl,
     &                 intcfl,      intgeomfl, intforcfl, geomnewfl,
     &                 intgeomchfl, intgradfl, bmatfl,    gconvfl
c congrad thresholds
      integer maxiter_congrad
      real*8 small_congrad
      common/congrad_thr/maxiter_congrad, small_congrad
c
c     # bummer error types.
      integer   wrnerr,  nfterr,  faterr
      parameter(wrnerr=0,nfterr=1,faterr=2)
c
c
      integer ierr
c
      namelist /input/
     &    calctype,    gradfl,       gradfm,        intgradout,
     &    intgradfl,   intgradfm,    intforceout,   intforcfl,
     &    intforcfm,   geomfl,       geomfm,        intgeomfl,
     &    intgeomfm,   intgeomchfl,  geomnewfl,     linonly,
     &    linallow,    maxstep,      maxiter,       rgfcoord,
     &    fixc,        geomch,       forceinau,     maxiter_congrad,
     &    small_congrad

c
c----------------------------------------------------
c
c  ##  set default values
c
      calctype='xx'
      intgradout=0
      intforceout=0
      maxstep=0.1d+00
      maxiter=50
      linonly=0
      linallow=0
      rgfcoord=0
      geomch = 0
      forceinau = 0
      maxiter_congrad = 1000
      small_congrad = 1d-12
      call izero_wr(mq,fixc,1)
c
c  ##  readin namelist
c
      open(unit=inp, file=c2iinp,status = 'old')
c
      write(nlist,'(1x,a)' )
     & '======== echo of the cart2int input ========'
      call echoin( inp, nlist, ierr )
      if ( ierr .ne. 0 ) then
         call bummer('readin: from echoin, ierr=',ierr,faterr)
      endif
c
      rewind(inp)
      read( inp, input, end=3 )
3     continue
      close(unit=inp)
c
c  ## consistency check for the input values
c
      if ((intgradout.lt.0).or.(intgradout.gt.1)) call bummer
     & ('INPUT ERROR: intgradout not allowed=',intgradout,faterr)
c
      if ((intforceout.lt.0).or.(intforceout.gt.1)) call bummer
     & ('INPUT ERROR: intforceout not allowed=',intforceout,faterr)
c
      if (geomfl.eq.geomnewfl) call bummer
     & ('INPUT ERROR: geomfl and geomnew cannot be the same',0,faterr)
c
      if ((calctype.eq.'cart2int').and.(rgfcoord.gt.0)) call bummer
     & ('INPUT ERROR: rgfcoord<>0 allowed for calctype=int2cart',
     & 0,faterr)
c
      if (calctype.eq.'xx') call bummer
     & ('INPUT ERROR: calctype not set',0,faterr)
c
      if (rgfcoord.lt.0) call bummer
     & ('INPUT ERROR: rgfcoord must be positive',0,faterr)
c
      if ((calctype.ne.'cart2int').and.(calctype.ne.'int2cart'))
     & call bummer
     & ('INPUT ERROR: wrong calctype, allowed: cart2int or int2cart',0,
     & faterr)
c
      if (maxiter.le.0) call bummer
     & ('maxiter value not allowed, maxiter=',maxiter,faterr)
c
      if (maxstep.le.0.0d+00) call bummer
     & ('nonpositive maxstep value not allowed',0,faterr)
c---     & ('maxstep value not allowed, maxstep=',maxstep,faterr)
c
      if ((linonly.lt.0).or.(linonly.gt.1)) call bummer
     & ('linonly value not allowed, linonly=',linonly,faterr)
c
      if ((linallow.lt.0).or.(linallow.gt.1)) call bummer
     & ('linallow value not allowed, linallow=',linallow,faterr)
c
      if (fixc(1).ne.0) write(nlist,*)'set fixc'
c
c  #  now read geometry
      call readgeom
c
c  #  now read gradient
      call readgrad
c
      return
      end
c
c***************************************************************
c
      subroutine absmax(n,a,i,xmax)
c  this routine returns the highest absolute value in the
c array a, from a(1) to a(n)
       implicit none 
       integer i,j,n
      real*8 a(n), xmax
c
      i=1
      xmax=abs(a(1))
      do 100 j=2,n
        if (xmax.lt.abs(a(j))) then
           xmax=abs(a(j))
           i=j
        end if
 100  continue
      return
      end
c
c***************************************************************
c
      real*8 function arc1 (x,y)
      implicit real*8 (a-h,o-z)
      real*8     small
      parameter( small=1d-11 )
      real*8     pi
      parameter( pi=3.14159265358979323846264338328d0 )
      real*8     two
      parameter( two=2d0 )
      real*8     pis2
      parameter( pis2 = pi / two )
c
      if (abs(x) .lt. small) go to 10
      s = atan(y/x)
      if (x. lt. 0.0) s = s + pi
      arc1 = s
      return
   10 arc1 = pis2
      return
      end
c
c***************************************************************
c
      subroutine bbt(nq,nek,b,ibcontr,x,y,t)
c
c  This subroutine computes y= B(Bt)x. t is just storage
c  parameters: INPUT
c              nq=number of internal coordinates
c              nek=3*na, number of Cartesians
c              b(54,nq): contains the non-zero elements of B
c              ibcontr(20,nq): coding info for B
c              x(nq): input vector
c              OUTPUT
c              y: B(Bt)x
      implicit real*8 (a-h,o-z)
      dimension b(54,nq),ibcontr(20,nq),x(nq),y(nq),t(nek)
      call btxv(nq,nek,b,ibcontr,x,t)
      call bxv(nq,nek,b,ibcontr,t,y)
      end
c
c***************************************************************
c
      subroutine bdiag(nq,b,ibcontr,dm1)
c  parameters: INPUT
c              nq=number of internal coordinates
c              b(54,nq): contains the non-zero elements of B
c              ibcontr(20,nq): coding info for B
c              OUTPUT
c              dm1=diag(BBtranspose)**-1 (inverse)
      implicit real*8 (a-h,o-z)
      parameter(zero=0.0d0,one=1.0d0)
      dimension b(54,nq),ibcontr(20,nq),dm1(nq)
      do i=1,nq
         s=zero
         natom=ibcontr(2,i)
         k3=0
         do k=1,natom
            k3=k3+3
            s=s+b(k3-2,i)**2
            s=s+b(k3-1,i)**2
            s=s+b(k3,i)**2
         end do
         if (s .ne. zero) then
            dm1(i) = one / s
         else
            dm1(i) = zero
         endif
      end do
      return
      end
c
c***************************************************************
c
      subroutine bdiag1(nq,nek,b,ibcontr,dm1)
c  parameters: INPUT
c              nq=number of internal coordinates
c              nek=3*na, number of Cartesians
c              b(54,nq): contains the non-zero elements of B
c              ibcontr(20,nq): coding info for B
c              OUTPUT
c              dm1=diag(BtransposeB)**-1 (inverse)
      implicit real*8 (a-h,o-z)
      parameter(one=1.0d0,eps=1.0d-9)
      dimension b(54,nq),ibcontr(20,nq),dm1(nek)
      call wzero(nek,dm1,1)
      do i=1,nq
         natom=ibcontr(2,i)
         k3=0
         do k=1,natom
            k3=k3+3
            iatom=ibcontr(k+2,i)
            iat3=iatom*3
            dm1(iat3-2)=dm1(iat3-2)+b(k3-2,i)**2
            dm1(iat3-1)=dm1(iat3-1)+b(k3-1,i)**2
            dm1(iat3)=dm1(iat3)+b(k3,i)**2
         end do
      end do
      do k=1,nek
         if(abs(dm1(k)).gt.eps) then
            dm1(k)=one/dm1(k)
         else
            dm1(k)=eps
         end if
      end do
      end
c
c********************************************************************
C JAN 11, 1995 PP: The definition of INV6 is now 1000/r**6
C SEP 13 PP CHANGED THE LIMIT DETERMINANT IN GDIIS TO 1E-10;
C  PUT IN THE RESTRICTION ON THE DIIS COEFFICIENTS WHICH WERE COMM. OUT
C JUN 28 PP CHANGED IT SO THAT A CARD BEGINNING WITH A BLANK SPACE
C    IS NO MORE MISTAKEN FOR A CONTINUATION LINE IN MACHB
c
c
      subroutine bread(na,inp,nlist,nq,nprim,ibcode,ibcontr,wri)
c
c
c   this routine reads the b matrix information and converts it to the
c   integer arrays ibcode(6,nprim) and ibcontr(20,nq).
c   nprim is the number of primitive internal coordinates,
c   ibcode(1,ip) is the type: 1=stretching,
c                             2=invr (1/r),
c                             3=inv6 (1/r**6),
c                             4=bending,
c                             5=out of plane,
c                             6=torsion,
c                             7=linear1,
c                             8=linear2
c   ibcode(2,ip) is the coefficient of the coordinate in the composite
c    internal coordinate. It is the product of the overall scale factor
c    for the coordinate and the (normalized coefficiient. It
c    is expressed as an integer, in units of 10**(-8) (to avoid
c    defining a structure containing both real and integer elements)
c   ibcode(3..6,ip) are the atoms participating in the coordinate
c    for stretchings only ibcode(3,ip) and ibcode(4,ip) are used etc.
c   ibcontr(nprim+1) is an array describing the contraction:
c    in goes from ibcontr(1,k-1)+1 to ibcontr(1,k), except when k=1.
c    ibcontr(2,k) holds the number of different atoms appearing in
c    the contracted internal coordinate (max. 10), and ibcontr(3..20,k)
c    gives these atoms.

C                                   INPUT DATA
C        EACH ELEMENTARY VALENCE COORDINATE ON A SEPARATE CARD
C        COL. 1 'K' OR' ' (BLANK). IF 'K' A NEW COORDINATE BEGINS, IF
C        BLANK THEN THE COMPOSITE INTERNAL COORDIMNATE BEGUN EARLIER IS
C         CONTINUED. ANY OTHER CHARACTER TERMINATES THE INPUT
C        COLS. 2-9 SCALE FACTOR FOR THE TOTAL COORDINATE (ONLY IF THERE
C        IS K IN COLUMN 1. BLANK OR ZERO IS INTERPRETED AS 1.0
C        COLS. 21-24 COORDINATE TYPE STRE,INVR,BEND,OUT ,TORS,LIN1,LIN2
C        COLS. 31-40,41-50,51-60,61-70 PARTICIPATING ATOMS A,B,C,D
C         FORMAT 4F10.X
C         A AND B ARE GIVEN FOR 'STRE' ORDER ARBITRARY
C         A AND B ARE GIVEN FOR -INVR- ORDER ARBITRARY
C        A,B,C FOR BEND - A AND B ARE END ATOMS, C IS THE APEX ATOM A-C-
C        ATOM A OUT OF THE BCD PLANE - C IS THE CENTRAL ATOM - COORDINA
C        TE POSITIVE IF A IS DISPLACED TOWARD THE VECTOR PRODUCT DB*DC
C        TORSION A-B-C-D, POSITIVE AS IN THE WILSON-DECIUS-CROSS BOOK
C        NOTE THAT THE VALUE OF THE COORDINATE VARIES BETWEEN -PI/2 TO
C        3PI/2   N O T  BETWEEN -PI/2 TO +PI/2
C        LIN1 L  COLLINEAR BENDING A-B-C DISTORTED IN THE PLANR OF ABD
C        POSITIVE IF A MOVES TOWARD D
C         LIN2 LINEAR BENDING A-C-B DISTORTED PERPENDICULAR TO THE PLANE
C        ABD - POSITIVE IF A MOVES TOWARD THE VECTOR CROSS PRODUCT CD*CA
C        THE LINEAR BENDINGS ARE A-C-B, I. E. THE APEX ATOM IS THIRD
c
c
c  ## parameter & common block section
c
c     # bummer error types.
      integer   wrnerr,  nfterr,  faterr
      parameter(wrnerr=0,nfterr=1,faterr=2)
c
c   ma is the max. number of atoms
      integer ma
      parameter (ma=200)
c   lq is the max. number of cartesisans
      integer lq
      parameter (lq=3*ma)
c   mq is the max. number of internal coordinates,
      integer mq, maxcard
      parameter (mq=lq-6)
c
c  # number of types of internal coordinates
      integer ntyp
      parameter(ntyp=8)
c
c  # defines the internal coordinate type
      integer             intctype(mq)
      CHARACTER*4         tipus(ntyp)
      common /icoorclass/ intctype
      common /ccoorclass/ tipus
c
c  ##  character section
c
      CHARACTER*1 WEW(80),WORT(3),WE
      CHARACTER*4 BLANK,TYP,TLAST
      CHARACTER*10 TENCHA,LBLANK
c
      dimension ibcode(6,*),ibcontr(20,*),a(4),ia(4)
c
      EQUIVALENCE (KA,IA(1)), (KB,IA(2)), (KC,IA(3)), (KD,IA(4))
c
c  ##  datq section
c
      DATA WORT/'K',' ','C'/
      DATA BLANK/'    '/
      DATA LBLANK/'          '/
c
c  ##  logical section
c
      logical wri, already
c
      real*8     zero,     one,     big
      parameter( zero=0d0, one=1d0, big=1d8 )
c
c     # initialization of common array moved from data. 20-oct-01 -rls
      tipus(1) = 'STRE'
      tipus(2) = 'INVR'
      tipus(3) = 'INV6'
      tipus(4) = 'BEND'
      tipus(5) = 'OUT '
      tipus(6) = 'TORS'
      tipus(7) = 'LIN1'
      tipus(8) = 'LIN2'
C
c-----------------------------------------------------------------
c
      NCARD=0
c  this is the (contracted) coordinate counter
      I=0
      C1=zero
      CCX=zero
      TLAST=BLANK
      IF (WRI) WRITE (nlist,10)
   10 FORMAT (/,1X,34HDEFINITION OF INTERNAL COORDINATES,/)
   20 READ (INP,30,END=45) WEW
      WE=WEW(1)
   30 FORMAT (80A1)
      BACKSPACE INP
      DO 40 K=1,3
         IF (WE.EQ.WORT(K)) GO TO 60
   40 CONTINUE
   45 CONTINUE
      IF (I.EQ.0) THEN
        WRITE(nlist,*) ' NO INTERNAL COORDINATES WERE FOUND ON',INP
      END IF
      C1 = SQRT(one / C1) * CCX
c  put in the last value for ibcontr
      nq=i
      ibcontr(1,nq)=ncard
      if(i.eq.1) then
        ilow=1
      else
        ilow=ibcontr(1,i-1)+1
      end if
      ihigh=ibcontr(1,i)
      DO 50 K=ilow,ihigh
         ibcode(2,k)=ibcode(2,k)*c1
   50 CONTINUE
      GO TO 370
   60 CONTINUE
C     CHECK IF THE BLANK IS REALLY BLANK
      IF (WE.EQ.WORT(2)) THEN
         READ (INP,65) TENCHA
         BACKSPACE INP
 65      FORMAT(A10)
         IF (TENCHA.NE.LBLANK) GO TO 45
      END IF
      NCARD=NCARD+1
      READ (INP,70) WE,CC,C,TYP,A
   70 FORMAT (A1,F9.5,F10.4,A4,6X,4F10.4)
      if(typ.eq.blank.and.ncard.eq.1) then
        write(nlist,*) 'first coordinate must be defined'
        call bummer('first coordinate must be defined',0,faterr)
      end if
      IF (TYP.EQ.BLANK) TYP=TLAST
      TLAST=TYP
      IF (CC .EQ. zero) CC = one
      IF (WE.EQ.WORT(1))then
          CCC=CCX
          CCX=CC
          ibcontr(1,i)=ncard-1
      end if
      IF (C .EQ. zero) C = one
      IF (WE.EQ.WORT(2).OR.WE.EQ.WORT(3)) C1=C1+C**2
      IF (WE.EQ.WORT(1)) then
         IF (I.ne.0) then
            IF (WRI) WRITE (nlist,80)
80          FORMAT (1X)
            C1 = SQRT(one / C1) * CCC
            if(i.eq.1) then
               ilow=1
            else
               ilow=ibcontr(1,i-1)+1
            end if
            ihigh=ibcontr(1,i)
            DO 90 K=ilow,ihigh
               ibcode(2,k)=ibcode(2,k)*c1
90          CONTINUE
         end if
         I=I+1
         C1=C**2
      end if
  130 DO 140 K=1,4
         ia(k)=a(k)+0.01d0
         ibcode(k+2,ncard)=ia(k)
  140 CONTINUE
      DO 150 K=1,ntyp
         IF (TYP.EQ.TIPUS(K)) GO TO 170
  150 CONTINUE
      WRITE (nlist,160) I
  160 FORMAT(/' UNDEFINED INT. COORDINATE TYPE AT NO. ',I3/1x,40('*'))
      GO TO 380
  170 IF (WRI) WRITE (nlist,180) I,TYP,IA,C,CCX
  180 FORMAT (1X,I3,1H.,A8,4I10,F12.5,F12.6)
      IF (KA.LT.1.OR.KA.GT.NA.OR.KB.LT.1.OR.KB.GT.NA) GO TO 350
      IF (K.GT.3.AND.(KC.LT.1.OR.KC.GT.NA)) GO TO 350
      IF (K.GT.4.AND.(KD.LT.1.OR.KD.GT.NA)) GO TO 350
      ibcode(1,ncard) = k
      ibcode(2,ncard) = c * big
      go to 20
  350 continue
      WRITE (nlist,360) I
  360 FORMAT(/' ATOMS ERRONOUSLY DEFINED,COORDINATE NO. ',I3/1X,
     & 40('*'))
  370 continue
c     fill in incontr
      do 500 i=1,nq
        if(i.eq.1) then
          ilow=1
        else
          ilow=ibcontr(1,i-1)+1
        end if
        ihigh=ibcontr(1,i)
c zero out ibcontr
        do 460 kk=3,12
          ibcontr(kk,i)=0
 460    continue
c  number of different primitives
        mprim=0
        DO 490 k=ilow,ihigh
          do 480 l=3,6
            iatom=ibcode(l,k)
            if(iatom.eq.0) go to 480
            already=.false.
            do 470 kk=3,20
              if(iatom.eq.ibcontr(kk,i)) already=.true.
              if(ibcontr(kk+1,i).eq.0) go to 475
 470        continue
 475        if(.not.already) then
              mprim=mprim+1
              if(mprim.gt.18)
     &          call bummer('mprim .gt.',18,faterr)
              ibcontr(mprim+2,i)=iatom
              ibcontr(2,i)=mprim
            end if
 480      continue
 490    continue
 500  continue
      return
  380 call bummer('error found',0,faterr)
C
      END
c
c***************************************************************
c
      subroutine btb(nek,nq,b,ibcontr,x,y,t)
c   calculates y=bt*b*x
       implicit none 
       integer nek,nq,ibcontr
      real*8  x(nek),y(nek),t(nq),b
c
      call bxv(nq,nek,b,ibcontr,x,t)
      call btxv(nq,nek,b,ibcontr,t,y)
c
      return
      end
c
c***************************************************************
c
      subroutine btxv(nq,nek,b,ibcontr,v,w)
c
c  this routine multiplies a vector v by B transpose and puts
c  the result in w
c  parameters: INPUT
c              nq=number of internal coordinates
c              nek=3*na, number of Cartesians
c              b(54,nq): contains the non-zero elements of B
c              ibcontr(20,nq): coding info for B
c              v(nq): input vector
c              OUTPUT
c              w=B(transpose)*v
      implicit real*8 (a-h,o-z)
      dimension b(54,nq),v(nq),w(nek),ibcontr(20,nq)
      call wzero(nek,w,1)
      do i=1,nq
        natom=ibcontr(2,i)
        k3=0
        do k=1,natom
          k3=k3+3
          iatom=ibcontr(k+2,i)
          iat3=iatom*3
          w(iat3-2)=w(iat3-2)+b(k3-2,i)*v(i)
          w(iat3-1)=w(iat3-1)+b(k3-1,i)*v(i)
          w(iat3)=w(iat3)+b(k3,i)*v(i)
        end do
      end do
      return
      end
c
c***************************************************************
c
      subroutine bxv(nq,nek,b,ibcontr,w,v)
c  parameters: INPUT
c              nq=number of internal coordinates
c              nek=3*na, number of Cartesians
c              b(54,nq): contains the non-zero elements of B
c              ibcontr(20,nq): coding info for B
c              w(nek): input vector
c              OUTPUT
c              v=B*w
      implicit real*8 (a-h,o-z)
      parameter(zero=0.0d0)
      dimension b(54,nq),v(nq),w(nek),ibcontr(20,nq)
      do i=1,nq
        s=zero
        natom=ibcontr(2,i)
        k3=0
        do k=1,natom
          k3=k3+3
          iatom=ibcontr(k+2,i)
          iat3=iatom*3
          s=s+b(k3-2,i)*w(iat3-2)
          s=s+b(k3-1,i)*w(iat3-1)
          s=s+b(k3,i)*w(iat3)
        end do
        v(i)=s
      end do
      return
      end
C
c************************************************************
c
      subroutine congrad(
     & nq,   nek,   b,    ibcontr,   x,
     & y,    dm1,   c,    t,         p,
     & r,    z,     oper, condamp)
c
c   this routine solves iteratively the equation Ax=c
c   inek,b,ibcontr,t are strictly arguments for "oper"
c   A is the operation = oper can be bbt or btb;
c   bbt forms y= (BBtranspose)x; t is a temporary variable, nek long
c   Btb forms y= Btranspose*B*x ; t is a  temp. variable
c   Algorithm describe e.g. by A. Greenbaum, C. Li and H. Z. Chao,
c   "Practical Iterative Methods for Large-Scale Computations",
c   ed.  D.L. Boley, D.G. Truhlar, Y. Saad, R.E. Wyatt, L.E. Collins
c   North-Holland, Amsterdam 1989.
c
c   The righ-hand side is c, and the final solution is x.
c   dm1 will hold the inverse diagonal of the matrix A=BBt
c
c   For general implementation, replace oneiter and its arguments
c   and also bdiag. The latter forms and inverts diag(A)**-1
c   precon multiplies a vector with diag(A)**-1
c   y,dm1,t,p,r,z are just workspaces
c
c   Parameters: INPUT
c   nq: number of internal coordinates
c   nek: number of Cartesians (3*na)
c   b: B matrix
c   ibcontr: contraction info for the compactly stored B matrix
c   c: right hand side of the equation
c   dm1: inverse diagonal of the matrix (in this case BBt)
c      OUTPUT
c   x: solution of the equation
c   oper: the external subroutine
c   new parameters added by ph
c   condamp:  damping factor
c   nlist: output filnumber
c
       implicit none
c     ##  parameter & common section
c
c     # constants
      real*8          zero,  half,  one,   two, three,  four,  five,
     & ten,   ten6,  tenm8, pi,  acc
      common /number/ zero,  half,  one,   two, three,  four,  five,
     & ten,   ten6,  tenm8, pi,  acc
c
c     # unit numbers
      integer       inp,       nlist,    ngeom,    ngrad,      nintc,
     &              nintgeom,  nintforc, ngeomnew, nintgeomch, nintgrad,
     &              nbmat,     ngconv
      common /tape/ inp,       nlist,    ngeom,    ngrad,      nintc,
     &              nintgeom,  nintforc, ngeomnew, nintgeomch, nintgrad,
     &              nbmat,     ngconv
      integer maxiter_congrad
      real*8 small_congrad
      common/congrad_thr/maxiter_congrad, small_congrad
c
c     # bummer error types.
      integer   wrnerr,  nfterr,  faterr
      parameter(wrnerr=0,nfterr=1,faterr=2)
c
      real*8 onem
      parameter(onem=-1.0d0)
      integer maxiter
      real*8    small
c
c

c     ##  integer section
c
      integer i0, i1, iter, ii, imx
      integer nq, nek


      integer ibcontr(20,nq)
c
c     ##  real*8 section
c
      real*8 aa(2),a(2),absmx
      real*8 b(54,nq),bb
      real*8 c(nq),condamp
      real*8 dm1(nq), damp
      real*8 p(nq)
      real*8 r(nq)
      real*8 sum,s
      real*8 t(nek)
      real*8 x(nq)
      real*8 y(nq)
      real*8 z(nq)
c
c     ##  external section
c
      real*8   ddot_wr
      external ddot_wr
      external oper
c
c-------------------------------
c      maxiter = maxiter_congrad
c      small = small_congrad
      maxiter = 1000
      small = 1d-12
      one = 1d0
c     # set up circular counter
      i0=1
      i1=2
      damp=0d0
c     # set the first colums of x to D**-1*c   x0=DM*c
      call precon(c,dm1,x,nq)
c     # y=A*x0
      call oper(nq,nek,b,ibcontr,x,y,t)
c     # put in a little damping
      call daxpy_wr(nq,damp,x,1,y,1)
c     # r0=c-Ax0
      call dcopy_wr(nq,c,1,r,1)
      call daxpy_wr(nq,onem,y,1,r,1)
c     # z=DM*r0
      call precon(r,dm1,z,nq)
c     # p0=z
      call dcopy_wr(nq,z,1,p,1)
c     # aa0=z(t)r0
      sum=ddot_wr(nq,z,1,r,1)
      aa(1)=sum
c
      iter=0
100   continue
      iter=iter+1
c     # y=Ap(k-1)
      call oper(nq,nek,b,ibcontr,p,y,t)
c     # put in a little damping
      call daxpy_wr(nq,damp,p,1,y,1)
c     # a(k-1)=aa(k-1)/p(k-1)(t)y
      sum=ddot_wr(nq,p,1,y,1)
      a(i0)=aa(i0)/sum
c&ps
      if(sum .eq. zero) a(i0) = aa(i0)
c&ps
c     # x(k)=x(k-1)+a(k-1)p(k-1)
      call daxpy_wr(nq,a(i0),p,1,x,1)
c     # r(k)=r(k-1)-a(k-1)*y
      s=-a(i0)
      call daxpy_wr(nq,s,y,1,r,1)
c     # check for convergence
      call absmax(nq,r,imx,absmx)
      if((absmx .lt. small) .or. (iter .gt. maxiter)) go to 1000
c     # z=DM*r(k)
      call precon(r,dm1,z,nq)
c     # aa(k)=z(t)*r(k);   bb=aa(k)/aa(k-1)
      sum=ddot_wr(nq,z,1,r,1)
      aa(i1)=sum
      bb=aa(i1)/aa(i0)
c     # p(k)=z+bb*p(k-1)
      call dscal_wr(nq,bb,p,1)
      call daxpy_wr(nq,one,z,1,p,1)
c
      ii=i1
      i1=i0
      i0=ii
      go to 100
1000  continue
      write(nlist,
     & '(15x,''linear equation solved in '',i3,'' itererations'')')
     & iter
      write(nlist,'(15x,''lin. eq. sol. residuum= '',e13.8)') absmx
c
      if (iter .gt. maxiter) then
         call bummer(
     &    'error in congrad, no convergence, iter=', iter, faterr)
      endif
      return
      end
c
c*********************************************************************
c
      real*8 function darcos (x)
      implicit real*8 (a-h,o-z)
      real*8    zero,     one,     onem,      two,     small
      parameter(zero=0d0, one=1d0, onem=-one, two=2d0, small=1d-11 )
      real*8    pi
      parameter(pi=3.14159265358979323846264338328d0)
      real*8    pis2
      parameter(pis2=pi/two)
      if (x .ge. one) go to 10
      if (x .le. onem) go to 20
      x1 = sqrt(one - x**2)
      if (abs(x) .lt. small) go to 30
      s = atan( x1 / x )
      if (x .lt. zero) s = s + pi
      darcos = s
      return
10    darcos = zero
      return
20    darcos = pi
      return
30    darcos = pis2
      return
c
      end
c
c**********************************************************
c
      subroutine distort (
     & nek,     nq,    na,      qc,      qq,
     & qq1,     xa,    ibcode,  ibcontr, b,
     & shftpt,  xy)
c
c   this routine distorts the molecule along given internal
c   coordinates and gives the Cartesians as a result
c
c   INPUT
c   nek=3*na, (input) number of cartesians
c   na (input) number of nuclei
c   nq (input) number of internal coordinates
c   qc (input) new internal coordinates
c   qq (input) original internal coordinates
c   qq1 (workspace, contains the final coordinates on entry)
c   xa(3,na)(input) cartesians in angstrom units
c   ibcode(6,*): coding info for the B matrix
c   ibcontr(20,nq): contraction for the B matrix
c   b (workspace) b matrix
c   shftpt (input) see machb
c   OUTPUT
c   xy(3,na)  (output) distorted cartesians
c
       implicit none
c  ## parameter & common block section
c
c
c     # bummer error types.
      integer   wrnerr,  nfterr,  faterr
      parameter(wrnerr=0,nfterr=1,faterr=2)
c
c  # unit numbers
      integer       inp,       nlist,    ngeom,    ngrad,      nintc,
     &              nintgeom,  nintforc, ngeomnew, nintgeomch, nintgrad,
     &              nbmat,     ngconv
      common /tape/ inp,       nlist,    ngeom,    ngrad,      nintc,
     &              nintgeom,  nintforc, ngeomnew, nintgeomch, nintgrad,
     &              nbmat,     ngconv
c
c
      integer           linonly,    linallow,    maxiter, forceinau
      real*8            maxstep
      common/transform/ linonly,    linallow,    maxstep,  maxiter,
     &                  forceinau
c
c   ma is the max. number of atoms
      integer ma
      parameter (ma=200)
c   lq is the max. number of cartesisans
      integer lq
      parameter (lq=3*ma)
c   mq is the max. number of internal coordinates,
      integer mq, maxcard
      parameter (mq=lq-6)
c
c   this is the maximum internal coordinate change -otherwise
c   use scaling
      real*8  cartmax
      parameter (cartmax=0.25d0)
c
      real*8 onem
      parameter (onem=-1.0d+00)
c
      real*8 condamp
      common /options/condamp
c
      real*8 zero,half,one,two,three,four,five,ten,ten6,tenm8,
     &                pi,acc
      common /number/ zero,half,one,two,three,four,five,ten,ten6,tenm8,
     &                pi,acc
c
      integer lcx
      parameter (lcx=10 000)
c
      real*8 bl
      common bl(lcx)
c
c  ##  integer section
c
      integer ibcode(6,*), ix, iy, idm1, it, ip, ir, iz, ibcontr(20,*),
     &        iu, irep, i, imax
      integer lincalc
      integer na, nq, nek
c
c  ##  real*8 section
c
      real*8 b(54,nq)
      real*8 delmax
      real*8 qc(nq),qq(nq),qq1(nq),qcsave(mq)
      real*8 sc1, shftpt, scale
      real*8 toler
      real*8 xa(3,na),xy(3,na)
c
c  ##  external section
c
      external btb, bbt
c
      save qcsave
c
c-----------------------------------------------
      one = 1d0
      maxstep = 1d-1
      maxiter = 1000
      linonly = 0
c # qq1 is the current set of internal coordinates,
c # qq is the original set
c # qq1 will be replaced by the current error
      toler=1.0d-9
      if (linonly.eq.1)maxiter=2
c
c ##  header information
c
      write(nlist,'(//,5x,60(1h*)/)')
      write(nlist,
     & '(10x,''Beginning the iterative coordinate transformation:'')')
      write(nlist,'(12x,''cartesian --> internal coordinates''//)')
      if (linonly.eq.1) write(nlist,
     & '(10x,''Performing linear calculation only.'')')
      write(nlist,'(5x,''Main parameters:'')')
      write(nlist,'(10x,''Number of iterations:'',t45,i3)')maxiter
      write(nlist,'(10x,''Convergance criterium:'',t45,e13.8)')toler
      write(nlist,
     & '(10x,''Max. int. coor. change in 1 iter:'',t45,f6.3)')
     & maxstep
c
c   define temporary storage
      call getmem(nek,ix)
      call getmem(nek,iy)
      call getmem(nek,idm1)
      call getmem(nek,it)
      call getmem(nek,ip)
      call getmem(nek,ir)
      call getmem(nek,iz)
      call getmem(nek,iu)
c
 500  lincalc=linonly
c
c  in the redundant case, first project the internal coordinate
c  changes
C  THIS DOES NOT GO - it projects OK but that changes the other
C  directions if there is redundancy. I do not understand!
c     call bdiag(nq,nek,b,ibcontr,bl(idm1))
c     call bbt(nq,nek,b,ibcontr,qc,bl(iu),t)
c     call congrad(nq,nek,b,ibcontr,qc,bl(iy),bl(idm1),bl(iu),
c    1  bl(it),bl(ip),bl(ir),bl(iz),bbt)
c  end projecting delta q
c  this should be done only if the coordinate system is redundant
c  calculate the diagonal of the matrix BtB
      call bdiag1(nq,nek,b,ibcontr,bl(idm1))
c qq1 is the actual geometry, on start set to qq
      call dcopy_wr(nq,qq,1,qq1,1)
c qc holds the desired set:
      call dcopy_wr(nq,qc,1,qcsave,1)
      call daxpy_wr(nq,one,qq,1,qc,1)
c xy is the current set of Cartesians
      call dcopy_wr(nek,xa,1,xy,1)
c   set the first value of the change
      call dcopy_wr(nq,qc,1,bl(ix),1)
c
      write(nlist,'(/,10x,''Geometries:'')')
      write(nlist,'(12x,''starting'',9x,''change'',11x,''desired'')')
      write(nlist,'(8x,50(1h-))')
      do  i=1,nq
      write(nlist,'(8x,f13.8,2(4x,f13.8))')qq(i),qcsave(i),qc(i)
      enddo ! do i=1,nq
      write(nlist,*)
c
c  #------------------------------------------------
c  ##  beginning the main loop
c
      do irep=1,maxiter
c
        if (lincalc.ne.1) then
         write(nlist,
     &    '(/,2x,5(1h-),2x,''Beginning iteration: '',i3,2x,5(1h-)/)')
     &    irep
        endif ! if (lincalc.ne.1)
       write(nlist,'(10x,''Current geometry:'')')
       write(nlist,'(8x,f16.10)')(qq1(i),i=1,nq)
c
c  calculate the error in the internal coordinates, put it in qq1
      call dscal_wr(nq,onem,qq1,1)
      call daxpy_wr(nq,one,qc,1,qq1,1)
c  do not calculate actual - should values
c  rather, look at the change
      call absmax(nq,bl(ix),imax,delmax)
      write(nlist,'(10x,''Current geometry deviation: '',e13.8)') delmax
c
c  # for irep=2
      if ((lincalc.eq.1).and.(irep.eq.2)) then
       write(nlist,'(10x,''Performed linear calculation'')')
       go to 1000
      endif ! if ((lincalc.eq.1).and.(irep.eq.2))
c
c  # check the convergence and quit
        if(delmax.lt. toler) then
         write(nlist,'(/,10x,''* * *  Convergence reached  * * *'')')
         go to 1000
        endif ! if(delmax.lt. toler)
      call absmax(nq,qq1,imax,delmax)
c
c  # scale back the change if neccessery
        if((delmax.gt.maxstep).and.(lincalc.ne.1)) then
         scale=maxstep/delmax
         call dscal_wr(nq,scale,qq1,1)
         write(nlist,
     &    '(10x,''Step scalling active, actual steps size:'')')
         do  i=1,nq
         write(nlist,'(28x,f13.8)')qq1(i)
         enddo ! do i=1,nq
        else
         scale=one
        endif ! if((delmax.gt.maxstep).and.(lincalc.ne.1))
c
c  delta X = Bt * (BBt)**-1 * Delta Q  to first order
c  first evaluate  x=[(BBt)**-1]Delta Q
c  Interface   subroutine congrad(nq,nek,b,ibcontr,x,y,dm1,c,t,p,r,z,op)
c   alllocate memory for the temporary stuff
c     call congrad(nq,nek,b,ibcontr,bl(ix),bl(iy),bl(idm1),qq1,
c    1  bl(it),bl(ip),bl(ir),bl(iz))
c   bl(ix) holds (BBt)**-1*qq1
c interface subroutine btxv(nq,nek,b,ibcontr,v,w): gives w=Bt*v
c     call btxv(nq,nek,b,ibcontr,bl(ix),bl(iy))
c  bl(iy) now holds the first-order change in the Cartesians
c     only for checking - calc. B*bl(iy)
c     call bxv(nq,nek,b,ibcontr,bl(iy),bl(it))
c     call add1(qq1,onem,bl(it),nq)
c     call absmax(nq,bl(it),ii,dmax)
c     print *, 'devi=', ii,dmax
c   end of testing
c
c   new method: x is the conj. grad solution of Bt*B*x=Bt*q
      call btxv(nq,nek,b,ibcontr,qq1,bl(iu))
      call congrad(nek,nq,b,ibcontr,bl(ix),bl(iy),bl(idm1),bl(iu),
     & bl(it),bl(ip),bl(ir),bl(iz),btb,condamp)
c     do 150 jj=1,nek
c 150   write(nlist,('f18.10'))bl(ix-1+jj)
c     write(nlist,*)
c   re-introduce the max. Cartesian idea - it slows things down
c     call absmax(nek,bl(iy),ii,cmax)
c       if(cmax.gt.cartmax) then
c        scale=cartmax/cmax
c        call mult(bl(iy),scale,nek)
c       end if
      call daxpy_wr(nek,one,bl(ix),1,xy,1)
c  calculate the "should" value of the internal coordinates
c  this is  the original value plus the scaled change
      sc1=one-one/scale
      call dcopy_wr(nq,qc,1,bl(ix),1)
      call daxpy_wr(nq,sc1,qq1,1,bl(ix),1)
c     this call to the b matrix routine only serves to calculate the
c    coordinates
c   save the old values of the int. coord. in bl(ix) with negative sign
      call dcopy_wr(nq,qc,1,bl(ix),1)
      call daxpy_wr(nq,onem,qq1,1,bl(ix),1)
      call machbnew(
     & na,      xy,       nq,  .true.,    shftpt,
     & ibcode,  ibcontr,  b,   qq1)
      do 720 i=1,nq
        call fixtor(bl(ix-1+i),qq1(i))
  720 continue
c   calculate change
      call daxpy_wr(nq,onem,qq1,1,bl(ix),1)
c
      enddo ! do irep=1,maxiter
c    end of the main loop
c---------------------------------------
c
c  ##  convergence not reached
c
c  #  if allowed, perform linear calculation:
        if (linallow.eq.1) then
         linonly=1
         call dcopy_wr(nq,qcsave,1,qc,1)
         write(nlist,
     & '(//5x,'' * * *  Calculation not converged  * * *'')')
         write(nlist,'(10x,''performing linear calculation'')')
         write(nlist,'(10x,''Set flag: linonly= '',i1)')linonly
         go to 500
        endif ! if (linallow.eq.1)
c
c  #  nothing helps => quit
      call bummer(
     & 'Calculation not converged, iter=',irep,faterr)
c
 1000 continue
      write(nlist,
     & '(/10x,''End of the iterative coordinate transformation'')')
      write(nlist,'(/,5x,60(1h*)/)')
c
      call retmem(8)
c
      return
      end
c
c*****************************************************************+
c
      subroutine fixtor (should,found)
c      fix torsional angles if they turn around fully
      implicit real*8 (a-h,o-z)
      dimension xnfo(8)
      common /number/ zero,half,one,two,three,four,five,ten,ten6,tenm8,
     & pi,acc
c
      real*8     scalup,       scaldn
      parameter( scalup=1.1d0, scaldn=0.9d0)
c
      data xnfo/1.0d0,2.0d0,3.0d0,4.0d0,6.0d0,9.0d0,4.5d0,1.5d0/
      deviat=should-found
c       deviat=qq(i)+ccin(i)-qq1(i)
c
c      write(*,*) 'orig, change,present,diff', qq(i),ccin(i),qq1(i),
c     1 deviat
c
c pp
c        the purpose of the following code is to identify cases when a
c        torsional angle changes suddenly by 2*pi
c        because of the 1/n normalization, the change may be 2*pi,
c        2*pi/2, /3, /4, /6 or /9 or /4.5 or /1.5
      adevi=abs(deviat)
      if(adevi.gt.half) then
         do  1285 nf=1,8
            xnfold=xnfo(nf)
c             this is the number of torsional angles linearly combined
c             it can be 1,2,3,4,6 or 9
            offset=two*pi/xnfold
            upper=scalup*offset
            xlower=scaldn*offset
            if(adevi.lt.upper .and. adevi.gt.xlower) then
               if(deviat.gt.zero) then
c                  qq1(i)=qq1(i)+offset
                  found=found+offset
                  go to 1286
               else
c                  qq1(i)=qq1(i)-offset
                  found=found-offset
                  go to 1286
               endif
            endif
1285     continue
      endif
1286  continue
      return
c
      end
c
c***************************************************************
c
      SUBROUTINE GETMEM(AMOUNT,ADDR)
C  THIS SUBROUTINE RESERVES AMOUNT SPACE IN THE COMMON BLOCK
C  BL, CHECKS FOR MAX. MEMORY, AND RETURNS AN INTEGER ADDRESS
C  BL(ADDR) TO BL(ADDR+AMOUNT-1) IS RESERVED.
C  USE: E.G. TO RESERVE 3 MATRICES,  A(N,N), B(N,3*K),C(N,3*K),
C  AND CALCULATE C=A*B. LET US ASSUME THAT ONLY C IS NEEDED LATER,
C  THEN GIVE C THE LOWER ADDRESS!!!
C  IN THE MAIN (CALLING) PROGRAM:
C     CALL GETMEM(N*3*K,IC)
C     CALL GETMEM(N**2,IA)
C     CALL GETMEM(N*3*K,IB)
C     CALL CALCCMX(BL(IA),BL(IB),BL(IC),N,3*K)
C  ....
C     IN THE CALLED PROGRAM CALCCMX:
C     SUBROUTINE CALCCMX(A,B,C,N,K3)
C     IMPLICIT REAL*8 (A-H,O-Z)
C     DIMENSION A(N,N),B(N,K3),C(N,K3)
C     ... GIVE VALUES TO A AND B
C     CALL SETA(A,N)
C     CALL SETB(B,N,K3)
C     CALL MTRXMUL(A,B,C,N,N,K3)
C     CALL RETMEM(2)
C      THESE LAST CALL FREES THE SPACE TAKEN UP BY A AND B
      INTEGER AMOUNT,ADDR
      COMMON /MMANAG/NREQ,MAXMEM,NADR(400)
      COMMON /GANZ/ LCORE,IOV,LAST,LFLAG(4),INUC,IBAS,NA,NBF,NSH,NCF,NCS
     1,NSY(4),NSYM,NGANZ(35),LOPT(30)
c
c     # bummer error types.
      integer   wrnerr,  nfterr,  faterr
      parameter(wrnerr=0,nfterr=1,faterr=2)
c
c-------------------------------------------------------------
C     NREQ, MAXMEM AND NADR(1) MUST BE INITIALIZED TO 0
C     SEE RETALL
C     LCORE IS THE MAXIMUM AVAILABLE MEMORY
      NREQ=NREQ+1
      ADDR=NADR(NREQ)+1
      IXX=ADDR+AMOUNT-1
      IF (IXX.GT.LCORE) THEN
         WRITE(*,*) 'MEMORY OVERFLOW,NEEDED',IXX,'AVAILABLE',LCORE
c         call bummer('MEMORY OVERFLOW,NEEDED',ixx,faterr)
      END IF
      IF (IXX.GT.MAXMEM) MAXMEM=IXX
      NADR(NREQ+1)=IXX
C
      END
c
c*****************************************************************
c
      SUBROUTINE INIT
C
C     .... PRELIMINARY INPUT PROCESSING,SET INITIAL VALUES
C
c      File description (default values)
c
c      file name       file name variable     file number
c    -----------------------------------------------------
c      cart2intls         c2ils             nlist =     8
c      cart2intin         c2iinp            inp   =     30
c      geom               geomfl            ngeom =     71
c      cartgrd            gradfl            ngrad =     72
c      intcfl             intcfl            nintc =     73
c      intgeom            intgeomfl         nintgeom =  20
c      intforce           intforcfl         nintforc =  21
c      geom.new           geomnewfl         ngeomnew =  22
c      intgeomch          intgeomdhfl       nintgeomch= 23
c      intgrad            intgradfl         nintgrad =  24
c      bmatrix            bmatfl            nbmat =     25
c
c

       implicit none
c   mq is the max. number of internal coordinates,
      integer mq
      parameter (mq=594)

c  ##  parameter & common section
c
c  # styles
      character*10 geomst,gradst
      common /style/ geomst,gradst
c
c  # formats
      character*30    geomfm, gradfm, intforcfm, intgeomfm, intgradfm
      common /format/ geomfm, gradfm, intforcfm, intgeomfm, intgradfm
c
c  # unit numbers
      integer       inp,       nlist,    ngeom,    ngrad,      nintc,
     &              nintgeom,  nintforc, ngeomnew, nintgeomch, nintgrad,
     &              nbmat,     ngconv
      common /tape/ inp,       nlist,    ngeom,    ngrad,      nintc,
     &              nintgeom,  nintforc, ngeomnew, nintgeomch, nintgrad,
     &              nbmat,     ngconv
c
c  # file names
      character*40     c2ils,       c2iinp,    geomfl,    gradfl,
     &                 intcfl,      intgeomfl, intforcfl, geomnewfl,
     &                 intgeomchfl, intgradfl, bmatfl,    gconvfl
      common /flnames/ c2ils,       c2iinp,    geomfl,    gradfl,
     &                 intcfl,      intgeomfl, intforcfl, geomnewfl,
     &                 intgeomchfl, intgradfl, bmatfl,    gconvfl
c global optimization parameters
c #MS added these so that shftpt could be initialized properly. 5/10/06
      logical writ,lgdi,lfdi,lmur
      integer itopti,mxopti,iangs,igeom,ifmat,isad,nfix,
     &        ifix,ngeo,ipu
      real*8  thract, shftpt,chmax,coormx
      common /opti/    itopti, mxopti,   chmax,  thract,  shftpt,
     &                 coormx, iangs,    igeom,  ifmat,   isad,
     &                 nfix,   ifix(mq), ngeo,   lgdi,    lfdi,
     &                 lmur,    writ,    ipu


      real*8 ang,debye,cbm,ajoule,evolt,ckalmo,dkel,cmm1,hertz
      common /unit/ ang,debye,cbm,ajoule,evolt,ckalmo,dkel,cmm1,hertz
c
c    constants
      real*8 zero,half,one,two,three,four,five,ten,ten6,tenm8,pi,acc
      common /number/ zero,half,one,two,three,four,five,ten,ten6,tenm8,
     &                pi,acc
c
c     common numbers
c
      zero  = 0.0d0
      half  = 0.5d0
      one   = 1.0d0
      two   = 2.0d0
      three = 3.0d0
      four  = 4.0d0
      five  = 5.0d0
      ten   = 10.0d0
      ten6  = 1.0d+6
      tenm8 = 1.0d-8
      pi    = 3.14159265358979323846264338328d0
      shftpt= pi/two
c
c     common tape
c
      nlist      = 8
      inp        = 30
      ngeom      = 71
      ngrad      = 72
      nintc      = 73
      nintgeom   = 20
      nintforc   = 21
      ngeomnew   = 22
      nintgeomch = 23
      nintgrad   = 24
      nbmat      = 25
      ngconv     = 26
c
c    default file names
c
      c2ils       = 'cart2intls'
      c2iinp      = 'cart2intin'
      geomfl      = 'geom'
      gradfl      = 'cartgrd'
      intcfl      = 'intcfl'
      intforcfl   = 'intforc'
      intgradfl   = 'intgrad'
      intgeomfl   = 'intgeom'
      geomnewfl   = 'geom.new'
      intgeomchfl = 'intgeomch'
      bmatfl      = 'bmatrix'
      gconvfl     = 'geoconvfl'
c
c  ##  formats
c
      geomfm    = '(1x,A2,2x,f5.1,4f14.8)'
      gradfm    = '(3d15.7)'
      intforcfm = '(f15.8)'
      intgradfm = '(f15.8)'
      intgeomfm = '(f15.8)'
c
c  ##  styles
c
      geomst = 'texas     '
      gradst = 'texas     '
c
c     common units
c
      ang    = 1.889725988d0
      debye  = 0.39342658d0
      cbm    = 0.117946332d30
      ajoule = 0.229371044d0
      evolt  = 0.036749026d0
      ckalmo = 1.5936018d-3
      dkel   = 3.1667911d-6
      cmm1   = 4.5563352d-6
      hertz  = 1.51982983d-16
c
      return
      end
c
c**********************************************************
c
      subroutine intforc(nq,nek,b,ibcontr,f,fi,qq,mxf,fmax1)
c    nq (input) number of coordinates
c    nek (input) 3 times the number of atoms (# of cartesiians)
c    b(54,nq) (input) b matrix stored compressed
c    ibcontr(20,nq) coding of b
c    lq,mq dimensions of b (in the main program)
c    f(1..nek) (input) cartesian forces
c    fi(1..nq) (output) descartes forces
c    nlist,(input)  output files
c    qq (input) coordinates
c    fmax1: maximum internal force
c
       implicit none
c  ##  parameter & common block section
c
c   ma is the max. number of atoms
      integer ma
      parameter (ma=200)
c   lq is the max. number of cartesisans
      integer lq
      parameter (lq=3*ma)
c   mq is the max. number of internal coordinates,
      integer mq, maxcard
      parameter (mq=lq-6)
c
      character*10 calctype
      integer          intgradout,intforceout,rgfcoord,fixc(mq),
     &                 geomch
      common /control/ intgradout,intforceout,rgfcoord,fixc,geomch,
     &                 calctype
c
c   constants
      real*8 zero,half,one,two,three,four,five,ten,ten6,tenm8,acc,pi
      common /number/ zero,half,one,two,three,four,five,ten,ten6,tenm8,
     &                pi, acc
c
c  # unit numbers
      integer       inp,       nlist,    ngeom,    ngrad,      nintc,
     &              nintgeom,  nintforc, ngeomnew, nintgeomch, nintgrad,
     &              nbmat,     ngconv
      common /tape/ inp,       nlist,    ngeom,    ngrad,      nintc,
     &              nintgeom,  nintforc, ngeomnew, nintgeomch, nintgrad,
     &              nbmat,     ngconv
c
      integer lcx
      parameter (lcx=10 000)
      real*8 bl
      common bl(lcx)
c
      real*8 condamp
      common /options/condamp
ctm
c
      integer nq,ibcontr(20,nq),i,ix,it,iy,idm,ic,ip,ir,iz,icount
      integer mxf
      integer nek
c
      real*8 b(54,nq)
      real*8 f(nek),fi(nq),fmax1
      real*8 qq(nq)
c
      external bbt
c
      call getmem(nq,ix)
      call getmem(nq,iy)
      call getmem(nq,idm)
      call getmem(nq,ic)
      call getmem(nek,it)
      call getmem(nq,ip)
      call getmem(nq,ir)
      call getmem(nq,iz)
      call bxv(nq,nek,b,ibcontr,f,bl(ic))
      call bdiag(nq,b,ibcontr,bl(idm))
c     do 100 i=1,nek
c100    write(nlist,('f16.10'))bl(ix-1+i)
      call congrad(nq,nek,b,ibcontr,bl(ix),bl(iy),bl(idm),bl(ic),
     1    bl(it),bl(ip),bl(ir),bl(iz),bbt,condamp)
      call dcopy_wr(nq,bl(ix),1,fi,1)
c
c ##  set fixc forces to zero
      icount = 0
 200  icount = icount + 1
      if ( fixc(icount).eq.0) go to 300
      fi(fixc(icount))=zero
      write(nlist,*)
     & '  Internal force on coord.',icount,
     & ' set to zero due to fixc input'
      go to 200
 300  continue
c
ctm   call twocol (nlist,nq,qq,fi,'internal coordinates and forces')
      call onecol (nlist,nq,qq,fi,'internal coordinates and forces')
      call absmax(nq,fi,mxf,fmax1)
      fmax1=sign(fmax1,fi(mxf))
      write(nlist,400) mxf,fmax1
 400  format(' maximum force on coord. ',i4, ' = ',f12.7)
 310  continue
      call retmem(8)
      end
c
c***************************************************************
c
      subroutine machbnew(
     &  na,       xa,       nq,    qonly,    shftpt,
     &  ibcode,   ibcontr,  bmat,  qq)
c
c     parameters: input
c     nek=3*na, na=number of nuclei
c     xa(3,*): nuclear coordinates (in Angstrom)
c     nq=number of internal coordinates
c     qonly: if .true., calculate only the coordinate values
c     shftpt: a constant which determines at which point does
c       a torsional coordinate change by 2*pi
c     ibcode(6,nprim): the encoding of primitive internal coordinates
c     ibcontr(20,nq): contains the contraction pattern, the
c     total number of atoms and the atoms participating in the
c     coordinate (max. 18)
c                 output
c     bmat(54,nq): contains the B matrix elements (for at most 18 atoms=
c         54 Cartesians). Note that this is in a way the transpose of
c         B since the SECOND subscript is the internal coordinate
c     qq(nq): the values of the internal coordinates
c
      IMPLICIT REAL*8 (A-H,O-Z)
c
c  ##  parameter and comon block section
c
c   ma is the max. number of atoms
      integer ma
      parameter (ma=200)
c   lq is the max. number of cartesisans
      integer lq
      parameter (lq=3*ma)
c   mq is the max. number of internal coordinates,
      integer mq, maxcard
      parameter (mq=lq-6)
c
c
c  # number of types of internal coordinates
      integer ntyp
      parameter(ntyp=8)
c
c  # defines the internal coordinate type
      integer             intctype(mq)
      CHARACTER*4         tipus(ntyp)
      common /icoorclass/ intctype
      common /ccoorclass/ tipus
c
      real*8    zero,     one,     two
      parameter(zero=0d0, one=1d0, two=2d0)
      real*8    small,      thous,     six
      parameter(small=1d-8, thous=1d3, six=6d0)
      real*8    pi
      parameter(pi=3.14159265358979323846264338328d0)
      real*8    twopi
      parameter(twopi=two*pi)
c
      integer na, nq
      integer ia(4),ibcode(6,*),ibcontr(20,nq)
c
      real*8 xa(3,*),qq(nq),bmat(54,nq)
      real*8 U(3), V(3), W(3), Z(3), X(3), UU(3), VV(3),
     &       WW(3), ZZ(3), UV(12)
c
      EQUIVALENCE (KA,IA(1)), (KB,IA(2)), (KC,IA(3)), (KD,IA(4))
      EQUIVALENCE (UV(1),UU(1)), (UV(4),VV(1)), (UV(7),WW(1)),
     &            (UV(10),ZZ(1))
c
      LOGICAL QONLY
c
c-----------------------------------------------------------
      call izero_wr(mq,intctype,1)
      call wzero(nq,qq,1)
      do 1000 i=1,nq
         if(i.eq.1) then
            iprim1=1
         else
            iprim1=ibcontr(1,i-1)+1
         end if
         iprim2=ibcontr(1,i)
         iatoms=ibcontr(2,i)
         if(.not.qonly) then
            call wzero(3*iatoms,bmat(1,i),1)
         end if
         do 900  ipr=iprim1,iprim2
            itype=ibcode(1,ipr)
            c = ibcode(2,ipr) * small
            ka=ibcode(3,ipr)
            kb=ibcode(4,ipr)
            kc=ibcode(5,ipr)
            kd=ibcode(6,ipr)
            intctype(i)=itype
            GO TO (190,200,205,210,230,260,280,300), itype
C
C..... STRETCH
C
190         CALL VEKTOR1(UU,R1,KA,KB,XA)
            VV(1)=-UU(1)
            VV(2)=-UU(2)
            VV(3)=-UU(3)
            QQ(I)=QQ(I)+R1*C
            GO TO 320
C
C.....INVERSE
C
200         CALL VEKTOR1(UU,R1,KA,KB,XA)
            RM1=ONE/R1
            RM2=RM1**2
            UU(1)=-RM2*UU(1)
            UU(2)=-RM2*UU(2)
            UU(3)=-RM2*UU(3)
            VV(1)=-UU(1)
            VV(2)=-UU(2)
            VV(3)=-UU(3)
            IA(3)=0
            IA(4)=0
            QQ(I)=QQ(I)+RM1*C
            GO TO 320
C  ... inverse sixth power multiplied by 100
205         CONTINUE
            UU(1)=XA(1,KA)-XA(1,KB)
            UU(2)=XA(2,KA)-XA(2,KB)
            UU(3)=XA(3,KA)-XA(3,KB)
            RM2=ONE/(UU(1)**2+UU(2)**2+UU(3)**2)
            RM6=(RM2**3) * thous
            RM8= -RM6 * RM2 * six
            UU(1)=UU(1)*RM8
            UU(2)=UU(2)*RM8
            UU(3)=UU(3)*RM8
            VV(1)=-UU(1)
            VV(2)=-UU(2)
            VV(3)=-UU(3)
            IA(3)=0
            IA(4)=0
            QQ(I)=QQ(I)+RM6*C
            GO TO 320
C.....BENDING
C
210         CALL VEKTOR1(U,R1,KA,KC,XA)
            CALL VEKTOR1(V,R2,KB,KC,XA)
            CO=ddot_wr(3,u,1,v,1)
            SI=S2(CO)
            DO 220 L=1,3
               UU(L)=(CO*U(L)-V(L))/(SI*R1)
               VV(L)=(CO*V(L)-U(L))/(SI*R2)
               WW(L)=-UU(L)-VV(L)
220         CONTINUE
            IA(4)=0
            QQ(I)=QQ(I)+C*DARCOS(CO)
            GO TO 320
C
C.....OUT OF PLANE
C
230         CALL VEKTOR1(U,R1,KA,KD,XA)
            CALL VEKTOR1(V,R2,KB,KD,XA)
            CALL VEKTOR1(W,R3,KC,KD,XA)
            CALL NORMAL (V,W,Z)
            STETA=ddot_wr(3,u,1,z,1)
            CTETA=S2(STETA)
            CFI1=ddot_wr(3,V,1,W,1)
            SFI1=S2(CFI1)
            CFI2=ddot_wr(3,W,1,U,1)
            CFI3=ddot_wr(3,V,1,U,1)
            DEN=CTETA*SFI1**2
            ST2=(CFI1*CFI2-CFI3)/(R2*DEN)
            ST3=(CFI1*CFI3-CFI2)/(R3*DEN)
            DO 240 L=1,3
               VV(L)=Z(L)*ST2
               WW(L)=Z(L)*ST3
240         CONTINUE
            CALL NORMAL (Z,U,X)
            CALL NORMAL (U,X,Z)
            DO 250 L=1,3
               UU(L)=Z(L)/R1
               ZZ(L)=-UU(L)-VV(L)-WW(L)
250         CONTINUE
            CX=-C
            IF (STETA .LT. zero) CX = C
            QQ(I) = QQ(I) - CX * DARCOS(CTETA)
            GO TO 320
C
C..... TORSION
C
260         CALL VEKTOR1(U,R1,KA,KB,XA)
            CALL VEKTOR1(V,R2,KC,KB,XA)
            CALL VEKTOR1(W,R3,KC,KD,XA)
            CALL NORMAL (U,V,Z)
            CALL NORMAL (W,V,X)
            CO=ddot_wr(3,U,1,V,1)
            CO2=ddot_wr(3,V,1,W,1)
            SI=S2(CO)
            SI2=S2(CO2)
            DO 270 L=1,3
               UU(L)=Z(L)/(R1*SI)
               ZZ(L)=X(L)/(R3*SI2)
               VV(L)=(R1 * CO/R2 - one) * UU(L) - R3 * CO2/R2 * ZZ(L)
               WW(L)=-UU(L)-VV(L)-ZZ(L)
270         CONTINUE
            CO=ddot_wr(3,Z,1,X,1)
            U(1)=Z(2)*X(3)-Z(3)*X(2)
            U(2)=Z(3)*X(1)-Z(1)*X(3)
            U(3)=Z(1)*X(2)-Z(2)*X(1)
            SI3=SQRT(U(1)**2+U(2)**2+U(3)**2)
            CO2=ddot_wr(3,U,1,V,1)
            S=ARC1(-CO,SI3)
            IF (CO2 .LT. zero) S = -S
            IF (S .GT. SHFTPT) S = S - twopi
            QQ(I) = QQ(I) - C * S
C
C.... REMEMBER THAT THE RANGE OF THIS COORDINATE IS -PI/2 TO 3*PI/2
C.... IN ORDER TO SHIFT THE DISCONTINUITY OFF THE PLANAR POSITION
C
            GO TO 320
C
C.....LINEAR COPLANAR BENDING
C
280         CALL VEKTOR1(U,R1,KA,KC,XA)
            CALL VEKTOR1(V,R2,KD,KC,XA)
            CALL VEKTOR1(X,R2,KB,KC,XA)
            CO=ddot_wr(3,V,1,U,1)
            CO2=ddot_wr(3,X,1,V,1)
            QQ(I) = QQ(I) + C * (pi - DARCOS(CO) - DARCOS(CO2))
            CALL NORMAL (V,U,W)
            CALL NORMAL (U,W,Z)
            CALL NORMAL (X,V,W)
            CALL NORMAL (W,X,U)
C
C..... COORDINATE POSITIV IF ATOM A MOVES TOWARDS ATOM D
C
            DO 290 L=1,3
               UU(L)=Z(L)/R1
               VV(L)=U(L)/R2
               WW(L)=-UU(L)-VV(L)
290         CONTINUE
            IA(4)=0
            GO TO 320
C
C.....LINEAR PERPENDICULAR BENDING
C
300         CALL VEKTOR1(U,R1,KA,KC,XA)
            CALL VEKTOR1(V,R2,KD,KC,XA)
            CALL VEKTOR1(Z,R2,KB,KC,XA)
            CALL NORMAL (V,U,W)
            CALL NORMAL (Z,V,X)
            DO 310 L=1,3
               UU(L)=W(L)/R1
               VV(L)=X(L)/R2
               WW(L)=-UU(L)-VV(L)
310         CONTINUE
            IA(4)=0
            CO=ddot_wr(3,U,1,W,1)
            CO2=ddot_wr(3,Z,1,W,1)
            QQ(I) = QQ(I) + C * (pi - DARCOS(CO) - DARCOS(CO2))
320         IF (QONLY) GO TO 900
            DO 340 J=1,4
               M=IA(J)
               IF (M.LE.0) GO TO 340
               iatoms=ibcontr(2,i)
               do 330 kk=1,iatoms
c   search the atoms participating in the coordinate. Is one of them
c   (ikk) the same as the atom
                  ikk=ibcontr(kk+2,i)
                  if(ikk.eq.m) then
                     imm=kk
                     go to 335
                  end if
330            continue
               print *, 'this should not happen, i,ibcontr,ia(1..4),j',
     1          i,(ibcontr(kk,i),kk=1,iatoms+2),ia,j
335            continue
               j3=j*3
*@ifdef debug
*      write(*,*) 'bmat(imm*3..,',i,')=',bmat(imm*3-2,i),
*     .   bmat(imm*3-1,i),bmat(imm*3,i)
*      write(*,*) 'uv(*)=',uv(j3-2),uv(j3-1),uv(j3), 'c=',c,' j3=',j3
*@endif
               bmat(imm*3-2,i)=bmat(imm*3-2,i)+uv(j3-2)*c
               bmat(imm*3-1,i)=bmat(imm*3-1,i)+uv(j3-1)*c
               bmat(imm*3,i)=bmat(imm*3,i)+uv(j3)*c
340         CONTINUE
900      continue
1000  continue
c     write(*,*)'coor, type'
c     do i=1,nq
c     write(*,*)i,intctype(i)
c     enddo
      return
      END
c
c********************************************************************
c
      subroutine nerror(noer,routine,message,n1,n2)
      character*(*) routine
      character*(*) message
c  # unit numbers
      integer       inp,       nlist,    ngeom,    ngrad,      nintc,
     &              nintgeom,  nintforc, ngeomnew, nintgeomch, nintgrad,
     &              nbmat,     ngconv
      common /tape/ inp,       nlist,    ngeom,    ngrad,      nintc,
     &              nintgeom,  nintforc, ngeomnew, nintgeomch, nintgrad,
     &              nbmat,     ngconv
c
c---------------------------------------------------------
c
      write(nlist,*) 'Error no. ',noer,' in ',routine,' ',message,n1,n2
c100  format(' Error No.',i4,' in ',a20,2x,/,a80,/,' variables ',2i6)
c     new line by ph
c     bummmer instead of stop 20
      call bummer('error in cart2int, error no.=',noer,2)
      end
c
c********************************************************************
c
      SUBROUTINE NOM (U)
c
       implicit none
      real*8 one
      parameter (one=1.0d+00)
      real*8 u(3),x
      integer i
      real*8 ddot_wr
      external ddot_wr
c
      X=one/SQRT(ddot_wr(3,U,1,U,1))
         DO I=1,3
         u(i)=u(i)*x
         enddo
      RETURN
      end
c
c********************************************************************
c
      SUBROUTINE NORMAL (U,V,W)
      IMPLICIT none
      real*8 U(3), V(3), W(3)
C
C     99999...  W WIRD EIN SENKRECHT AUF DIE EBENE(U,V) STEHENDER EINHEI
C      TOR
C
      W(1)=U(2)*V(3)-U(3)*V(2)
      W(2)=U(3)*V(1)-U(1)*V(3)
      W(3)=U(1)*V(2)-U(2)*V(1)
      CALL NOM (W)
      RETURN
C
      END
c
c*********************************************************
c
      subroutine precon(x,dm1,y,nq)
c  calculates y=dm1*x
       implicit none 
       integer nq,i
      real*8 x(nq),dm1(nq),y(nq)
c
      do i=1,nq
        y(i)=dm1(i)*x(i)
      end do
      end
c
c*********************************************************
c
      subroutine readgrad
c
c     this routine reads the gradient file
c
c     input:
c          na         number of atoms
c          gradfl      file name
c          gradst      what form the file containes geometry
c                          currently known:  TEXAS
c          gradfm      format of the gradient input
c
c     output:
c          f(*)         cartesian forces
c
c   written by Peter Szalay jun. 95
c
       implicit none
c  ##  paramenter & common section
c
c   ma is the max. number of atoms
      integer ma
      parameter (ma=200)
c   lq is the max. number of cartesisans
      integer lq
      parameter (lq=3*ma)
c
c  common related to input data
      character*3 symb
      real*8 xa,etot,f,mass
      integer ia,na
      common /geom/ xa(3,ma),  etot,     f(lq),  mass(ma),  ia(ma),
     & na
      common /cgeom/ symb(ma)
c
      integer           linonly,    linallow,    maxiter, forceinau
      real*8            maxstep
      common/transform/ linonly,    linallow,    maxstep,  maxiter,
     &                  forceinau
c  # styles
      character*10 geomst,gradst
      common /style/ geomst,gradst
c
c  # formats
      character*30    geomfm, gradfm, intforcfm, intgeomfm, intgradfm
      common /format/ geomfm, gradfm, intforcfm, intgeomfm, intgradfm
c
c  # unit numbers
      integer       inp,       nlist,    ngeom,    ngrad,      nintc,
     & nintgeom,  nintforc, ngeomnew, nintgeomch, nintgrad,
     & nbmat,     ngconv
      common /tape/ inp,       nlist,    ngeom,    ngrad,      nintc,
     & nintgeom,  nintforc, ngeomnew, nintgeomch, nintgrad,
     & nbmat,     ngconv
c
c  # file names
      character*40     c2ils,       c2iinp,    geomfl,    gradfl,
     & intcfl,      intgeomfl, intforcfl, geomnewfl,
     & intgeomchfl, intgradfl, bmatfl,    gconvfl
      common /flnames/ c2ils,       c2iinp,    geomfl,    gradfl,
     & intcfl,      intgeomfl, intforcfl, geomnewfl,
     & intgeomchfl, intgradfl, bmatfl,    gconvfl
c
c     # bummer error types.
      integer   wrnerr,  nfterr,  faterr
      parameter(wrnerr=0,nfterr=1,faterr=2)
c
c  ##  integer section
c
      integer i, ii, icount
      integer nek, nunit
c
c  ## real*8 section
c
      real*8 xfqx
c
      real*8     smallf
      parameter( smallf=1d-5 )
c
      logical yesno
c
c-------------------------------------------------------------
c
      nek=3*na
c
      inquire(file=gradfl,opened=yesno,number=nunit)
      if(.not.yesno) then
         nunit=ngrad
         open(unit=nunit,file=gradfl,status='unknown',form='formatted')
      endif
c
      write (nlist,10)
10    format (//,10x,'Cartesian forces',/)
      write(nlist,
     & '(6x,1h#,2x,''Symb'',8x,''fx'',10x,''fy'',10x,''fz'')')
      write(nlist,'(4x,49(1h-))')
c
      if(gradst.eq.'texas     '.or.gradst.eq.'TEXAS     ') then
         read (nunit,gradfm) (f(i),i=1,nek)
         do ii=1,nek
ctm       f(ii)=f(ii)*8.2387295d0
         if (forceinau.eq.1) then
            f(ii)=-f(ii)/0.529177249d+00
         else
            f(ii)=-f(ii)*8.2387295d0
         endif
         enddo                  ! do ii=1,nek
         icount = 0
         do ii=1,na
            write (nlist,101)ii,symb(ii),(f(icount+i),i=1,3)
            icount = icount + 3
         enddo                  ! do ii=1,na
c
101      format (5x,i2,3x,a3,2x,3f12.6)
      else
         write(nlist,*) ' error in readgrad: style ',gradst,
     +    ' is not known'
         write(nlist,*) ' program terminates'
         call bummer('error in readgrad: style ',0,faterr)
      end if
c
      xfqx=f(1)
      do 1100 i=2,nek
         xfqx=xfqx+f(i)
1100  continue
      if (abs(xfqx) .gt. smallf) write (nlist,1101) xfqx
1101  format (/,1x,26hforces do not vanish, sum=,f15.7,/)
c
      return
      end
c
c**************************************************************
c
      subroutine readgeom
c
c     this routine reads the geom geomfl
c
c     input:
c          geomfl       geomfl name
c          geomst       what form the geomfl containes geometry
c                          currently known:  TEXAS
c          geomfm       format of the geometry input
c          iangs        =1   input in angstroms
c                       =0   input in bohr
c          nlist         unit for the listings
c          na           number of atoms (in some cases)
c
c     output:
c          x(*),y(*),z(*)      cartesian coordinates
c          symb(*)             element symbols
c          ia                  atomic number
c          na                  number of atoms (in some cases)
c
c   written by Peter Szalay jun. 95
c
      implicit real*8 (a-h,o-z)
c
c  # styles
      character*10 geomst,gradst
      common /style/ geomst,gradst
c
c   ma is the max. number of atoms
      integer ma
      parameter (ma=200)
c   lq is the max. number of cartesisans
      integer lq
      parameter (lq=3*ma)
c   mq is the max. number of internal coordinates,
      integer mq
      parameter (mq=lq-6)
c
c  common related to input data
      character*3 symb
      real*8 xa,etot,f,mass
      integer ia,na
      common /geom/ xa(3,ma),  etot,     f(lq),  mass(ma),  ia(ma),
     &              na
      common /cgeom/ symb(ma)
c
c
c  # formats
      character*30    geomfm, gradfm, intforcfm, intgeomfm, intgradfm
      common /format/ geomfm, gradfm, intforcfm, intgeomfm, intgradfm
c
c  # unit numbers
      integer       inp,       nlist,    ngeom,    ngrad,      nintc,
     &              nintgeom,  nintforc, ngeomnew, nintgeomch, nintgrad,
     &              nbmat,     ngconv
      common /tape/ inp,       nlist,    ngeom,    ngrad,      nintc,
     &              nintgeom,  nintforc, ngeomnew, nintgeomch, nintgrad,
     &              nbmat,     ngconv
c
c  # file names
      character*40     c2ils,       c2iinp,    geomfl,    gradfl,
     &                 intcfl,      intgeomfl, intforcfl, geomnewfl,
     &                 intgeomchfl, intgradfl, bmatfl,    gconvfl
      common /flnames/ c2ils,       c2iinp,    geomfl,    gradfl,
     &                 intcfl,      intgeomfl, intforcfl, geomnewfl,
     &                 intgeomchfl, intgradfl, bmatfl,    gconvfl
c
c    these are the parameters for the global geometry optimization
      logical writ,lgdi,lfdi,lmur
      integer itopti,mxopti,icoormx,iangs,igeom,ifmat,isad,nfix,
     &        ifix,ngeo,ipu
      real*8  thract, shftpt,chmax,coormx
      common /opti/itopti,mxopti,chmax,thract,shftpt,coormx,iangs,
     +     igeom,ifmat,isad,nfix,ifix(mq),ngeo,lgdi,lfdi,lmur,writ,ipu
c
      common /unit/ ang,debye,cbm,ajoule,evolt,ckalmo,dkel,cmm1,hertz
c   unit numbers
c
      logical yesno
c
c     # bummer error types.
      integer   wrnerr,  nfterr,  faterr
      parameter(wrnerr=0,nfterr=1,faterr=2)
c
c---------------------------------------------------------------
c
      inquire(file=geomfl,opened=yesno,number=nunit)
      if(.not.yesno) then
         nunit=ngeom
         open(unit=nunit,file=geomfl,status='unknown',form='formatted')
      endif
c
      write (nlist,10)
 10   format
     & (//10x,'Initial Cartesian nuclear coordinates in angstrom'/)
      write(nlist,
     & '(6x,1h#,2x,''Symb'',2x,''Nukl'',5x,1hx,11x,1hy,11x,1hz,'//
     & '10x,''mass'')')
      write(nlist,'(4x,65(1h-))')
c
      if(geomst.eq.'texas     '.or.geomst.eq.'TEXAS     ') then
        if(na.eq.0) then
          na=1
 101      continue
          read(nunit,geomfm,end=103) symb(na),xiat,xa(1,na),xa(2,na),
     +         xa(3,na),mass(na)
          ia(na)=xiat+0.1
          if (iangs.ne.1) then
            xa(1,na)=xa(1,na)/ang
            xa(2,na)=xa(2,na)/ang
            xa(3,na)=xa(3,na)/ang
          endif
          write (nlist,102) na,symb(na),ia(na),xa(1,na),xa(2,na),
     &     xa(3,na),mass(na)
 102      format (5x,i2,3x,a3,3x,i2,4f12.7)
          na=na+1
          go to 101
 103      continue
          na=na-1
        else
          do 201 i=1,na
          read(nunit,geomfm) symb(i),xiat,xa(1,i),xa(2,i),xa(3,i)
          ia(i)=xiat+0.1
            if (iangs.ne.1) then
             xa(1,i)=xa(1,i)/ang
             xa(2,i)=xa(2,i)/ang
             xa(3,i)=xa(3,i)/ang
            endif
          write (nlist,102) i,symb(i),ia(i),xa(1,i),xa(2,i),xa(3,i),
     &     mass(i)
  201   continue
        endif
      else
         write(nlist,*) ' error in readgeom: geomst ',geomst,
     +        ' is not known'
         write(nlist,*) ' program terminates'
         call bummer('error in readgeom: geomst',0,faterr)
      end if
c
      if(na.eq.0.or.na.gt.ma) then
         write(nlist,*) ' too many or too few atoms. na=',na,
     +        ' maxatom=',ma
         write(nlist,*) ' program terminates'
         call bummer(' too many or too few atoms. na=',na,faterr)
      endif
      write(nlist,*)
c
      return
      end
c
c*****************************************************************
c
      SUBROUTINE RETALL
C     RESETS THE OCCUPIED MEMORY TO NOTHING
      COMMON /MMANAG/NREQ,MAXMEM,NADR(400)
      NREQ=0
      MAXMEM=0
      NADR(1)=0
      return
      END
c
c*****************************************************************
c
      SUBROUTINE RETMEM(N)
C     REMOVES THE RESERVATION FOR THE LAST OCCUPIED BLOCK
      COMMON /MMANAG/NREQ,MAXMEM,NADR(400)
      DO 100 I=1,N
        IF (NREQ.GT.0) NREQ=NREQ-1
 100  CONTINUE
      END
c
c*****************************************************************
c
      REAL*8 FUNCTION S2 (X)
      IMPLICIT REAL*8 (A-H,O-Z)
      parameter(one=1d0)
      S2 = SQRT(one-X**2)
      RETURN
      END
c
c*****************************************************************
c
      SUBROUTINE TWOCOL(IFILE,N,Q,F,TITLE)
      IMPLICIT REAL*8 (A-H,O-Z)
C     THIS SUBROUTINE PRINTS INTERNAL COORDINATES AND FORCES, OR
C     SIMILAR QUANTITIES, IN TWO COLUMNS
C     PARAMETERS:
C     IFILE: OUTPUT FILE
C     N:(INPUT) NUMBER OF COORDINATES
C     Q(1:N),F(1:N) (INPUT) COORDINATES & FORCES
C     TITLE (CHARACTER STRING), (INPUT) TITLE
      CHARACTER*(*) TITLE
      DIMENSION Q(*),F(*)
      WRITE(IFILE,100) TITLE
 100  FORMAT(1X,A60,/)
      NN=(N+1)/2
      DO 200 I=1,NN
      II=I+NN
      IF (II.LE.N) THEN
        WRITE (IFILE,300) I,Q(I),F(I),II,Q(II),F(II)
      ELSE
        WRITE(IFILE,300) I,Q(I),F(I)
      END IF
 200  CONTINUE
 300  FORMAT(2(1X,I3,2X,F12.7,1X,F12.7))
      END
c
c*****************************************************************
c
      SUBROUTINE ONECOL(IFILE,N,Q,F,TITLE)
      IMPLICIT REAL*8 (A-H,O-Z)
C     THIS SUBROUTINE PRINTS INTERNAL COORDINATES AND FORCES, OR
C     SIMILAR QUANTITIES, IN ONE COLUMNS
C     PARAMETERS:
C     IFILE: OUTPUT FILE
C     N:(INPUT) NUMBER OF COORDINATES
C     Q(1:N),F(1:N) (INPUT) COORDINATES & FORCES
C     TITLE (CHARACTER STRING), (INPUT) TITLE
      CHARACTER*(*) TITLE
      DIMENSION Q(*),F(*)
      WRITE(IFILE,100) TITLE
 100  FORMAT(//1X,A35,/)
      DO 200 I=1,N
        WRITE(IFILE,300) I,Q(I),F(I)
 200  CONTINUE
 300  FORMAT((1X,I3,2X,F12.7,1X,F12.7))
      END



      SUBROUTINE VEKTOR1(U,R,I,J,XA)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION U(3), XA(3,*)
      parameter(one=1d0)
C
C       BILDET DEN NORMIERTEN ENTFERNUNGSVEKTOR VOM KERN J NACH KERN I
C        UND DIE ENTFERNUNG R
C
      U(1)=XA(1,I)-XA(1,J)
      U(2)=XA(2,I)-XA(2,J)
      U(3)=XA(3,I)-XA(3,J)
      R=SQRT(u(1)**2+u(2)**2+u(3)**2)
      rr=one/R
      u(1)=u(1)*rr
      u(2)=u(2)*rr
      u(3)=u(3)*rr
C
      END

      subroutine print_revision(name,nlist)
      implicit none
      integer nlist
      character*12 name
      character*70 string
  50  format('* ',a,t12,a,t40,a,t68,' *')
      write(string,50) name(1:len_trim(name)) // '.f',
     .   '$Revision: 2.7.10.2 $','$Date: 2014/10/09 17:21:33 $'
      call substitute(string)
      write(nlist,'(a)') string
      return
      end


       subroutine progheader(nlist)
       implicit none
       integer nlist
       character*70 str
c
c    Conversion of nuclear coordinates between cartesian and internal
c    coordinates
c
c    Columbus distribution version 5.9
c
c
c    This program converts geometry from->to:
c    cartesian -> internal
c    internal -> cartesian (using small corretion to a reference geom)
c
c    b(54,nq): contains the non-zero elements of B
c    bmat(lq,mq): contains the complete B matrix

      write (nlist,"(30x,'program ""cart2int"" 7.0'/
     &  28x,'columbus program system'//
     &  13x,'this program converts the molecular geometry '/            
     &  13x,'from cartesian to internal coordinates and vice versa'/)")
c
c     # print out the name and address of the local programmer
c     # responsible for maintaining this program.
c
      call who2c( 'ARGOS', nlist )
 20   format(5('****'),'*** File revision status: ***',5('****'))
      write(nlist,20)
      call print_revision('cart2int    ',nlist)
 21   format(17('****'))
      write(nlist,21)
      end

      subroutine substitute(str)
      implicit none
      character*70 str
      integer i
       do i=1, len_trim(str)
         if (str(i:i).eq.'$') str(i:i)=' '
       enddo
      return
      end


