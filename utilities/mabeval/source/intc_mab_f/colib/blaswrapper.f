
*@ifdef essl
*        subroutine dgemul_wr(a,n1,c1,  b,n2,c2,  c,n3,  n4,n5,n6)
*        implicit none
*        real*8 a(*),b(*),c(*)
*        character*1  c1,c2
*@if defined( int64 ) && defined ( noblas64 )
*         integer*8 n1,n2,n3,n4,n5,n6
*         integer*4 ln1,ln2,ln3,ln4,ln5,ln6
*         ln1=n1
*         ln2=n2
*         ln3=n3
*         ln4=n4
*         ln5=n5
*         ln6=n6
*         call dgemul(a,ln1,c1,b,ln2,c2,c,ln3,ln4,ln5,ln6) 
*@else
*         integer n1,n2,n3,n4,n5,n6
*         call dgemul(a,n1,c1,b,n2,c2,c,n3,n4,n5,n6) 
*@endif
*         return
*         end
*@endif




        real*8 FUNCTION DNRM2_WR ( N, X, INCX )
         implicit none
         real*8 dnrm2,X
*@if defined( int64 ) && defined ( noblas64 )
*         integer*8 n,incx
*         integer*4 ln,lincx
*         integer*8 maxval
*         parameter (maxval=2**31-1)
*         if (n.gt.maxval) then
*           call bummer('dnrm2_wr, n exceeds maxval',n,2)
*         endif
*         ln=n
*         lincx=incx 
*         dnrm2_wr = dnrm2(ln,x,lincx)
*@else
         integer n,incx
         dnrm2_wr= dnrm2(n,x,incx)
*@endif
        return
        end

      SUBROUTINE DGETRF_WR( M, N, A, LDA, IPIV, INFO )
c*
c  -- LAPACK routine (version 3.0) --
c     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
c     Courant Institute, Argonne National Lab, and Rice University
c     March 31, 1993
c*
*@if defined( int64 ) && defined ( noblas64 )
*      INTEGER*8            INFO, LDA, M, N
*      INTEGER*8            IPIV( * )
*      INTEGER*4        LINFO,LLDA,LM,LN
*c     INTEGER*4, allocatable :: LIPIV(*) 
*      real*8   A( LDA, * )
*       LM=M
*       LN=N
*       LLDA=LDA
*       call dgetrf(lm,ln,a,llda,ipiv,linfo)
*       info=linfo
*@else
      INTEGER            INFO, LDA, M, N
      INTEGER            IPIV( * )
      real*8   A( LDA, * )
      call dgetrf(m,n,a,lda,ipiv,info)
*@endif

       return
       end

      SUBROUTINE DGETRI_WR( N, A, LDA, IPIV, WORK, LWORK, INFO )
c*
c  -- LAPACK routine (version 3.0) --
c     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
c     Courant Institute, Argonne National Lab, and Rice University
c     June 30, 1999
*
*@if defined( int64 ) && defined ( noblas64 )
*      INTEGER*8            INFO, LDA, LWORK, N
*      INTEGER*8            IPIV( * )
*      real*8   A( LDA, * ), WORK( * )
*      INTEGER*4 linfo,llda,llwork,ln
*       ln=n
*       llda=lda
*       llwork=lwork
*      call dgetri(ln,a,llda,ipiv,work,llwork,linfo)
*      info=linfo
*@else
      INTEGER            INFO, LDA, LWORK, N
      INTEGER            IPIV( * )
      real*8   A( LDA, * ), WORK( * )
      
      call dgetri(n,a,lda,ipiv,work,lwork,info)
*@endif
      return
      end




      SUBROUTINE DGEMM_WR ( TRANSA, TRANSB, M, N, K, ALPHA, 
     .    A, LDA, B, LDB,
     $                   BETA, C, LDC )
c     .. Scalar Arguments ..
      implicit none
      CHARACTER*1        TRANSA, TRANSB
      real*8   ALPHA, BETA
c     .. Array Arguments ..

*@if defined( int64 ) && defined ( noblas64 )
*      INTEGER*8            M, N, K, LDA, LDB, LDC
*      INTEGER*4            LM, LN, LK, LLDA, LLDB, LLDC
*      real*8   A( LDA, * ), B( LDB, * ), C( LDC, * )
*      
*       lm=m
*       ln=n
*       lk=k
*       llda=lda
*       lldb=ldb
*       lldc=ldc
*      
*       call dgemm(TRANSA, TRANSB, lM, lN, lK, ALPHA, A, lLDA, B, lLDB,
*     .    BETA, C, lLDC )
* 
*@else
      INTEGER            M, N, K, LDA, LDB, LDC
       real*8   A( LDA, * ), B( LDB, * ), C( LDC, * )
       call dgemm(TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB,
     .    BETA, C, LDC )
*@endif
       
       return
       end

    

      SUBROUTINE DGEMV_WR ( TRANS, M, N, ALPHA, A, LDA, X, INCX,
     $                   BETA, Y, INCY )
c     .. Scalar Arguments ..
      implicit none
      real*8   ALPHA, BETA
      CHARACTER*1        TRANS
c     .. Array Arguments ..
 
*@if defined( int64 ) && defined ( noblas64 )
*      INTEGER*8            INCX, INCY, LDA, M, N
*      INTEGER*4            lINCX, lINCY, lLDA, lM, lN
*      real*8   A( LDA, * ), X( * ), Y( * )
*       lincx=incx
*       lincy=incy
*       llda = lda
*       lm=m
*       ln=n
*       call  DGEMV ( TRANS, lM, lN, ALPHA, A, lLDA, X, lINCX,
*     $                   BETA, Y, lINCY )
*@else
       INTEGER            INCX, INCY, LDA, M, N
       real*8   A( LDA, * ), X( * ), Y( * )
       call  DGEMV ( TRANS, M, N, ALPHA, A, LDA, X, INCX,
     $                   BETA, Y, INCY )
*@endif
       return
       end


      real*8 function ddot_wr(n,dx,incx,dy,incy)
      implicit none
      real*8 ddot, dx(*),dy(*),dtemp
*@if defined( int64 ) && defined ( noblas64 )
*      integer*8 n,incx,incy
*      integer*4 ln,lincx,lincy
*       
*      ln=n
*      lincx=incx
*      lincy=incy 
*      ddot_wr=ddot(ln,dx,lincx,dy,lincy)
*@else
      integer n,incx,incy
      ddot_wr=ddot(n,dx,incx,dy,incy)
*@endif
      return
      end



      subroutine  dcopy_wr(n,dx,incx,dy,incy)
c
c     copies a vector, x, to a vector, y.
c     uses unrolled loops for increments equal to one.
c     jack dongarra, linpack, 3/11/78.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
      implicit none
      real*8 dx(n),dy(n)
*@if defined( int64 ) && defined ( noblas64 )
*      integer*8 n,incx,incy
*      integer*4 ln,lincx,lincy
*      ln=n
*      lincx=incx
*      lincy=incy 
*      call  dcopy(ln,dx,lincx,dy,lincy)
*@else
      integer n,incx,incy
       call  dcopy(n,dx,incx,dy,incy)
*@endif
      return
      end



      subroutine  dscal_wr(n,da,dx,incx)
c
c     scales a vector by a constant.
c     uses unrolled loops for increment equal to one.
c     jack dongarra, linpack, 3/11/78.
c     modified 3/93 to return if incx .le. 0.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
      implicit none
      real*8 da,dx(*)
*@if defined( int64 ) && defined ( noblas64 )
*      integer*8 n,incx
*      integer*4 ln,lincx
*      ln=n
*      lincx=incx
*      call dscal(ln,da,dx,lincx)
*@else
      integer n,incx
      call dscal(n,da,dx,incx)
*@endif

      return
      end






      SUBROUTINE DGER_WR  ( M, N, ALPHA, X, INCX, Y, INCY, A, LDA )
c     .. Scalar Arguments ..
      implicit none
      real*8   ALPHA
c     .. Array Arguments ..
*@if defined( int64 ) && defined ( noblas64 )
*      INTEGER*8          INCX, INCY, LDA, M, N
*      INTEGER*4          lINCX, lINCY, lLDA, lM, lN
*      real*8   A( LDA, * ), X( * ), Y( * )
*       lincx=incx
*       lincy=incy 
*       llda = lda
*       lm = m
*       ln = n
*       call dger( lM, lN, ALPHA, X, lINCX, Y, lINCY, A, lLDA )
*     
*@else
      INTEGER          INCX, INCY, LDA, M, N
       real*8   A( LDA, * ), X( * ), Y( * )
      call dger( M, N, ALPHA, X, INCX, Y, INCY, A, LDA )
*@endif
      return
      end

      real*8 function dasum_wr(n,dx,incx)
      implicit none
c
c     takes the sum of the absolute values.
c     jack dongarra, linpack, 3/11/78.
c     modified 3/93 to return if incx .le. 0.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
      real*8 dasum,dx(*),dtemp
*@if defined( int64 ) && defined ( noblas64 )
*      integer*8 n,incx
*      integer*4 ln,lincx
*       ln=n
*       lincx=incx
*      dasum_wr=dasum(ln,dx,lincx)
*@else
      integer n,incx
      dasum_wr=dasum(n,dx,incx)
*@endif
      return
      end



      subroutine daxpy_wr(n,da,dx,incx,dy,incy)
c
c     constant times a vector plus a vector.
c     uses unrolled loops for increments equal to one.
c     jack dongarra, linpack, 3/11/78.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
      implicit none
      real*8 dx(*),dy(*),da
*@if defined( int64 ) && defined ( noblas64 )
*      integer*8 n,incx,incy
*      integer*4 ln,lincx,lincy
*       ln=n
*       lincx=incx
*       lincy=incy
*       call daxpy(ln,da,dx,lincx,dy,lincy)
*@else
      integer n,incx,incy
      call daxpy(n,da,dx,incx,dy,incy)
*@endif
      return
      end

c


      subroutine  dswap_wr (n,dx,incx,dy,incy)
c
c     interchanges two vectors.
c     uses unrolled loops for increments equal one.
c     jack dongarra, linpack, 3/11/78.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
      implicit none
      real*8 dx(*),dy(*),dtemp
*@if defined( int64 ) && defined ( noblas64 )
*      integer*8 n,incx,incy
*      integer*4 ln,lincx,lincy
*       ln=n
*       lincx=incx
*       lincy=incy
*       call dswap(ln,dx,lincx,dy,lincy)
*@else
      integer n,incx,incy
      call dswap (n,dx,incx,dy,incy)
*@endif
      return
      end 

      integer function idamax_wr(n,dx,incx)
c
c     finds the index of element having max. absolute value.
c     jack dongarra, linpack, 3/11/78.
c     modified 3/93 to return if incx .le. 0.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
      real*8 dx(*),dmax
      integer i,incx,ix,n
c
      idamax_wr = 0
      if( n.lt.1 .or. incx.le.0 ) return
      idamax_wr = 1
      if(n.eq.1)return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      ix = 1
      dmax = dabs(dx(1))
      ix = ix + incx
      do 10 i = 2,n
         if(dabs(dx(ix)).le.dmax) go to 5
         idamax_wr = i
         dmax = dabs(dx(ix))
    5    ix = ix + incx
   10 continue
      return
c
c        code for increment equal to 1
c
   20 dmax = dabs(dx(1))
      do 30 i = 2,n
         if(dabs(dx(i)).le.dmax) go to 30
         idamax_wr = i
         dmax = dabs(dx(i))
   30 continue
      return
      end


      SUBROUTINE DSPR_WR  ( UPLO, N, ALPHA, X, INCX, AP )
       implicit none
c     .. Scalar Arguments ..
      real*8   ALPHA
c     .. Array Arguments ..
      real*8   AP( * ), X( * )
      CHARACTER*1        UPLO
*@if defined( int64 ) && defined ( noblas64 )
*      INTEGER*8          INCX, N
*      INTEGER*4         LINCX,LN
*      LINCX=INCX
*      LN=N
*      call DSPR  ( UPLO, lN, ALPHA, X, lINCX, AP )
*@else
      INTEGER          INCX, N
      call DSPR  ( UPLO, N, ALPHA, X, INCX, AP )
*@endif
      return
      end

*@ifdef molcas_ext
*      subroutine  icopy_(n,ivx,incx,ivy,incy)
*@else
      subroutine  icopy(n,ivx,incx,ivy,incy)
*@endif
      
       call icopy_WR(n,ivx,incx,ivy,incy)
      return
      end


      subroutine  icopy_WR(n,ivx,incx,ivy,incy)
c
c  blas-style integer vector copy.
c
      implicit integer(a-z)
      integer ivx(*),ivy(*)
      integer i,incx,incy,ix,iy,n
c
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)then
c
c        code for both increments equal to 1
c
         do 30 i = 1,n
            ivy(i) = ivx(i)
30       continue
c
      else
c
c        code for unequal increments or equal increments
c          not equal to 1
c
         ix = 1
         iy = 1
         if(incx.lt.0)ix = (-n+1)*incx + 1
         if(incy.lt.0)iy = (-n+1)*incy + 1
         do 10 i = 1, n
            ivy(iy) = ivx(ix)
            ix = ix + incx
            iy = iy + incy
10       continue
c
      endif
c
      return
      end

c deck izero
      subroutine izero_wr(n,ia,inca)
c
c zero the elements of ia
c
      integer n, inca
      integer ia(inca,*)
c
      integer i
      do 10 i = 1, n
         ia(1,i) = 0
10    continue
c
      return
      end
*@ifdef essl
**
*      subroutine dgemx_wr(m,n,alpha,a,lda,x,incx,y,incy)
*      implicit none
*      real*8 a(*),x(*),y(*), alpha
*@if defined( int64 ) && defined ( noblas64 )
*      integer m,n,lda,incx,incy
*      integer*4 lm,ln,llda,lincx,lincy 
*       lm=m
*       ln=n
*       llda=lda
*       lincx=incx
*       lincy=incy
*      call dgemx(lm,ln,alpha,a,llda,x,lincx,y,lincy)
*@else
*      integer m,n,lda,incx,incy
*      call dgemx(m,n,alpha,a,lda,x,incx,y,incy)
*@endif
*      return
*      end
**
*      subroutine dgemtx_wr(m,n,alpha,a,lda,x,incx,y,incy)
*      implicit none
*      real*8 a(*),x(*),y(*), alpha
*@if defined( int64 ) && defined ( noblas64 )
*      integer m,n,lda,incx,incy
*      integer*4 lm,ln,llda,lincx,lincy 
*       lm=m
*       ln=n
*       llda=lda
*       lincx=incx
*       lincy=incy
*      call dgemtx(lm,ln,alpha,a,llda,x,lincx,y,lincy)
*@else
*      integer m,n,lda,incx,incy
*      call dgemtx(m,n,alpha,a,lda,x,incx,y,incy)
*@endif
*      return
*      end
*@endif
        subroutine dsymm_wr( SIDE, UPLO, M, N, ALPHA, A, LDA, B, LDB,
     $                   BETA, C, LDC )
        implicit none
        character*(*) side,uplo
        integer m,n,lda,ldb,ldc
        real*8 A(*), B(*), C(*), alpha,beta
*@if defined ( int64) && defined ( noblas64 )
*        integer*4 lm,ln,llda,lldb,lldc
*        lm=m
*        ln=n
*        llda=lda
*        lldb=ldb
*        lldc=ldc
*        call dsymm(side,uplo,lm,ln,alpha,a,llda,b,lldb,beta,c,lldc)
*@else
        call dsymm(side,uplo,m,n,alpha,a,lda,b,ldb,beta,c,ldc)
*@endif
        return
        end
 
         subroutine dpotrf_wr( UPLO, N, A, LDA, INFO )
         implicit none
         character*(*) uplo
         integer n,lda,info
         real*8 A(*) 
*@if defined ( int64) && defined ( noblas64 )
*         integer*4 ln,llda,linfo
*          ln=n
*          llda=lda
*          call dpotrf(uplo,ln,a,llda,linfo)
*          info=linfo
*@else
         call dpotrf(uplo,n,a,lda,info)
*@endif
         return
         end


         subroutine dpotri_wr( UPLO, N, A, LDA, INFO )
         implicit none
         character*(*) uplo
         integer n,lda,info
         real*8 A(*)
*@if defined ( int64) && defined ( noblas64 )
*         integer*4 ln,llda,linfo
*          ln=n
*          llda=lda
*          call dpotri(uplo,ln,a,llda,linfo)
*          info=linfo
*@else
         call dpotri(uplo,n,a,lda,info)
*@endif
         return
         end

      SUBROUTINE DSYEVD_WR( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK,IWORK,
     .                   LIWORK, INFO )
      CHARACTER          JOBZ, UPLO
      real*8   A( LDA, * ), W( * ), WORK( * )
      INTEGER            INFO, LDA, LIWORK, LWORK, N
      INTEGER            IWORK( * )
*@if defined (int64) && defined (noblas64) 
*      INTEGER*4          LINFO,LLDA,LLIWORK,LLWORK,LN
*      LLDA=LDA
*      LLIWORK=LIWORK
*      LLWORK=LWORK
*      LN = LN
*      LINFO=INFO
**
*      CALL DSYEVD(JOBZ,UPLO,LN,A,LLDA,W,WORK,LLWORK,IWORK,LLIWORK,LINFO)
**
*      INFO=LINFO
*@else
      CALL DSYEVD(JOBZ,UPLO,N,A,LDA,W,WORK,LWORK,IWORK,LIWORK,INFO)
*@endif
      return
      end


      SUBROUTINE DSYEVX_WR( JOBZ, RANGE, UPLO, N, A, LDA, VL, VU, IL,IU,
     .                   ABSTOL, M, W, Z, LDZ, WORK, LWORK, IWORK,
     .                   IFAIL, INFO )
*
      CHARACTER          JOBZ, RANGE, UPLO
      INTEGER            IL, INFO, IU, LDA, LDZ, LWORK, M, N
      real*8   ABSTOL, VL, VU
      INTEGER            IFAIL( * ), IWORK( * )
      real*8   A( LDA, * ), W( * ), WORK( * ), Z( LDZ, * )
*@if defined (int64) && defined (noblas64)
*      INTEGER*4   LIL, LINFO, LIU, LLDA, LLDZ, LLWORK, LM, LN
*      INTEGER*4, allocatable :: LIFAIL(:)
*      
*      LIL=IL
*      LINFO=INFO
*      LIU=IU
*      LLDA=LDA
*      LLDZ=LDZ
*      LLWORK=LWORK
*      LM=M
*      LN=N
*      allocate LIFAIL(LN)
*      CALL DSYEVX( JOBZ, RANGE, UPLO, LN, A, LLDA, VL, VU, LIL, LIU,
*     .                   ABSTOL, LM, W, Z, LLDZ, WORK, LLWORK, IWORK,
*     .                   LIFAIL, LINFO )
*      IF (LINFO.gt.0) then
*       DO I=1,LN
*         IFAIL(I)=LIFAIL(I)
*       ENDDO
*      ENDIF
*      INFO=LINFO
*@else
      CALL DSYEVX( JOBZ, RANGE, UPLO, N, A, LDA, VL, VU, IL, IU,
     .                   ABSTOL, M, W, Z, LDZ, WORK, LWORK, IWORK,
     .                   IFAIL, INFO )
*@endif
       return
       end





       

       



