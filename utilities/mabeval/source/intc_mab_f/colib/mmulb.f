c
c  main program for testing matrix multipliciation speeds
c  and tuning the colib library routines.
c
      implicit none 
c
      integer nlist, n
c
      integer   ndim
      parameter(ndim=2000*2000)
c
      real*8 a(ndim), b(ndim), c(ndim), right(ndim)
c
      nlist = 6
c
      call header( nlist )
c
      do 10 n = 1, 7
c        call comp( a, b, c, right, n,10 )
10    continue
c
      do 100 n = 8, 128, 8
         call comp( a, b, c, right, n,100)
100   continue
c
      do 200 n = 200, 500, 100
         call comp( a, b, c, right, n,10)
200   continue

      do 300 n = 750, 2000, 250
         call comp( a, b, c, right, n,1 )
300   continue


c
9999  stop
      end
deck comp
      subroutine comp( a, b, c, right, n ,repeatm)
c
c  compare matrix multiplication execution times for various
c  algorithms and library routines.
c
      implicit none
c
c     for better timing accuracy repeat each measurement repeatm times
      integer   irepeat,repeatm
c
      integer n, iunit
      real*8 a(n,n), b(n,n), c(n,n), right(n,n)
c
      integer i, j, itype, k, nerror
      real*8 cij, bkj, dcpu, diffmx
      real*8 dnflop, cpu0(1), cpu1(1)
c
      real*8   ddot_wr
      external ddot_wr
c
      real*8    zero,     one,     e6,     floor
      parameter(zero=0d0, one=1d0, e6=1d6, floor=1d-5)
c
      integer nlist
      save nlist
c
c     # ntypes = number of types of products to compute and measure.
c     #          the first 9 are colib matrix-matrix product routines.
c     #          The next 4 are colib matrix-vector product routines.
c     #          additional calls should be added as necessary using
c     #          vendor-specific routines.
c
      integer    ncolib
      parameter( ncolib=9+4+4 )
c
      integer    ntypes
      parameter (ntypes=17)
      character*12 clabs(ntypes)
      real*8 mflops(ntypes)
c
      data clabs /
     & 'fortran-dot', 'fortran-axpy', 'ddot()',    'daxpy()',
     & 'mxm()',       'mxma()',       'mtxm()',    'mxmt()',    
     & 'mtxmt()',     'm*v(mxma)',    'm*v(mxv)',  'm*v(mtxv)',  
     & 'm*v(mxva)'
     & , 'dgemm()',   'm*v(dgemm)', 'm*v(dgemv)', 'matmul()'
     & /
c
c     # initialize a and b to have non-zero, non-trivial entries.
c     # a(*) and b(*) are symmetric to facilitate later timings
c     # involving various transpositions.
c
      do 20 j = 1, n
         do 10 i = 1, n
            b(i,j) = (i+j)
            a(i,j) = (i+j)**2
10       continue
20    continue
c
      itype = 0
c
c     # number of floating point operations, scaled to measure mflops.
c
      dnflop = 2*n*n*n*repeatm
      dnflop = dnflop/e6
c
c     # matrix multiplication using fortran dot products.
c     # the right result is saved in right(*,*)
c
      call getima( 1, cpu0 )
      do irepeat=1,repeatm
      do 50 i =1, n
         do 40 j =1, n
            cij = zero
            do 30 k = 1, n
               cij = cij + a(i,k) * b(k,j)
30          continue
            right(i,j) = cij
40       continue
50    continue
      enddo 
      call getima( 1, cpu1 )
      dcpu   = cpu1(1)-cpu0(1)
      itype  = itype + 1
      mflops(itype) = dnflop/max(dcpu,floor)
c
c     # matrix multiplication with fortran daxpy.
c
      call getima( 1, cpu0 )
      do irepeat=1,repeatm
      call wzero( n*n, c, 1 )
      do 80 j = 1, n
         do 70 k = 1, n
            bkj = b(k,j)
            do 60 i = 1, n
               c(i,j) = c(i,j) + a(i,k) * bkj
60          continue
70       continue
80    continue
      enddo
      call getima( 1, cpu1 )
      dcpu   = cpu1(1)-cpu0(1)
      itype  = itype + 1
      mflops(itype) = dnflop/max(dcpu,floor)
      call check( n, c, right, diffmx, nerror )
      if ( nerror .ne. 0 ) write(nlist,6060) nerror, itype, diffmx
c
c     # matix multiplication with ddot.
c
      call getima( 1, cpu0 )
      do irepeat=1,repeatm
      do 120 i = 1, n
         do 110 j = 1, n
            c(i,j) = ddot_wr( n, a(i,1), n,   b(1,j), int(1,8) )
110      continue
120   continue
      enddo
      call getima( 1, cpu1 )
      dcpu   = cpu1(1)-cpu0(1)
      itype  = itype + 1
      mflops(itype) = dnflop/max(dcpu,floor)
      call check( n, c, right, diffmx, nerror )
      if ( nerror .ne. 0 ) write(nlist,6060) nerror, itype, diffmx
c
c     # matrix multiplication with daxpy.
c
      call getima( 1, cpu0 )
      do irepeat=1,repeatm
      call wzero( n*n, c, 1 )
      do 140 j = 1, n
         do 130 k = 1, n
            call daxpy_wr( n, b(k,j),   a(1,k), 1,    c(1,j), 1 )
130      continue
140   continue
      enddo 
      call getima( 1, cpu1 )
      dcpu   = cpu1(1)-cpu0(1)
      itype  = itype + 1
      mflops(itype) = dnflop/max(dcpu,floor)
      call check( n, c, right, diffmx, nerror )
      if ( nerror .ne. 0 ) write(nlist,6060) nerror, itype, diffmx
c
c     # matrix multiplication with mxm() subroutine call.
c
      call getima( 1, cpu0 )
      do irepeat=1,repeatm
      call mxm( a, n,  b, n,  c, n )
      enddo 
      call getima( 1, cpu1 )
      dcpu   = cpu1(1)-cpu0(1)
      itype  = itype + 1
      mflops(itype) = dnflop/max(dcpu,floor)
      call check( n, c, right, diffmx, nerror )
      if ( nerror .ne. 0 ) write(nlist,6060) nerror, itype, diffmx
c
c     # matrix multiplication with mxma() subroutine call.
c
      call getima( 1, cpu0 )
      do irepeat=1,repeatm
      call mxma( a,1,n,  b,1,n,  c,1,n,  n,n,n )
      enddo 
      call getima( 1, cpu1 )
      dcpu   = cpu1(1)-cpu0(1)
      itype  = itype + 1
      mflops(itype) = dnflop/max(dcpu,floor)
      call check( n, c, right, diffmx, nerror )
      if ( nerror .ne. 0 ) write(nlist,6060) nerror, itype, diffmx
c
c     # matrix(transpose)*matrix subroutine call.
c
      call getima( 1, cpu0 )
      do irepeat=1,repeatm
      call mtxm( a, n,    b, n,    c, n )
      enddo 
      call getima( 1, cpu1 )
      dcpu   = cpu1(1)-cpu0(1)
      itype  = itype + 1
      mflops(itype) = dnflop/max(dcpu,floor)
      call check( n, c, right, diffmx, nerror )
      if ( nerror .ne. 0 ) write(nlist,6060) nerror, itype, diffmx
c
c     # matrix*matrix(transpose) subroutine call.
c
      call getima( 1, cpu0 )
      do irepeat=1,repeatm
      call mxmt( a, n,    b, n,    c, n )
      enddo 
      call getima( 1, cpu1 )
      dcpu   = cpu1(1)-cpu0(1)
      itype  = itype + 1
      mflops(itype) = dnflop/max(dcpu,floor)
      call check( n, c, right, diffmx, nerror )
      if ( nerror .ne. 0 ) write(nlist,6060) nerror, itype, diffmx
c
c     # matrix(transpose)*matrix(transpose) subroutine call.
c
      call getima( 1, cpu0 )
      do irepeat=1,repeatm
      call mtxmt( a, n,    b, n,    c, n )
      enddo 
      call getima( 1, cpu1 )
      dcpu   = cpu1(1)-cpu0(1)
      itype  = itype + 1
      mflops(itype) = dnflop/max(dcpu,floor)
      call check( n, c, right, diffmx, nerror )
      if ( nerror .ne. 0 ) write(nlist,6060) nerror, itype, diffmx
c
c     # matrix multiplication with series of matrix-vector products.
c
c     # mxma() calls.
c
      call getima( 1, cpu0 )
      do irepeat=1,repeatm
      do 150 i = 1, n
         call mxma(a,1,n, b(1,i),1,1, c(1,i),1,1, n,n,1)
150   continue
      enddo 
      call getima( 1, cpu1 )
      dcpu   = cpu1(1)-cpu0(1)
      itype  = itype + 1
      mflops(itype) = dnflop/max(dcpu,floor)
      call check( n, c, right, diffmx, nerror )
      if ( nerror .ne. 0 ) write(nlist,6060) nerror, itype, diffmx
c
c     # mxv() calls.
c
      call getima( 1, cpu0 )
      do irepeat=1,repeatm
      do 160 i = 1, n
         call mxv(a, n, b(1,i), n, c(1,i) )
160   continue
       enddo 
      call getima( 1, cpu1 )
      dcpu   = cpu1(1)-cpu0(1)
      itype  = itype + 1
      mflops(itype) = dnflop/max(dcpu,floor)
      call check( n, c, right, diffmx, nerror )
      if ( nerror .ne. 0 ) write(nlist,6060) nerror, itype, diffmx
c
c     # mtxv() calls.
c
      call getima( 1, cpu0 )
      do irepeat=1,repeatm
      do 170 i = 1, n
         call mtxv( a,n,    b(1,i), n,    c(1,i) )
170   continue
      enddo 
      call getima( 1, cpu1 )
      dcpu   = cpu1(1)-cpu0(1)
      itype  = itype + 1
      mflops(itype) = dnflop/max(dcpu,floor)
      call check( n, c, right, diffmx, nerror )
      if ( nerror .ne. 0 ) write(nlist,6060) nerror, itype, diffmx
c
c     # mxva() calls.
c
      call getima( 1, cpu0 )
      do irepeat=1,repeatm
      do 180 i = 1, n
         call mxva( a,1,n, b(1,i),1, c(1,i),1, n,n)
180   continue
      enddo 
      call getima( 1, cpu1 )
      dcpu   = cpu1(1)-cpu0(1)
      itype  = itype + 1
      mflops(itype) = dnflop/max(dcpu,floor)
      call check( n, c, right, diffmx, nerror )
      if ( nerror .ne. 0 ) write(nlist,6060) nerror, itype, diffmx
c
c
c     # dgemm() call.
c
      call getima( 1, cpu0 )
      do irepeat=1,repeatm
      call dgemm_wr( 'n', 'n', n,n,n, one, a,n, b,n, zero, c,n )
      enddo 
      call getima( 1, cpu1 )
      dcpu   = cpu1(1)-cpu0(1)
      itype  = itype + 1
      mflops(itype) = dnflop/max(dcpu,floor)
      call check( n, c, right, diffmx, nerror )
      if ( nerror .ne. 0 ) write(nlist,6060) nerror, itype, diffmx
c
c     # series of dgemm() calls.
c
      call getima( 1, cpu0 )
      do irepeat=1,repeatm
      do 200 i = 1, n
         call dgemm_wr( 'n', 'n', n,1,n,
     &    one, a,n, b(1,i),n, zero, c(1,i),n )
200   continue
      enddo 
      call getima( 1, cpu1 )
      dcpu   = cpu1(1)-cpu0(1)
      itype  = itype + 1
      mflops(itype) = dnflop/max(dcpu,floor)
      call check( n, c, right, diffmx, nerror )
      if ( nerror .ne. 0 ) write(nlist,6060) nerror, itype, diffmx
c
c     # series of dgemv() calls.
c
      call getima( 1, cpu0 )
      do irepeat=1,repeatm
      do 210 i = 1, n
         call dgemv_wr( 'n', n,n, one, a,n, b(1,i),1, zero, c(1,i),1 )
210   continue
      enddo 
      call getima( 1, cpu1 )
      dcpu   = cpu1(1)-cpu0(1)
      itype  = itype + 1
      mflops(itype) = dnflop/max(dcpu,floor)
      call check( n, c, right, diffmx, nerror )
      if ( nerror .ne. 0 ) write(nlist,6060) nerror, itype, diffmx
c
c     # intrinsic matmul() call.
c
      call getima( 1, cpu0 )
      do irepeat=1,repeatm
      c = matmul( a, b )
      enddo 
      call getima( 1, cpu1 )
      dcpu   = cpu1(1)-cpu0(1)
      itype  = itype + 1
      mflops(itype) = dnflop/max(dcpu,floor)
      call check( n, c, right, diffmx, nerror )
      if ( nerror .ne. 0 ) write(nlist,6060) nerror, itype, diffmx
c
      write(nlist,6040) n, (mflops(i), i=1,itype)
c
      return
c
      entry header( iunit )
c
c     # print out a header to label the columns.
c
      nlist = iunit
c
      write(nlist,6010) (i, clabs(i), i=1,ntypes)
      write(nlist,6020) ' n ', (i, i=1,ntypes)
      write(nlist,6030) '---', (' -----', i=1,ntypes)
c
      return
c
6010  format(' mmulb product type legend:'/ (i3,': ',a) )
6020  format(/20x,'mflops for each product type:'/1x,a4,20(i6,2x))
6030  format(1x,a4,20a8)
6040  format(1x,i3,20f8.1)
6060  format(/' *** nerror=',i8,' itype=',i3,' diffmx=',1pe12.4,' ***')
c
      end
deck check
      subroutine check( n, c, right, diffmx, nerror )
c
c  compare c(*) to right(*).
c  return: diffmx = max( abs( c(i,j)-right(i,j)))
c          nerror = number of elements that differ by more than
c                   2 digits.
      implicit none 
c
      integer n, nerror
      real*8 c(n,n), right(n,n), diffmx
c
      integer i, j
      real*8 diff, term
c
      real*8     zero,     one,     r100
      parameter( zero=0d0, one=1d0, r100=100d0 )
c
      diffmx = zero
      nerror = 0
      do 20 j = 1, n
         do 10 i = 1, n
            diff = c(i,j) -right(i,j)
            if ( diff .ne. zero ) then
               diffmx = max( diffmx, abs(diff) )
               term = abs(c(i,j)) + abs(right(i,j))
               if ( (one + diff/(term*r100)) .ne. one ) then
                  nerror = nerror + 1
               endif
            endif
10       continue
20    continue
c
      return
      end


