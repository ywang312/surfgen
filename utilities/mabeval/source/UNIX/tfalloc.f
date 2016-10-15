c     # test the falloc() memory allocation routine.
c     # see also $COLUMBUS/source/replace_me/main.f
      real*8 work(1)
      integer nelem, clean, addr, offset, ifirst
c
      nelem = 100000
      clean = 0
      write(*,*) 'input nelem, clean'
      read(*,*) nelem, clean
c
      call falloc( nelem, 8, clean, work, addr, offset )
      ifirst = offset + 1
c
      if ( addr .ne. 0 ) then
c        # work(ifirst : ifirst+nelem-1) should be useable.
         call driver( work(ifirst), nelem, addr, ifirst )
      else
         write(*,*) 'falloc() error, addr=0'
      endif
c
      stop
      end
      subroutine driver( a, lcore, mem1, ifirst )
      implicit integer(a-z)
      integer lcore, mem1, ifirst
      real*8 a(lcore)
c
      write(*,6010) lcore, mem1, ifirst
6010  format(' workspace allocation information: lcore=',i10,
     & ' mem1=',i10,' ifirst=',i10)
c
c     # write to the first and last elements of a(*).
c
      a(1)     = (1)
      a(lcore) = (lcore)
      write(*,6020) a(1), a(lcore)
6020  format(' a(1)=',f4.1,5x,'a(lcore)=',f12.1)
c
      return
      end
