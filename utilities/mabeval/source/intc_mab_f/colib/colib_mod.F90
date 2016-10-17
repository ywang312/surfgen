    MODULE colib 
#ifdef INT64
      integer,parameter :: kint=8
#else
      integer,parameter :: kint=4
#endif


      CONTAINS

          subroutine little_endian(ll)
          implicit none
          integer*2 longw
          integer*1 shortw(2)
          equivalence (longw,shortw)
          logical ll
          longw=1
          ll=shortw(1).eq.1
          write(6,*) 'little_endian=',ll
          return
          if (shortw(1).eq.1) then
              write(6,*) 'little endian'
          elseif (shortw(2).eq.1) then
              write(6,*) 'big endian'
          else
              stop 'failed check_endian'
          endif
          return
          end subroutine little_endian



      subroutine plab16( p, u, nuw )
      implicit none
!     maxbuf = must be a multiple of 4 
      integer, parameter ::  maxbuf = 128
      integer nuw
      real*8 p(nuw),p1
      integer u(nuw)
      integer*2 u2(maxbuf),u1
      integer iup,plen,forbyt,i,imax
!
!  pack integral labels from u(*) into p(*,*).
!
!  p(*) = packed array (working precision in the calling program).
!  u(*) = unpacked array.  u( 1 : ((nuw+3)/4)*4 ) are referenced.
!  nuw  = number of unpacked integral labels.

        if (nuw.le.0) return
        if (nuw.le.maxbuf) then
          iup=(nuw-1)/4+1 
!    maps to the range of +-32767
          u2(1:nuw)=u(1:nuw)
100    format(a,20i4)
!    IBM AIX/XLF does not work properly with size argument??
          p(1:iup)=transfer(u2(1:nuw),p)
          return
        else
          imax=0
          do i=1,((nuw-1)/maxbuf)*maxbuf,maxbuf
            u2(1:maxbuf)=u(i:i+maxbuf-1) 
            p(i/4+1:(i+maxbuf)/4)=transfer(u2(1:maxbuf),p)
            imax=i
          enddo
           imax=imax+maxbuf
           iup=(nuw-1)/4+1 
            u2=u(imax:nuw) 
            p(imax/4+1:iup)=transfer(u2(1:nuw-imax+1),p)
          return
         endif
         return
        end subroutine plab16
         
      subroutine ulab16( p, u, nuw )
      implicit none
      integer, parameter ::  maxbuf = 128
      integer, parameter ::  mask = z'ffff'
      integer nuw
      real*8 p(nuw),p1(maxbuf/4)
      integer u(nuw)
      integer*2 u2(maxbuf),u1
      integer iup,plen,forbyt,i,imax,j
!
!  unpack integral labels from p(*) into u(*,*).
!
!  p(*) = packed array (working precision in the calling program).
!  u(*) = unpacked array.  u( 1 : ((nuw+3)/4)*4 ) are referenced.
!  nuw  = number of unpacked integral labels.

        if (nuw.le.0) return
        iup=(nuw-1)/4+1 
        if (nuw.le.maxbuf) then
          u2(1:nuw)=transfer(p(1:iup),u2(1:nuw))
!    maps to the range of +-32767 so use iand
!    little/big endian of no consequence
!  g95 does not allow the usage of arrays
           u(1:nuw)=iand(mask,int(u2(1:nuw),kint))
          return
        else
          imax=0
          do i=1,((nuw-1)/maxbuf)*maxbuf,maxbuf
            u2(1:maxbuf)=transfer(p(i/4+1:i/4+maxbuf/4),u2)
            u(i:i+maxbuf-1)=iand(mask,int(u2(1:maxbuf),kint))
            imax=i
          enddo
           imax=imax+maxbuf
            u2(1:nuw-imax+1)=transfer(p(imax/4+1:iup),u2)
            u(imax:nuw)=iand(mask,int(u2(1:nuw-imax+1),kint))
          return
         endif
         return
        end subroutine ulab16



      subroutine plab32( p, u, nuw )
      implicit none
!     maxbuf = must be a multiple of 4
      integer, parameter ::  maxbuf = 16
      integer nuw
      real*8 p(nuw),p1
      integer u(nuw)
      integer*4 u2(maxbuf),u1
      integer iup,plen,forbyt,i,imax
!
!  pack integral labels from u(*) into p(*,*).
!
!  p(*) = packed array (working precision in the calling program).
!  u(*) = unpacked array.  u( 1 : ((nuw+3)/4)*4 ) are referenced.
!  nuw  = number of unpacked integral labels.

        if (nuw.le.0) return 
        if (nuw.le.maxbuf) then
          iup=(nuw-1)/2+1
!    maps to the range of +-(2**31-1)
          u2(1:nuw)=u(1:nuw)
100    format(a,20i4)
          p(1:iup)=transfer(u2(1:nuw),p)
          return
        else
          imax=0
          do i=1,((nuw-1)/maxbuf)*maxbuf,maxbuf
            u2(1:maxbuf)=u(i:i+maxbuf-1)
            p(i/2+1:(i+maxbuf)/2)=transfer(u2(1:maxbuf),p)
            imax=i
          enddo
           imax=imax+maxbuf
           iup=(nuw-1)/2+1 
            u2=u(imax:nuw) 
            p(imax/2+1:iup)=transfer(u2(1:nuw-imax+1),p)
          return
         endif
         return
        end subroutine plab32

      subroutine ulab32( p, u, nuw )
      implicit none
      integer, parameter ::  maxbuf = 16
      integer, parameter ::  mask = z'ffffffff'
      integer nuw
      real*8 p(nuw),p1
      integer u(nuw)
      integer*4 u2(maxbuf),u1
      integer iup,plen,forbyt,i,imax,j
!
!  unpack integral labels from p(*) into u(*,*).
!
!  p(*) = packed array (working precision in the calling program).
!  u(*) = unpacked array.  u( 1 : ((nuw+3)/4)*4 ) are referenced.
!  nuw  = number of unpacked integral labels.

        if (nuw.le.0) return
        iup=(nuw-1)/2+1 
        if (nuw.le.maxbuf) then
          u2(1:nuw)=transfer(p(1:iup),u2)
!    maps to the range of +-(2**31-1) so use iand
!    little/big endian of no consequence
          u(1:nuw)=iand(mask,int(u2(1:nuw),kint))
          return
        else
          imax=0
          do i=1,((nuw-1)/maxbuf)*maxbuf,maxbuf
            u2(1:maxbuf)=transfer(p(i/2+1:i/2+maxbuf/2),u2)
            u(i:i+maxbuf)=iand(mask,int(u2(1:maxbuf),kint))
            imax=i
          enddo
           imax=imax+maxbuf
            u2(1:nuw-imax+1)=transfer(p(imax/2+1:iup),u2)
            u(imax:nuw)=iand(mask,int(u2(1:nuw-imax+1),kint))
          return
         endif
         return
        end subroutine ulab32


      integer function forbyt( i )
!
!  compute the number of working precision words required to hold
!  an (unpacked) integer array.
!
      integer i
!
       if (bit_size(I).eq.32) then
          forbyt = (i+1)/2
       elseif (bit_size(I).eq.64) then
          forbyt= i
!      else
!       call bummer('invalid bit_size=',bit_size(I),2)
       ENDIF
      return
      end function forbyt


 
      subroutine plab8( p, u, nuw )
      implicit none
!     maxbuf = must be a multiple of 8
      integer, parameter ::  maxbuf = 16
      integer nuw
      real*8 p(nuw),p1
      integer u(nuw)
      integer*1 u2(maxbuf),u1
      integer iup,plen,forbyt,i,imax
!
!  pack integral labels from u(*) into p(*,*).
!
!  p(*) = packed array (working precision in the calling program).
!  u(*) = unpacked array.  u( 1 : ((nuw+3)/4)*4 ) are referenced.
!  nuw  = number of unpacked integral labels.

        if (nuw.le.0) return
        if (nuw.le.maxbuf) then
          iup=(nuw-1)/8+1
!
!     this conversion removes the little/big endian problem
!     however, integer*1 range  -127:127  (one bit for sign)
!     so decompressing must be done by bit operations and
!     little/big endian has to be considered
!
          u2(1:nuw)=u(1:nuw)
100    format(a,20i4)
          p(1:iup)=transfer(u2(1:nuw),p)
!         write(6,'10(z16)') p(1:5)
          return
        else
          imax=0
          do i=1,((nuw-1)/maxbuf)*maxbuf,maxbuf
            u2(1:maxbuf)=u(i:i+maxbuf-1)
            p(i/8+1:(i+maxbuf)/8)=transfer(u2(1:maxbuf),p)
            imax=i
          enddo
           imax=imax+maxbuf
           iup=(nuw-1)/8+1
            u2=u(imax:nuw)
            p(imax/8+1:iup)=transfer(u2(1:nuw-imax+1),p)
          return
         endif
         return
        end subroutine plab8

      subroutine ulab8( p, u, nuw )
      implicit none
      integer, parameter ::  maxbuf = 16
      integer, parameter ::  mask = z'ff'
      integer nuw
      real*8 p(nuw),p1
      integer u(nuw)
      integer*1 u2(maxbuf),u1
      integer iup,plen,forbyt,i,imax,j
      logical first,ll
      integer endianshift
      data first /.true./
      save first,endianshift

!  on both little_endian and big_endian machines this
!  code works without additional shift
#ifdef obsolete
       if (first) then
         call little_endian(ll)
         if(ll) then
          endianshift=0
         else
          endianshift=0
         endif
        first=.false.
       endif  
#endif
  
!
!  unpack integral labels from p(*) into u(*,*).
!
!  p(*) = packed array (working precision in the calling program).
!  u(*) = unpacked array.  u( 1 : ((nuw+3)/4)*4 ) are referenced.
!  nuw  = number of unpacked integral labels.

        if (nuw.le.0) return
        iup=(nuw-1)/8+1 
        if (nuw.le.maxbuf) then
!         write(6,'a,10(z16)') 'ulab p:',p(1:5)
          u2(1:nuw)=transfer(p(1:iup),u2)
!
!       bit operations restore range 1:255
!
#ifdef obsolte
          do j=1,nuw
           u(j)=ishft(iand(mask,int(u2(j),kint)),endianshift)
          enddo 
#else
           u(1:nuw)=(iand(mask,int(u2(1:nuw),kint)))
#endif
          return
        else
          imax=0
          do i=1,((nuw-1)/maxbuf)*maxbuf,maxbuf
            u2(1:maxbuf)=transfer(p(i/8+1:i/8+maxbuf/8),u2)
#ifdef obsolete
            do j=i,i+maxbuf-1
             u(j)=ishft(iand(mask,int(u2(j-i+1),kint)),endianshift)
            enddo
#else
           u(i:i+maxbuf-1)=(iand(mask,int(u2(1:maxbuf),kint)))
#endif
            imax=i
          enddo
           imax=imax+maxbuf
            u2(1:nuw-imax+1)=transfer(p(imax/8+1:iup),u2)
#ifdef obsolete
            do j=imax,nuw
             u(j)=ishft(iand(mask,int(u2(j-imax+1),kint)),endianshift)
            enddo
#else
           u(imax:nuw)=(iand(mask,int(u2(1:nuw-imax+1),kint)))
#endif
          return
         endif
         return
        end subroutine ulab8



      subroutine plab1( p, u, nuw )
      implicit none
      integer, parameter ::  maxbuf = 64
      integer, parameter ::  mask = 1 
      integer nuw
      real*8 p(nuw),pp
      integer*8 p1
      integer u(nuw)
      integer*1 u2(maxbuf),u1
      integer icur,iup,plen,forbyt,i,imax,j


        if (nuw.le.0) return
        if (nuw.le.maxbuf) then
          iup=(nuw-1)/64+1
          u2(1:nuw)=u(1:nuw)
100    format(a,20i4)
          p1=0
          do i=0,nuw-1
            p1=ior(p1, ishft(iand(mask,u(i+1)),i))
          enddo 
          p(1)=transfer(p1,pp)
          return
        else
          imax=0
          icur=1
          do j=1,((nuw-1)/maxbuf)*maxbuf,maxbuf
            u2(1:maxbuf)=u(j:j+maxbuf-1)
             p1=0
            do i=0,maxbuf-1
            p1=ior(p1, ishft(iand(mask,int(u(j+i),kint)),i))
            enddo 
            p(icur)=transfer(p1,pp)
            imax=j
            icur=icur+1
          enddo
           imax=imax+maxbuf
           iup=(nuw-1)/64+1
            u2=u(imax:nuw)
             p1=0
            do i=0,nuw-imax
            p1=ior(p1, ishft(iand(mask,int(u(imax+i),kint)),i))
            enddo 
            p(icur)=transfer(p1,pp)
          return
         endif
         return
        end subroutine plab1

      subroutine ulab1( p, u, nuw )
      implicit none
      integer, parameter ::  maxbuf = 64
      integer nuw
      real*8 p(nuw)
      integer*8 p1
      integer u(nuw)
      integer*1 u2(maxbuf),u1
      integer, parameter :: ione=1
      integer icur,iup,plen,forbyt,i,imax,j

        if (nuw.le.0) return
        if (nuw.le.maxbuf) then
          p1=transfer(p(1),p1)
          do i=0,nuw-1
            u1=iand(ione,ishft(p1,-i))
!   the copy step avoids little/big endian
!   dependencies by implicit transformation
            u(i+1)=u1
          enddo
          return
        else
          imax=0
          icur=1
          do j=1,((nuw-1)/maxbuf)*maxbuf,maxbuf
            p1=transfer(p(icur),p1)
            do i=0,maxbuf-1
              u1=iand(ione,ishft(p1,-i))
              u(imax+i+1)=u1
            enddo 
            imax=j-1
            icur=icur+1
          enddo
           imax=imax+maxbuf
           iup=(nuw-1)/64+1
            p1=transfer(p(icur),p1)
            do i=0,nuw-imax
              u1=iand(ione,ishft(p1,-i))
              u(imax+i+1)=u1
            enddo
          return
         endif
         return
        end subroutine ulab1

   end module colib


          
 
