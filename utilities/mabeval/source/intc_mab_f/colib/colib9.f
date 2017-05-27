
ctm  stems from cisrt2.f

      subroutine code3x (
     & ibuffer, intblk, stn121, stn2,
     & stn3,   nintkl, ksym,   k1s,
     & k1f,    l1,     iscr )
c
c  pack the 3-external header info into two working precision words.
c
c  03-jul-01 modified to allow large basis sets (>170). -rls
c
      implicit none 
c     # dummy:
      integer intblk, stn121, stn2, stn3, nintkl,
     & ksym, k1s, k1f, l1, iscr
      integer i,j, stn3hi,stn3lo
      integer*8 intblk8,stn1218,stn28,stn38,nintkl8,
     .  ksym8,k1s8,k1f8,l18,iscr8,stn3hi8,stn3lo8

c     # bummer error types.
      integer   wrnerr,  nfterr,  faterr
      parameter(wrnerr=0,nfterr=1,faterr=2)

c     # to simplify matters and since header packing is not
c     # time critical at all(tm)  

c     # maximum number of basis functions 1023
      integer maxbfn
      parameter (maxbfn=2**10-1)
c     #maximum buffer length  262143 
      integer maxsiz
      parameter (maxsiz=2**18-1)
      integer*8 m46,m28,m10,m56,m38,m33,m23,m13,m3
      parameter (m46=46,m28=28,m10=10,m56=56,m38=38,m33=33,
     .           m23=23,m13=13,m3=3)

      integer*8 ibuffer(2)



c     intblk           <   2**18 -1     18bit
c     stn121           <   2**18 -1     18bit
c     stn2             <   2**18 -1     18bit
c     stn3             <   2**18 -1     18bit
c     nintkl           <   2**18 -1     18bit
c     ksym             <   9             5bit
c     k1s              <   2**10 -1     10bit
c     k1f              <   2**10 -1     10bit
c     l1               <   2**10 -1     10bit
c     iscr             <   4             3bit
c-----------------------------------------------
c        total                          128 bit


      if ( intblk .gt. maxsiz .or. intblk .lt. 0 ) then
         call bummer(': intblk=', intblk, faterr)
      elseif ( stn121 .gt. maxsiz .or. stn121 .lt. 0 ) then
         call bummer(': stn121=', stn121, faterr)
      elseif ( stn2 .gt. maxsiz .or. stn2 .lt. 0 ) then
         call bummer(': stn2=', stn2, faterr)
      elseif ( stn3 .gt. maxsiz .or. stn3 .lt. 0 ) then
         call bummer(': stn3=', stn3, faterr)
      elseif ( nintkl .gt. maxsiz .or. nintkl .lt. 0 ) then
         call bummer(': nintkl=', nintkl, faterr)
      elseif ( ksym .gt. 8 .or. ksym .lt. 0 ) then
         call bummer(': ksym=', ksym, faterr)
      elseif ( k1s .gt. maxbfn .or. k1s .lt. 0 ) then
         call bummer(': k1s=', k1s, faterr)
      elseif ( k1f .gt. maxbfn .or. k1f .lt. 0 ) then
         call bummer(': k1f=', k1f, faterr)
      elseif ( l1 .gt. maxbfn .or. l1 .lt. 0 ) then
         call bummer(': l1=', l1, faterr)
      elseif ( iscr .gt. 3 .or. iscr .lt. 0 ) then
         call bummer(': iscr=', iscr, faterr)
      endif

    
      stn3hi=ishft(stn3,-8)
      stn3lo=iand(stn3,2**8-1)

*@ifdef int64
      ibuffer(1)=ior(ishft(intblk,m46), ior(ishft(stn121,m28),
     &          ior(ishft(stn2,  m10), stn3hi)))
      ibuffer(2)=ior(ishft(stn3lo,m56),ior(ishft(nintkl,m38),
     .          ior(ishft(ksym,m33), ior(ishft(k1s,m23),
     .          ior(ishft(k1f,m13), ior(ishft(l1,m3),iscr))))))
*@else
*C     intrinsics usually require explict type conversions
**
*      intblk8=intblk
*      stn1218=stn121
*      stn28=stn2
*      stn3lo8=stn3lo
*      stn3hi8=stn3hi
*      nintkl8=nintkl
*      ksym8=ksym
*      k1s8=k1s
*      k1f8=k1f
*      l18=l1
*      iscr8=iscr
*      ibuffer(1)=ior(ishft(intblk8,m46), 
*     .           ior(ishft(stn1218,m28),
*     &          ior(ishft(stn28, m10), stn3hi8)))
*      ibuffer(2)=ior(ishft(stn3lo8,m56),
*     .          ior(ishft(nintkl8,m38),
*     .          ior(ishft(ksym8,m33), 
*     .          ior(ishft(k1s8,m23),
*     .          ior(ishft(k1f8,m13), 
*     .          ior(ishft(l18,m3),iscr8))))))
*@endif
*@ifdef debug
*       call printbinary(ibuffer(1),'code3x ibuffer(1)')
*       call printbinary(ibuffer(2),'code3x ibuffer(2)')
*@endif
      return
      end

ctm code4x stems from cisrt2.f


      subroutine code4x(
     & ibuffer, strtnf, stn121, stn3,
     & nintkl, lsym,   ksym,   l1s,
     & l1f,    k1s,    k1f,    k1l1 )
       implicit none
c
c  pack the 4-external header info into three working precision words.
c
c     # to simplify matters and since header packing is not
c     # time critical at all(tm)

c     # dummy:
      integer strtnf, stn121, stn3, nintkl, lsym, ksym,
     & l1s, l1f, k1s, k1f, k1l1
      integer i,j
c
c     # bummer error types.
      integer   wrnerr,  nfterr,  faterr
      parameter(wrnerr=0,nfterr=1,faterr=2)


c     # maximum number of basis functions 1023
      integer maxbfn
      parameter (maxbfn=2**10-1)
c     #maximum buffer length  1048575
      integer maxsiz
      parameter (maxsiz=2**20-1)

      integer nintklhi,nintkllo
      integer*8 ibuffer(3)
      integer*8 m44,m24,m4,m48,m40,m20,m10
      parameter (m44=44,m24=24,m4=4,m48=48,m20=20,m10=10,
     .  m40=40)
      integer*8 strtnf8,stn1218,stn38,nintkl8,lsym8,ksym8,
     .  l1s8,l1f8,k1s8,k1f8,k1l18,nintklhi8,nintkllo8

c     strtnf           <   2**20 -1     20bit
c     stn121           <   2**20 -1     20bit
c     stn3             <   2**20 -1     20bit
c     nintkl           <   2**20 -1     20bit
c     lsym             <   9             4bit
c     Ksym             <   9             4bit
c     l1s              <   2**10 -1     10bit
c     l1f              <   2**10 -1     10bit
c     k1s              <   2**10 -1     10bit
c     k1f              <   2**10 -1     10bit
c     k1l1             <   2**20 -1     20bit
c-----------------------------------------------
c        total                          148 bit

c
c
      if ( strtnf .gt. maxsiz .or. strtnf .lt. 0 ) then
         call bummer('code4x: strtnf=', strtnf, faterr)
      elseif ( stn121 .gt. maxsiz .or. stn121 .lt. 0 ) then
         call bummer('code4x: stn121=', stn121, faterr)
      elseif ( stn3 .gt. maxsiz .or. stn3 .lt. 0 ) then
         call bummer('code4x: stn3=', stn3, faterr)
      elseif ( nintkl .gt. maxsiz .or. nintkl .lt. 0 ) then
         call bummer('code4x: nintkl=', nintkl, faterr)
      elseif ( lsym .gt. 8 .or. lsym .lt. 1 ) then
         call bummer('code4x: lsym=', lsym, faterr)
      elseif ( ksym .gt. 8 .or. ksym .lt. 1 ) then
         call bummer('code4x: ksym=', ksym, faterr)
      elseif ( l1s .gt. maxbfn .or. l1s .lt. 0 ) then
         call bummer('code4x: l1s=', l1s, faterr)
      elseif ( l1f .gt. maxbfn .or. l1f .lt. 0 ) then
         call bummer('code4x: l1f=', l1f, faterr)
      elseif ( k1s .gt. maxbfn .or. k1s .lt. 0 ) then
         call bummer('code4x: k1s=', k1s, faterr)
      elseif ( k1f .gt. maxbfn .or. k1f .lt. 0 ) then
         call bummer('code4x: k1f=', k1f, faterr)
      elseif ( k1l1 .gt. maxbfn*maxbfn .or. k1l1 .lt. 0 ) then
         call bummer('code4x: k1l1=', k1l1, faterr)
      endif

      nintklhi=iand(ishft(nintkl,-16),2**4-1)
      nintkllo=iand(nintkl,2**16-1)
*@ifdef int64
      ibuffer(1)=ior(ishft(strtnf,m44), ior(ishft(stn121,m24),
     &          ior(ishft(stn3, m4), nintklhi)))
      ibuffer(2)=ior(ishft(nintkllo,m48),ior(ishft(lsym,m44),
     .          ior(ishft(ksym,m40), ior(ishft(l1s,m20),
     .          ior(ishft(l1f,m10), k1s)))))
      ibuffer(3)=ior(ishft(k1f,m20),k1l1)
*@else
*      strtnf8=strtnf
*      stn1218=stn121
*      stn38=stn3
*      nintklhi8=nintklhi
*C     intrinsics usually require explict type conversions
*      ibuffer(1)=ior(ishft(strtnf8,m44), 
*     .           ior(ishft(stn1218,m24),
*     &           ior(ishft(stn38, m4), nintklhi8)))
*      nintkllo8=nintkllo
*       lsym8=lsym
*       ksym8=ksym
*       l1s8=l1s
*       l1f8=l1f
*       k1s8=k1s
*       k1f8=k1f
*       k1l18=k1l1
*      ibuffer(2)=ior(ishft(nintkllo8,m48),
*     .           ior(ishft(lsym8,m44),
*     .           ior(ishft(ksym8,m40), 
*     .           ior(ishft(l1s8,m20),
*     .           ior(ishft(l1f8,m10),k1s8)))))
*      ibuffer(3)=ior(ishft(k1f8,m20),k1l18)
*@endif
*@ifdef debug
*       call printbinary(ibuffer(1),'code4x ibuffer(1)')
*       call printbinary(ibuffer(2),'code4x ibuffer(2)')
*       call printbinary(ibuffer(3),'code4x ibuffer(3)')
*@endif
      return
      end



ctm dcod3x stems from ciudg7.f

      subroutine dcod3x(
     & ibuffer,  iblsym, ksym,   k1s,
     & k1f,    l1bl,   intblk, stn121,
     & stn3,   stn1,   stn2,   nintkl )
c
c  decodes integral labels prepared by cisrt.
c  see subroutine code3x () for the packing convention.
c
c
      implicit none 
c     # dummy:
      integer iblsym, ksym, k1s, k1f, l1bl, intblk, stn121, stn3,
     & stn1, stn2, nintkl
      integer*8 iblsym8,ksym8,k1s8,k1f8,l1bl8,intblk8,stn1218,
     . stn38,stn18,stn28,nintkl8
       integer i,j
       integer*8 ibuffer(2),iwork
c
c     # bummer error types.
       integer   wrnerr,  nfterr,  faterr
       parameter(wrnerr=0,nfterr=1,faterr=2)
       integer*8 m8,m3,m5,m10,m18
       parameter(m3=7,m5=2**5-1,m10=2**10-1, m18=2**18-1)
       parameter(m8=2**8-1)
       integer*8 mm8 
       parameter (mm8=8)
       integer*8 n33,n23,n13,n3,n56,n38,n46,n28,n10
       parameter (n33=-33,n23=-23,n13=-13,n3=-3,n56=-56,
     .          n38=-38,n46=-46,n28=-28,n10=-10)

c     intblk           <   2**18 -1     18bit
c     stn121           <   2**18 -1     18bit
c     stn2             <   2**18 -1     18bit
c     stn3             <   2**18 -1     18bit
c     nintkl           <   2**18 -1     18bit
c     ksym             <   9             5bit
c     k1s              <   2**10 -1     10bit
c     k1f              <   2**10 -1     10bit
c     l1               <   2**10 -1     10bit
c     iscr             <   4             3bit
c-----------------------------------------------
c        total                          128 bit

c


*@ifdef int64
       iwork=ibuffer(2)
       iblsym = iand(iwork,m3)
       ksym   = iand(ishft(iwork,n33),m5)
       k1s    = iand(ishft(iwork,n23),m10)
       k1f    = iand(ishft(iwork,n13),m10)
       l1bl   = iand(ishft(iwork,n3),m10)       
       stn3   = iand(ishft(iwork,n56),m8)
       nintkl   =iand(ishft(iwork,n38),m18)
       iwork=ibuffer(1)
       intblk = iand(ishft(iwork,n46),m18)
       stn121 = iand(ishft(iwork,n28),m18)
       stn3   = ior(ishft(iand(iwork,m10),mm8),stn3)
       if (iblsym.ne.1) then
           stn1 = stn121
           stn2 = iand(ishft(iwork,n10),m18)
       endif 
*@else
*       iwork=ibuffer(2)
*       iblsym8 = (iand(iwork,m3))
*       ksym8   = (iand(ishft(iwork,n33),m5))
*       k1s8    = (iand(ishft(iwork,n23),m10))
*       k1f8    = (iand(ishft(iwork,n13),m10))
*       l1bl8   = (iand(ishft(iwork,n3),m10))
*       stn38   = (iand(ishft(iwork,n56),m8))
*       nintkl8  =(iand(ishft(iwork,n38),m18))
*       iwork=ibuffer(1)
*       intblk8 = (iand(ishft(iwork,n46),m18))
*       stn1218 = (iand(ishft(iwork,n28),m18))
*       stn38   = (ior(ishft(iand(iwork,m10),mm8),stn38))
*       if (iblsym8.ne.1) then
*           stn1 = stn1218
*           stn28 = (iand(ishft(iwork,n10),m18))
*           stn2=stn28 
*       endif
*       iblsym=iblsym8
*       ksym=ksym8
*       k1s=k1s8 
*       k1f=k1f8
*       l1bl=l1bl8
*       stn3=stn38
*       nintkl=nintkl8
*       intblk=intblk8
*       stn121=stn1218
*       stn3=stn38
*@endif
*@ifdef debug
*       call printbinary(ibuffer(1),'dcod3x ibuffer(1)')
*       call printbinary(ibuffer(2),'dcod3x ibuffer(2)')
*@endif
      return
      end


ctm dcod4x stems from ciudg8.f

      subroutine dcod4x(
     & ibuffer,  intblk, stn121, stn3,
     & nintkl, lsym,   ksym,   l1s,
     & l1f,    k1s,    k1f,    k1l1 )
c
c  decodes integral labels prepared by cisrt
c  see subroutine code4x() for the packing convention.
c
c  03-jul-01 modified to allow large basis sets (>170). -rls
c
      implicit none
c     # dummy:
      integer intblk, stn121, stn3, nintkl, lsym, ksym, l1s, l1f,
     & k1s, k1f, k1l1
c
c     # bummer error types.
       integer*8 ibuffer(3),iwork
c
c     # bummer error types.
       integer   wrnerr,  nfterr,  faterr
       parameter(wrnerr=0,nfterr=1,faterr=2)
       integer*8 m4,m10,m16,m20
       parameter(m4=2**4-1,m10=2**10-1,m16=2**16-1,m20=2**20-1)
      integer*8 n16,n44,n24,n4,n48,n40,n20,n10
      parameter (n44=-44,n24=-24,n4=-4,n48=-48,n20=-20,n10=-10,
     .   n40=-40,n16=16)
      integer*8 intblk8,stn1218,stn38,nintkl8,lsym8,ksym8,l1s8,l1f8,
     .   k1s8,k1f8,k1l18




c     strtnf           <   2**20 -1     20bit
c     stn121           <   2**20 -1     20bit
c     stn3             <   2**20 -1     20bit
c     nintkl           <   2**20 -1     20bit
c     lsym             <   9             4bit
c     ksym             <   9             4bit
c     l1s              <   2**10 -1     10bit
c     l1f              <   2**10 -1     10bit
c     k1s              <   2**10 -1     10bit
c     k1f              <   2**10 -1     10bit
c     k1l1             <   2**20 -1     20bit
c-----------------------------------------------
c        total                          148 bit

*@ifdef int64
       iwork=ibuffer(1)
       intblk = iand(ishft(iwork,n44),m20)
       stn121 = iand(ishft(iwork,n24),m20)
       stn3   = iand(ishft(iwork,n4),m20)
       nintkl = ishft(iand(iwork,m4),n16)
       iwork=ibuffer(2)
       nintkl = ior(iand(ishft(iwork,n48),m16),nintkl)
       lsym   = iand(ishft(iwork,n44),m4)
       ksym   = iand(ishft(iwork,n40),m4)
       l1s    = iand(ishft(iwork,n20),m10)
       l1f    = iand(ishft(iwork,n10),m10)
       k1s    = iand(iwork,m10)
       iwork=ibuffer(3)
       k1f    = iand(ishft(iwork,n20),m10)
       k1l1   = iand(iwork,m20)
*@else
*       iwork=ibuffer(1)
*       intblk8 = (iand(ishft(iwork,n44),m20))
*       stn1218 = (iand(ishft(iwork,n24),m20))
*       stn38   = (iand(ishft(iwork,n4),m20))
*       nintkl8 = (ishft(iand(iwork,m4),n16))
*       iwork=ibuffer(2)
*       nintkl8 = (ior(iand(ishft(iwork,n48),m16),nintkl8))
*       lsym8   = (iand(ishft(iwork,n44),m4))
*       ksym8   = (iand(ishft(iwork,n40),m4))
*       l1s8    = (iand(ishft(iwork,n20),m10))
*       l1f8    = (iand(ishft(iwork,n10),m10))
*       k1s8    = (iand(iwork,m10))
**
*       iwork=ibuffer(3)
*       k1f8    = (iand(ishft(iwork,n20),m10))
*       k1l18   = (iand(iwork,m20))
**
*      intblk=intblk8
*      stn121=stn1218
*      stn3=stn38
*      nintkl=nintkl8
*      lsym=lsym8
*      ksym=ksym8
*      l1s=l1s8
*      l1f=l1f8
*      k1s=k1s8
*      k1f=k1f8
*      k1l1=k1l18 
*@endif
*@ifdef debug
*       call printbinary(ibuffer(1),'dcode4x ibuffer(1)')
*       call printbinary(ibuffer(2),'dcode4x ibuffer(2)')
*       call printbinary(ibuffer(3),'dcode4x ibuffer(3)')
*@endif

      return
      end

      subroutine encodf(
     & nv,     nsv,    symw,   head,
     & tail,   icode,  val,    isv,
     & yb,     yk,     nw,     wword,
     & isv2,   imode )
c**********************************************************************
c
c  pack formula info for a loop into a structure.
c
c  input:
c  nv = number of loop values.
c  nsv = number of symmetry versions of the loop.
c  symw = loop symmetry.
c  head = head within the drt.
c  tail = tail within the drt.
c  icode = loop code.
c  val(1:nv) = loop values.
c  isv(1:nsv) = contributing symmetry versions (bra).
c  isv(1:nsv) = contributing symmetry versions (ket).
c  imode = 0: isv(bra) = isv(ket) (mcscf)
c  imode = 1: isv(bra) <> isv(ket) (transft,transmom)
c  yb(1:nsv) = bra loop weights for each symmetry.
c  yk(1:nsv) = ket loop weights for each symmetry.
c
c  output:
c  nw = number of packed entries.
c  wword(1:nw) = working precision buffer of packed entries.
c
      implicit none 
c     # dummy:
      integer nv, nsv, symw, head, tail, icode, nw
      integer isv(8), yb(8), yk(8) , isv2(8), imode
C
c     # local:
      integer*8 i, j, i1, i2, i3, i4, i5, i6, i7, i8, mask, ipt
c
c     # bummer error types.
      integer   wrnerr,  nfterr,  faterr
      parameter(wrnerr=0,nfterr=1,faterr=2)
c
c     # set up any required bit operators...
c
      integer*8 val(3), wword(8), wwords,iword
      integer*8   m4,        m8,        m16
      parameter(m4=2**4-1, m8=2**8-1, m16=2**16-1)

      integer*8  ip4, ip8, ip16, ipack8, ipack6, sand
*@ifdef int64
      ip4(i,j)  = ior(ishft(i,4),iand(j,m4))
      ip8(i,j)  = ior(ishft(i,8),iand(j,m8))
      ip16(i,j) = ior(ishft(i,16),iand(j,m16))
*@else
*      ip4(i,j)  = ior(ishft(i,int(4,8)),iand(j,m4))
*      ip8(i,j)  = ior(ishft(i,int(8,8)),iand(j,m8))
*      ip16(i,j) = ior(ishft(i,int(16,8)),iand(j,m16))
*@endif
      ipack8(i1,i2,i3,i4,i5,i6,i7,i8) = ip16(
     & ip8( ip4(i1,i2), ip4(i3,i4) ),
     & ip8( ip4(i5,i6), ip4(i7,i8) )  )
      ipack6(i1,i2,i3,i4,i5,i6) = ip16(
     & ip8( ip4(i1,i2),i3), ip8(i4,ip4(i5,i6)) )
      sand(i,j,mask) = iand( ishft(i,-j), mask )

c  error check

        if (imode.ne.0 .and. imode.ne.1)
     .      call bummer('encodf invalid imode, imode=',imode,2)
c

c     # pack the formula entries...
c

*@ifdef int64
      nw = (nv+1) + (nsv+1)/2 + imode
      wword(1) = ior(ishft(
     & ipack8(isv(1),isv(2),isv(3),isv(4),
     & isv(5),isv(6),isv(7),isv(8)),32),
     & ipack6(nsv,symw,head,tail,nw,icode))
       if (head.gt.255 .or. tail.gt.255) then
           call bummer('head/tail failure in encodf',255,2)
       endif 
      if (imode.eq.1) then
      wword(2) = ipack8(isv2(1),isv2(2),isv2(3),isv2(4),
     .                  isv2(5),isv2(6),isv2(7),isv2(8))
      endif
*@else
*      nw = (nv+1) + (nsv+1)/2 + imode
*      wword(1) = ior(ishft(
*     & ipack8(int(isv(1),8),int(isv(2),8),
*     .        int(isv(3),8),int(isv(4),8),
*     &        int(isv(5),8),int(isv(6),8),
*     .        int(isv(7),8),int(isv(8),8)),int(32,8)),
*     & ipack6(int(nsv,8),int(symw,8),
*     .        int(head,8),int(tail,8),
*     .        int(nw,8),int(icode,8)))
*      if (imode.eq.1) then
*      wword(2) = ipack8(int(isv2(1),8),int(isv2(2),8),
*     .                  int(isv2(3),8),int(isv2(4),8),
*     .                  int(isv2(5),8),int(isv2(6),8),
*     .                  int(isv2(7),8),int(isv2(8),8))
*      endif
*@endif
      do 10 i = 1, nv
         wword(i+1+imode) = val(i)
10    continue
CC
      ipt = 1
      do 20 i = 1, ((nsv+1)/2)
*@ifdef int64
         if (yb(isv(ipt)).gt.65535 .or. yb(isv(ipt+1)).gt.65535 
     .   .or. yk(isv2(ipt)).gt.65535 .or. yk(isv2(ipt+1)).gt.65535) then
          call bummer('yb/yk gt. 65535 failure',65535,2)
         endif 
         wword(1+nv+i+imode) = ior(ishft(
     &    ip16(yb(isv(ipt  )),yk(isv2(ipt  ))), 32),
     &    ip16(yb(isv(ipt+1)),yk(isv2(ipt+1)))      )
       
*@else
*         wword(1+nv+i+imode) = ior(ishft(
*     &    ip16(int(yb(isv(ipt  )),8),int(yk(isv2(ipt  )),8)), 
*     .    int(32,8)),
*     &    ip16(int(yb(isv(ipt+1)),8),int(yk(isv2(ipt+1)),8))      )
*@endif
         ipt = ipt + 2
20    continue
c
       return
c**********************************************************************
c deck pcodef
      entry pcodef( icode, nw, wwords )
c**********************************************************************
c
c     # pack a code value into the formula structure.
c     # this is used for icode=0 and icode=1.
c
      nw     = 1
*@ifdef int64
      wwords = ip4(nw,icode)
*@else
*      wwords = ip4(int(nw,8),int(icode,8))
*@endif
      return
c
c**********************************************************************
c deck decodf
      entry decodf(
     & nv,     nsv,    symw,   head,
     & tail,   icode,  val,    isv,
     & yb,     yk,     nw,     wword ,
     & isv2, imode )
c**********************************************************************
c
c     # unpack a packed formula structure.
c     #
c     # input:
c     # wword(*) = packed structure.
c     #
c     # output:
c     # icode = loop code.
c     # nw = number of packed words.
c     # nv,nsv,symw,head,tail,val(*),isv(*),yb(*),yk(*) = formula
c     #                           entries. (computed only if icode>1 )
c
      iword = wword(1)
      icode = int(sand(iword,int(0,8),m4),4)
      nw    = int(sand(iword,int(4,8),m4),4)
      if ( icode .gt. 1 ) then
         tail   = int(sand(iword,int(8,8),m8),4)
         head   = int(sand(iword,int(16,8),m8),4)
         symw   = int(sand(iword,int(24,8),m4),4)
         nsv    = int(sand(iword,int(28,8),m4),4)
         isv(1) = int(sand(iword,int(60,8),m4),4)
         isv(2) = int(sand(iword,int(56,8),m4),4)
         isv(3) = int(sand(iword,int(52,8),m4),4)
         isv(4) = int(sand(iword,int(48,8),m4),4)
         isv(5) = int(sand(iword,int(44,8),m4),4)
         isv(6) = int(sand(iword,int(40,8),m4),4)
         isv(7) = int(sand(iword,int(36,8),m4),4)
         isv(8) = int(sand(iword,int(32,8),m4),4)
CC
         if (imode.eq.1) then
         iword = wword(2)
         isv2(1) = int(sand(iword,int(28,8),m4),4)
         isv2(2) = int(sand(iword,int(24,8),m4),4)
         isv2(3) = int(sand(iword,int(20,8),m4),4)
         isv2(4) = int(sand(iword,int(16,8),m4),4)
         isv2(5) = int(sand(iword,int(12,8),m4),4)
         isv2(6) = int(sand(iword,int( 8,8),m4),4)
         isv2(7) = int(sand(iword,int( 4,8),m4),4)
         isv2(8) = int(sand(iword,int( 2,8),m4),4)
         endif

         nv = nw - 1 - (nsv + 1) / 2 -imode
         do 210 i = 1, nv
            val(i) = wword(i+1+imode)
210      continue
         ipt = 1
         do 220 i = 1, (nsv+1)/2
            iword     = wword(1+nv+i+imode)
            yb(ipt)   = int(sand( iword, int(48,8), m16 ),4)
            yk(ipt)   = int(sand( iword, int(32,8), m16 ),4)
            yb(ipt+1) = int(sand( iword, int(16,8), m16 ),4)
            yk(ipt+1) = int(sand( iword, int( 0,8), m16 ),4)
            ipt       = ipt + 2
220      continue
      endif
      return
      end

      function getstepventry(nbfn,occ,bfn,ref)

      integer getstepventry,nbfn,occ(*),bfn,ref
      integer ibit,iword,itmp
      integer wordlen
*@ifdef int64
        parameter(wordlen=64)
*@else
*      parameter(wordlen=32)
*@endif


      ibit=(ref-1)*nbfn*2 + (bfn)*2
      iword=(ibit-1)/wordlen
c     write(*,*) 'getstv: getting ibit,iword=',ibit,iword+1
      ibit = ibit-iword*wordlen
      iword=occ(iword+1)

c       extract the two bits no ibit and ibit-1 from iword


      itmp=iand(3,ishft(iword,-ibit+2))
c   write(*,*) 'getstepventry: extracting bfn,ref',bfn,ref,'=',itmp
      getstepventry=itmp

      return
      end

ctm in2 stems from ciden9.f


      function in2(i)
c     this function gives the number of full real words for a given
c     number of integer words
*@ifdef int64
      in2 = i
*@else
*       in2=(i-1)/2+1
*@endif
      return
      end

ctm in stems from ciden9.f

      function in(i)
c     this function ensures that for ibm machines an integer array
c     ends at a double word boundary
*@ifdef int64
      in = i
*@else
*       in=((i-1)/2)*2+2
*@endif
      return
      end

ctm labget stems from ciden2.f

      subroutine labget(word,lab,ipos)
c
c  this routine gets the integer label 'lab' from the 'ipos' position
c  of the working precision word 'word'.
c
       implicit integer(a-z)
*@ifdef int64
      integer iword(2),lab,ipos
      integer m32
      parameter(m32=2**32-1)
*@if defined t3e64 || defined cray
*CC
*      iword(1)=shiftr(word,32)
*      iword(2)=and(word,m32)
*      lab=iword(ipos)
*@else 
        iword(1)=ishft(word,-32)
        iword(2)=iand(word,m32)
        lab=iword(ipos)
*@endif
*@else
*      integer word(2)
*      lab=word(ipos)
*@endif
      return
      end

ctm labgt1 stems from ciden2.f

      subroutine labgt1(word,lab,ipos)
c
c  this routine gets the integer label 'lab' from the 'ipos' position
c  of the working precision word 'word'.
c
      implicit integer(a-z)
*@ifdef int64
      integer iw16, word
      parameter(m16=2**16-1)
*@if defined t3e64 || defined cray
*      iw16=and(shiftr(word,16*(4-ipos)),m16)
*      lab=iw16
*@else                                          
        iw16=iand(ishft(word,16*(ipos-4)),m16)
        lab=iw16
*@endif
*@else
*       real*8 word
*       real*8 dpword
*      integer*2 iw16(4)
*      equivalence (iw16(1),dpword)
**
*      dpword=word
*      lab=iw16(ipos)
*@endif
      return
      end

ctm labpt1 stems from ciden9.f

      subroutine labpt1(word,lab,ipos)
c
c  this routine puts the integer label 'lab' into the 'ipos' position
c  of the working precision word 'word'.
c
      implicit integer(a-z)
*@ifdef int64
      integer iw16(4),word
      integer m16
      parameter(m16=2**16-1)
*@if defined t3e64 || defined cray
*      integer dpack4
*      dpack4(i1,i2,i3,i4)=
*     +  or(shiftl(or(shiftl(i1,16),i2),32),or(shiftl(i3,16),i4))
*CC
*      iw16(1)=shiftr(word,48)
*      iw16(2)=and(shiftr(word,32),m16)
*      iw16(3)=and(shiftr(word,16),m16)
*      iw16(4)=and(word,m16)
*      iw16(ipos)=and(lab,m16)
*      word=dpack4(iw16(1),iw16(2),iw16(3),iw16(4))
*@else                                             
      integer dpack4
      dpack4(i1,i2,i3,i4)=
     +  ior(ishft(ior(ishft(i1,16),i2),32),ior(ishft(i3,16),i4))
CC
      iw16(1)=ishft(word,-48)
      iw16(2)=iand(ishft(word,-32),m16)
      iw16(3)=iand(ishft(word,-16),m16)
      iw16(4)=iand(word,m16)
      iw16(ipos)=iand(lab,m16)
      word=dpack4(iw16(1),iw16(2),iw16(3),iw16(4))
*@endif
*@else
*      real*8 word
*      integer*2 iw16(4)
*      real*8 dpword
*      equivalence (iw16(1),dpword)
*      dpword=word
*      iw16(ipos)=lab
*      word=dpword
*@endif
      return
      end

ctm labput stems from ciden9.f

      subroutine labput(word,lab,ipos)
c
c  this routine puts the integer label 'lab' into the 'ipos' position
c  of the working precision word 'word'.
c
      implicit integer(a-z)
*@ifdef int64
       integer word
       integer iword(2)
       integer m32
      parameter(m32=2**32-1)
CC
*@if defined t3e64 || defined cray
*       iword(1)=shiftr(word,32)
*       iword(2)=and(word,m32)
*       iword(ipos)=and(lab,m32)
*       word= or(shiftl(iword(1),32),iword(2))
*@else                                           
       iword(1)=ishft(word,-32)
       iword(2)=iand(word,m32)
       iword(ipos)=iand(lab,m32)
       word= ior(ishft(iword(1),32),iword(2))
*@endif
*@else
*       integer word(2)
*       integer lab,ipos
*       word(ipos)=lab
*@endif
      return
      end

ctm low32 stems from ciden9.f

      integer function low32( wkword )
c
c     returns the "lowest" 32 bits of the working precision word wkword
c
*@ifdef int64
*@if defined t3e64 || defined cray
*      integer wkword
*      low32 = and( wkword, mask(96) )
*@else
        integer wkword,m32
        parameter (m32=2**32-1)
        low32=iand(wkword,m32)
*@endif
*@else
*      integer*4 wkword(2)
*      low32=wkword(2)
*@endif
      return
      end

c changed 13.1.98
      subroutine packlp( word, code, head, bver, ybt, ykt)
ctm
ctm old version
ctm   subroutine packlp( word, code, head, xbt, ybt, ykt)
ctm
c
c  28-jun-90 wlkxmx->xbarmx, wlkmx->nwlkmx changes. -rls
c  09-apr-90 rwlphd->head change completed. -rls
c  17-apr-89 ft packing check moved here from wloop. -rls
c  modification of packing in  2-12-87 (p. szalay)
c                       (also in ci program!!!)
c
c  note:  if bits become scarce, then xbt may be eliminated in the
c         ft since it depends on head which is now explicitly
c         specified in each ft entry.
c         of course, the programs which read the ft will require
c         modification consistent with this change also. -rls
c
      implicit integer(a-z)
c
      integer   nrowmx,        xbarmx,         nwlkmx
      parameter(nrowmx=1023    , nwlkmx=2**20-1)

c      code  head  bver  ykt  ybt
c        4    10    2    24   24             64bit
c        4    10    2   16/ 8  24            2*32bit

      integer code, head, ybt, ykt
      integer bver

c     # for enforcement of ft packing limits. note that xbarmx
c     # may be 2 to 4 times smaller than nwlkmx since xbt is only
c     # addressed above the first level and only for a single chaining
c     # scheme at any one time.  Specifically xbt=xbar(head,ichain)
c     # whereas nwlkmx = sum( xbar(1:4,1:3) ). -rls
c
c
c     # bummer error types.
      integer   wrnerr,  nfterr,  faterr
      parameter(wrnerr=0,nfterr=1,faterr=2)
c
*@ifdef int64
      integer word,word1
       integer oder,links
       oder(i,j) = ior(i,j)
       links(i,j) = ishft(i,j)
*@else
*      real*8 word, word1
*c     # for 32-bit integer machines...
*      integer iword(2)
*      equivalence (word1,iword(1))
*@endif

c
c
c     # check items to be packed...
*@ifdef debug
*      if(head.gt.nrowmx)then
*         call bummer('packlp: head=',head,faterr)
*      elseif(ybt.gt.nwlkmx)then
*         call bummer('packlp: ybt=',ybt,faterr)
*      elseif(ykt.gt.nwlkmx)then
*         call bummer('packlp: ykt=',ykt,faterr)
*      endif
*@endif
c
*@ifdef int64
      word1=oder(oder(oder(oder(links(ybt,24),ykt),links(bver,48)),
     + links(head,50)),links(code,60))
*@else
*@if defined  (oldoldsun) || defined ( fps)
*      iword(1)=or(or(or(lshift(code,28),lshift(head,18)),
*     + lshift(bver,16)),rshift(ybt, 8))
*      iword(2)=or(lshift(ybt,24),ykt)
*@else
*c  32-bit integer machines with f90 bit operators...
*      iword(1)=ior(ior(ior(ishft(code,28),ishft(head,18)),
*     + ishft(bver,16)),ishft(ybt,-8))
*      iword(2)=ior(ishft(ybt,24),ykt)
*@endif
*@endif
c
      word = word1
c
      return
      end


      subroutine packstep( nintorb, istepv, output,stepno )
c
c ==============================================================
c compress each step vector istepv into nintorb*2 bits
c starting from (stepno-1)*nintorb*2+1 ... stepno*nintorb*2
c
c  input: nintorb   number of internal orbitals
c         istepv    step vector
c         stepno    number of the valid z walk
c  output: output   packed step vector
c
c          1 ... nintorb (ivalz=1)
c          1 ... nintorb (ivalz=2) ....
c
c tm october, 1st 1999
c===============================================================
      implicit logical(a-z)
       integer stepno, nintorb,istepv(nintorb),output(*)

      integer firstbit,firstwordm1,ii,intorb,ibit
      integer wordlen
*@ifdef int64
        parameter(wordlen=64)
*@else
*      parameter(wordlen=32)
*@endif


      firstbit=(stepno-1)*nintorb*2
      firstwordm1= (firstbit)/wordlen

      ibit=firstbit-firstwordm1*wordlen

      intorb=1
      ii=firstwordm1
      ii=ii+1

c  output must have been initialized with zeros!
100   output(ii) = ior(output(ii),ishft(istepv(intorb),ibit))
      ibit=ibit+2
      intorb=intorb+1
      if (intorb.gt.nintorb) goto 999
      if (ibit.gt.(wordlen-2)) then
         ibit=0
         ii=ii+1
      endif
      goto 100
999   return
      end
c
ctm This file contains a collection of various strongly machine
ctm dependent packing routines - on a long-term basis they
ctm have to generalized and to be included into colib3.f
ctm
c deck rl
      integer function rl(i)
c
c     this function maps real word lengths on integer word lengths
c
c**********************************************************************
c  the use of this function is dangerous, due to machine-dependent data
c  alignment constraints, and should be eliminated in the calling
c  programs. -rls
c**********************************************************************
c
      integer i
c
c      integer first
c      save    first
c      data first / 1 /
c
c      if ( first .eq. 1 ) then
cc        # print an annoying message the first time through. -rls
c         call bummer(' rl(): obsolete function called.'
c     &    //' the calling program should be updated', 0, 0)
c         first = 0
c      endif
c
*@ifdef int64
      rl= i
*@else
*c     # assuming i*4 declarations
*      rl=2*i
*@endif
      return
      end

*@ifdef ciformulatape
*      subroutine unplp(word, code, rwlphd, wtlphd, mlp2, mlp1, xbar)
*c
*c     modification of packing   2-12-87 (p. szalay)
*c                       (also in ftape program!!!)
*c
*      implicit logical(a-z)
*      integer nrowmx
*      parameter (nrowmx=1023)
*c
*      integer code, rwlphd, wtlphd, mlp2, mlp1
*      integer bver, xbar(nrowmx,3)
*c
*      integer ybta, ybtb,mask10,mask24
*      real*8 word1
*      integer mask4, mask8, mask14, mask7,mask16, mask18, maskex
*      parameter (mask4=2**4-1, mask8=2**8-1,mask10=2**10-1)
*      parameter (mask14=2**14-1, mask16=2**16-1, mask7=2**7-1)
*      parameter (mask24=2**24-1, maskex=2**24-2**8 )
*@ifdef int64
*       integer word
*        integer m4,m8,m16
*        parameter (m4=2**4-1, m8=2**8-1, m16=2**16-1)
*       code   = iand(ishft (word,-60),mask4)
*       rwlphd = iand(ishft (word,-50),mask10 )
*       bver   = iand(ishft (word,-48),3)
*       wtlphd = xbar(rwlphd,bver)
*       mlp2   = iand(ishft (word,-24),mask24  )
*       mlp1   = iand(word,mask24   )
*@else
*c     # 32-bit machines with vax/ibm/mil-std-1753/f90 bit operators.
*      integer   iand, ishft, ior
*      intrinsic iand, ishft, ior
*      real*8 word
*      integer iword(2)
*      equivalence (word1,iword(1))
*      word1 = word
*      code   = iand(ishft(iword(1),-28),mask4)
*      rwlphd = iand(ishft(iword(1),-18),mask10)
*      bver   = iand(ishft(iword(1),-16),3)
*      wtlphd= xbar(rwlphd,bver)
*      ybta   = iand(ishft(iword(1), 8),maskex)
*      ybtb   = iand(ishft(iword(2),-24),mask8 )
*      mlp2   = ior(ybta,ybtb)
*      mlp1   = iand(iword(2),mask24)
*@endif
*      return
*      end
*@endif
c

      subroutine uncmprlimvec(limvec,nwalk,valid,npos,maxval)
      implicit none
c     input: compressed limvec (length encoding)
c            nwalk   (length of uncompressed limvec)
c            npos    (length of compressed limvec)
c            valid   LAST compressed item corresponds to
c                    valid (1) or invalid (0) walk
      integer nwalk,limvec(nwalk),cpos,curpos2,maxval,lvalidold
      integer lvalid,j,i,ii,curval,curpos,valid,npos,lvalidnew
c     in order to operate in place we need to start from the end
c             the last entry: mod(npos,2)=1   valid
c                             mod(npos,2)=0  .not. valid
c
c there are two special cases:
c   (1) the entire limvec has to be decompressed
c   (2) the z-part has to be decompressed
c
c     valid= 0 or 1
      do i=npos+1,nwalk
          limvec(i)=0
      enddo
c     write(6,*) 'uncmprlimvec: nwalk,valid,npos,maxval=',
c    .  nwalk,valid,npos,maxval 
c     write(6,*) 'uncmprlimvec: initial limvec='
c     write(6,'(20i4)') (limvec(i),i=1,npos)
      curpos=nwalk
c     value of last valid walk
      lvalidnew=valid
      lvalidold=valid
      do i=npos,1,-1
          if (limvec(i).eq.maxval) then
             curpos2=curpos-maxval+1
c            take previous lvalid 
             lvalid=lvalidold
          else
             curpos2=curpos-limvec(i)
             lvalid =lvalidnew
c            switch lvalid 
          endif
         lvalidold=lvalid    
         lvalidnew=mod(lvalid+1,2)
c        write(0,111) i,curpos2+1,curpos,lvalid
c111     format('cpos=',i6,' overwriting ',2i7,' with ',i1)
         do j=curpos2+1,curpos
           limvec(j)=lvalid
         enddo
          curpos=curpos2
         if (curpos.lt.0) call bummer('uncmprlimvec failed',curpos,2)
c     write(6,*) 'uncmprlimvec: step ',i,'curpos=',curpos,'limvec='
c     write(6,'(80i1)') (limvec(ii),ii=1,nwalk)
      enddo
c      write(0,*) 'curpos,curpos2=',curpos,curpos2 
c     write(6,*) 'uncmprlimvec: limvec='
c     write(6,'(80i1)') (limvec(i),i=1,nwalk)

      return
      end

      subroutine cmprlimvec(limvec,nwalk,cpos,ifac,maxval)
c     #  use length encoding for limvec right now
c     #  negative number: n consecutive invalid walks
c     #  positive number: n consecutive valid walks
c     #  first number indicates sign of the first limvec entry
c     #  0: invalid   1:valid
c     #  subsequently we can avoid the sign
c     #  add a 0 termination
c     #  add maxval:  a value of maxval in the compressed
c                     array indicates maxval-1 values of the
c                     same sign to be continued   
c
c     input 10111000  maxval=10
c     output       1 1 3 3   cpos= 4 ifac=1 : contracted vector of length 4 
c                                             first value corresponds to 1's
c     input 10111000   maxval=3 
c     output  1 1 3 1 3 1    cpos = 6 ifac=1  entries with val=maxval 
c                                             indicate maxval-1 valid of current
c                                             sign to be continued 
c     CHANGE: ifac is the value of the LAST index vector entry (0 or 1)
c
c      let's try to find a good compression automatically
c
      implicit none
      integer nwalk,limvec(nwalk),cpos,maxval
      integer imaxval,iwalk,ll,i,l,curval,curpos,ifac
      integer  mmval(4),icheck
      real*8 cmprfactor,savecmprfactor
      data mmval /9,99,999,9999/

      ifac=limvec(nwalk)
c
c     find good compression factor 
c
 111  format(50i2)
      savecmprfactor=0.0d0 
      do l=1,4
      cpos=0
      curval=limvec(1)
      curpos=1
      icheck=0
c     write(0,*) 'initial limvec, l=',l
c     write(0,111) limvec(1:nwalk)
      do i=2,nwalk
        if (limvec(i).ne.curval) then
          cpos=cpos+1
c         limvec(cpos)=(i-curpos)
          icheck=icheck+i-curpos
          curpos=i
          curval=limvec(i)
c         write(0,112) 'cpos,curpos,i,icheck=',
c    .     cpos,curpos,i,icheck
 112  format(a,4i6)
        elseif (i-curpos.eq.mmval(l)) then
          cpos=cpos+1
c         limvec(cpos)=mmval(l)
          curpos=i-1
          icheck=icheck+mmval(l)-1
c         write(0,112) 'MM cpos,curpos,i,icheck=',
c    .     cpos,curpos,i,icheck
        endif
      enddo
      cpos=cpos+1
c     limvec(cpos)=(nwalk-curpos+1)
       icheck=icheck+nwalk-curpos+1
      if (cpos.gt.nwalk) call bummer('cmprlimvec failed',cpos,2)
      if (icheck.ne.nwalk) 
     .     call bummer('cmprlimvec failed icheck=',icheck,2)
c     compression factor (conservative) in bytes
c     100 - cpos*l*100/nwalk  since 01 index requires 1 byte per walk   
      cmprfactor=100d0- dble(cpos*(l))*100.d0/dble(nwalk)
      write(6,100) nwalk,cpos,mmval(l),cmprfactor
 100  format('nwalk=',i8,' cpos=',i8,' maxval=',i5,
     .    ' cmprfactor=',f8.2,' %.')
      if (cmprfactor-savecmprfactor.lt.1.0d0) then
        maxval=mmval(max(1,l-1))
        ll=max(1,l-1)
        goto 500
      else
        maxval= mmval(l)
        ll=l
      endif
       savecmprfactor=cmprfactor
      enddo
 500  continue 

c
c    compress index vector
c

      imaxval=0
      cpos=0
      curval=limvec(1)
      curpos=1
      icheck=0
      do i=2,nwalk
c         write(0,*) 'testing i,limvec(i)=',i,limvec(i)
        if (limvec(i).ne.curval) then
          cpos=cpos+1
          icheck=icheck+i-curpos
          limvec(cpos)=(i-curpos)
          curpos=i
          curval=limvec(i)
c         write(0,112) 'VNORM cpos,i-1,icheck=',
c    .     cpos,i-1,icheck
        elseif (i-curpos+1.eq.maxval) then
          cpos=cpos+1
          imaxval=imaxval+1
          icheck=icheck+maxval-1
          limvec(cpos)=maxval   
          curpos=i
c         write(0,112) 'VMAX cpos,curpos,i,imaxval=',
c    .     cpos,curpos,i,imaxval
        endif
      enddo
      cpos=cpos+1
      limvec(cpos)=(nwalk-curpos+1)
       icheck=icheck+nwalk-curpos+1
      if (cpos.gt.nwalk) call bummer('cmprlimvec failed',cpos,2)
      cmprfactor=100d0- dble(cpos*(ll))*100.d0/dble(nwalk)
      write(6,101) nwalk,cpos,maxval,cmprfactor
 101  format(' compressed with: nwalk=',i8,' cpos=',i8,' maxval=',i5,
     .    ' cmprfactor=',f8.2,' %.')
      return
      end


     
       subroutine printbinary (w,comment)
        integer*8 w,wloc(64),m,m1
        character*(*) comment

        integer*8 icnt

         m=0
         m1=1
         do icnt=1,64
          wloc(icnt)=iand(ishft(w,m),m1)
          m=m-1
         enddo
         write(6,50) comment,(wloc(icnt),icnt=1,64)
 50      format(a,2x,64i1)
         return
         end
 
