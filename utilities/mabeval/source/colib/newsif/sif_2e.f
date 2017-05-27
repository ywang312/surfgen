      subroutine sif_2e(rwflag,name,labels,values,maxint,numint,fhdle,
     .                  work,lwork,lastrec,lstflg)

      implicit none
      include 'sif_data.h'
      integer maxint,fhdle
      integer numint,lwork,labels(4,maxint)
      character*(*) name
      character*4  rwflag
      real*8 values(maxint),work(lwork)
      logical lastrec
      integer lstflg,numsave
      integer ipt(10),aoint2
      integer i,itypea,itypeb,ierr,atebyt,locnum,loclast,locitypea
      integer locitypeb,locifmt,locibvtyp,nmax2,nrec,num,last,reqnum

c
c     read/write  a chunk of a 2e array
c
c     rwflag: 'init','read','writ'                                  
c     name:  integral name
c     labels: label field
c     values: integral field
c     maxint: maximum number of integrals to be read
c             or maximum number of integrals to be written 
c     numint: number of integrals, read 
c             number of integrals remaining in arrays (1:numint)
c     work  : work space
c     lwork : length of work space
c     lastrec: last bucket of integrals to be read or to be written
c    
c
c     purpose: read/write a one-electron matrix on the sifs file
c
      do i=len_trim(name)+1,8
         name(i:i)=' '
      enddo
      call siftyp( itypea, itypeb, name, ierr )
      if (ierr.ne.0) then
         write(6,*) '#',name,'#'
         call bummer('sif_2e: siftyp failed',ierr,2)
      endif
      if (itypea.ne.3) 
     .   call bummer('sif_2e: invalid itypea',itypea,2)

c
c      ipt(1):  buffer 
c

      ipt(1)=1
      ipt(2)=ipt(1)+atebyt(status%info(4,fhdle))

      if (rwflag.eq.'init') then
          aoint2=status%unitno(2,fhdle)
c         rewind(aoint2)
          call sifsk1( aoint2, status%info(1,fhdle), ierr )
c         write(6,*) 'initialized fhdle=',fhdle,'aoint2=',aoint2
          return
      endif
      if (rwflag.eq.'read') then
       numint=0
       aoint2=status%unitno(2,fhdle)
       if (maxint.lt.status%info(5,fhdle)) 
     .  call bummer('maxint.lt.nmax2, maxint2=',maxint,2)
 100   continue
       if (status%info(5,fhdle).gt.maxint-numint) goto 150 
c      write(0,*) 'calling sifrd2'
       call sifrd2(
     & aoint2,  status%info(1,fhdle),  4,  0,
     & work(ipt(1)),  locnum, loclast,  locitypea,
     & locitypeb,  locifmt,locibvtyp,  values(numint+1),
     & labels(1,numint+1),  -1,   ierr )
       if (ierr.ne.0) call bummer('sif_2e: sifrd2 failed',ierr,2)
       if (locifmt.ne.status%info(7,fhdle))
     .     call bummer('sif_2e: locifmt.ne.ifmt2',locifmt,2)
       if (loclast.eq.2) then 
           numint=numint+locnum
           goto 150
       endif
       if (locitypeb.ne.itypeb .or. locitypea.ne.itypea) goto 100
       numint=numint+locnum
       if (loclast.eq.0) goto 100
 150   continue
       lastrec=loclast.ne.0
c      write(6,85) numint,name,lastrec
 85    format(' sif_2e: returning ',i8,2x,a8,' integrals ...',
     .    ' last record=',l8)
       return
       endif
       if (rwflag.eq.'writ') then 
       numint=0
       aoint2=status%unitno(2,fhdle)
       nmax2=status%info(5,fhdle)
       nrec=0
 200   continue
        num=min(status%info(5,fhdle),maxint-numint)
c
c    lastflag should be 1 or 2
c
        last=0
c       write(6,'(a,i6,i6,i6,i6)') 'LAB:',numint,num,maxint,last
        if (lastrec .and. numint+num.eq.maxint) last=lstflg
        if (num.lt.status%info(5,fhdle)) then
         do i=1,num 
           values(i)=values(maxint-num+i)
           labels(1,i)=labels(1,maxint-num+i)
           labels(2,i)=labels(2,maxint-num+i)
           labels(3,i)=labels(3,maxint-num+i)
           labels(4,i)=labels(4,maxint-num+i)
         enddo
         numint=0
c       write(6,*) ' before 250, last,num=',last,num,lstflg
        if (last.eq.0) then 
           numint=num
           goto 250
        endif
        endif 
       numsave=num
c       write(6,'(a,i4,a,i6)') 
c    .     'setting lastflag=',last,'offset=',numint+1
c      write(0,*) 'calling sifew2'
        call   sifew2(
     & aoint2,  status%info(1,fhdle),    4,    num,
     & last,    itypea,  itypeb,  status%info(7,fhdle),
     & 0,  values(numint+1),  labels(1,numint+1),  -1,
     & work(ipt(1)),  0,   nrec,    reqnum,
     & ierr )
       if (ierr.ne.0) call bummer('sif_2e: sifew2 failed',ierr,2)
       numint=numint+numsave
       if (.not. lastrec .or. last.eq.0 ) goto 200
 250   continue
c      write(6,86) maxint,name,numint  
 86    format(' sif_2e: writing   ',i8,2x,a8,' integrals ...',
     .    i8,'remaining in buffer=')
      endif 
      return
      end

