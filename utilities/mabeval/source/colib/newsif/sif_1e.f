      subroutine sif_1e(rwflag,name,array,mxarray,fhdle,lstflg,
     .    symoff,fcore,work,lwork)
c
c   INPUT:
c     name:    integral type to be written
c     array:   integral values
c     mxarray:  space in array
c     fhdle:   handle as returned by sif_init
c     lstflg:  last flag for the last record of this type
c     symoff:  offset for dense symmetry blocks in array
c              referenced for itypea>0 , only
c     fcore :  core contribution 
c     work  :  work space
c     lwork :  length of work space
c     rwflag:  0   read    1 writ  
c
c
      implicit none
       real*8 small
       parameter (small=1.d-12)
       character*(*) name
       character*4 rwflag
c      'read','writ'
       integer lwork,fhdle,mxarray
       real*8 array(mxarray),work(lwork),symoff(36),fcore
       integer itypea,itypeb
       integer lstflg,numtot,ierr,nrec
       integer kntin(36),kntout(36)
       integer ij,jsym,isym,mapout(255)
       integer ipt(10),aoints,i,forbyt,atebyt
       integer lasta,lastb,last

       include 'sif_data.h'



c
c     purpose: read/write a one-electron matrix on the sifs file
c
      do i=len_trim(name)+1,8
         name(i:i)=' '
      enddo 
      call siftyp( itypea, itypeb, name, ierr )
      if (ierr.ne.0) then 
         write(6,*) '#',name,'#'
         call bummer('sif_1e: siftyp failed',ierr,2)
      endif
c      status%info(1,fhdle)=2
c      status%info(2,fhdle)=lrec1
c      status%info(3,fhdle)=nmax1
c      status%info(4,fhdle)=lrec2
c      status%info(5,fhdle)=nmax2
c      status%info(6,fhdle)=ifmt1
c      status%info(7,fhdle)=ifmt2

       aoints=status%unitno(1,fhdle)
c
c     mapout() default: no mapping
c
      do i=1,status%nbft(fhdle)
         mapout(i)=i
      enddo

      if (rwflag.eq.'writ') then 
c
c     kntin : number of integrals per irrep-block
c             only referenced for itypea>0 
c     kntin(nsym*(nsym+1)/2) indicates whether there are non-vanashing 
c                            integrals in a irrep block
c     symoff(nsym*(nsym+1)/2) indicates the start of the nbpsy(isym)*nbpsy(jsym)
c                             isym/jsym irrep block in the array
c
c     symoff(i*(i-1)/2+j) = offset of irrep block (i,j)
c      if (symoff(i*(i-1)/2+j).lt.0)   no contributions 
c
       if (itypea.gt.0) then 
       ij=0
        do isym=1,status%nsym(fhdle)
          do jsym=1,isym
           ij=ij+1 
           kntin(ij)=0
           if (symoff(ij).gt.0) kntin(ij)=1
          enddo
        enddo
        else
        ij=0
        do isym=1,status%nsym(fhdle)
          ij=ij+isym
           kntin(ij)=1
        enddo
        endif
       else
        kntin(1:36)=0
       endif

       

c      ipt(1): buffer   lrec1=info(2)
c          2 : values   nmax1=info(3)
c          3 : labels   nmax1*2 
c          4 : symb     nbft     
       ipt(1)=1
       ipt(2)=ipt(1)+atebyt(status%info(2,fhdle))
       ipt(3)=ipt(2)+atebyt(status%info(3,fhdle))
       ipt(4)=ipt(3)+forbyt(status%info(3,fhdle)*2)
       ipt(5)=ipt(4)+forbyt(status%nbft(fhdle))
       if (ipt(5)-1.gt.lwork)
     .     call bummer('insufficient work mem',lwork,2)
           

       nrec=0
       if (rwflag.eq.'writ') then 
       call sifw1x(
     & aoints, status%info(1,fhdle),lstflg,itypea,
     & itypeb,           mapout,  array,
     & status%nsym(fhdle),status%nbpsy(1,fhdle),symoff,  kntin,
     & work(ipt(1)),work(ipt(2)), work(ipt(3)),  fcore,
     & small,   kntout,  numtot,  nrec,
     & ierr )
       elseif (rwflag.eq.'read') then 
       call sifr1x(
     & aoints, status%info(1,fhdle),itypea,  itypeb,
     & status%nsym(fhdle),status%nbpsy(1,fhdle),symoff, work(ipt(1)),
     & work(ipt(2)), work(ipt(3)),  mapout,   work(ipt(4)),
     & mxarray,  array,   fcore,   kntin,
     & lasta,   lastb,   last,    nrec,
     & ierr )
       else
        call bummer('sif_1e: invalid rwflag',rwflag,2)
       endif
       if (ierr.ne.0) call bummer('sif_1e failed with ierr=',ierr,2)
       if (rwflag.eq.'writ') then 
       write(6,85) numtot,name,small,nrec
 85    format(' Wrote ',i8,2x,a8,' integrals gt ',e10.4,' on ',
     .          i4,' records.')
       else
       numtot=0
       do i=1,36
          numtot=numtot+kntin(i)
       enddo
        if (numtot.gt.mxarray) 
     .   call bummer('sif_1e: numtot>mxarray ',numtot,2)

       write(6,86) numtot,name 
 86    format(' Read ',i8,2x,a8, ' integrals.')
       endif  
      return
      end

