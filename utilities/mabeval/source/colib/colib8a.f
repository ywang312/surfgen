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

      subroutine sif_baslab(nsym,nbpsy,slabel,baslab,prefix)
      implicit none
      integer nsym,nbpsy(nsym)
      character*8 temp,baslab(*)
      character*3 prefix
      character*4 slabel(8)
      integer isym,i,n,nn,k,j

c
c     constructs basis label vector
c     nsym : number of irreps
c     nbpsy: number of basis functions per irrep
c     baslab: output basis label vector
c     prefix: prefix for basis label vector
c            'SAO'  ->   SAO_NNNN
c            'MO '  ->   MO__NNNN
c            'SYM'  ->   SSS_NNNN
c

     
       n=0
       nn=0
       do isym=1,nsym
         if (prefix(1:3).eq.'SYM') nn=0
         do j=1,nbpsy(isym)
          n=n+1
          nn=nn+1
c       rechtsbuendig n
          write(temp,'(i8)') nn
c       substitute spaces by underscores
          do k=1,8
             if (temp(k:k).eq.' ') temp(k:k)='_'
          enddo
          baslab(n)=temp
         enddo
       enddo

       if (prefix(1:2).eq.'MO') then
          do i=1,n
           baslab(i)(1:2)='MO'
          enddo
       elseif (prefix(1:3).eq.'SAO') then
          do i=1,n
          baslab(i)(1:3)='SAO'
          enddo
       elseif (prefix(1:3).eq.'SYM') then
          n=0
          do isym=1,nsym
           temp(1:4)=adjustl(slabel(isym))
           do k=1,4
             if (temp(k:k).eq.' ') temp(k:k)='_'
           enddo
           do j=1,nbpsy(isym)
            n=n+1
            baslab(n)(1:4)=temp(1:4)  
           enddo
          enddo
         endif
 
         return
         end

       subroutine sif_energy(name,value,n,energy,ietype)
       implicit none
       integer i,itypea,itypeb,ierr,ietype(*), n
       character*(*) name
       real*8 value,energy(*)

       do i=len_trim(name)+1,8
           name(i:i)=' '
       enddo 

       call siftyp(itypea,itypeb,name,ierr)
       if (ierr.ne.0) then 
            write (6,*) 'name:',name,'#'
            call bummer('invalid name',ierr,2)
       endif 

       n=n+1
       energy(n)=value
       ietype(n)=itypea*1024+itypeb
       write(6,*) name,' >', itypea,itypeb,ietype(n)
       return
       end 
 
      subroutine sif_finalize(fhdle)
      implicit none
      integer fhdle
      include 'sif_data.h'
c
c     purpose: finalize an open SIFS file
c              clean status%
c 

      close (unit=status%unitno(1,fhdle))
      if (status%info(1,fhdle).eq.2)
     .    close (unit=status%unitno(2,fhdle))
      
      status%unitno(1:2,fhdle)=-1
      status%nbft(fhdle)=-1
      status%nsym(fhdle)=-1
      status%nbpsy(1:8,fhdle)=-1
      status%info(1:7,fhdle)=-1
      return
      end

      subroutine sif_init(rwflag,fname,fhdle,nsym,nbpsy,
     .                recordlen,map,nmap,imtype,
     .                energy,nenergy,ietype,
     .                title,ntitle,imode,
     .                slabel,
     .                baslab)
c
c       the last four lines of arguments are optional
c

c
c     Inititalize all parameters for a SIFS-File
c      * compute info header parameter
c      * compute title info 
c      * compute symmetry labels, basis function labels and map vectors
c      * open file and write header records
c      * maintain status of the SIFS-File
c      STATUS=nsym,nbpsy,info 
c
       implicit none
       character*(*) fname
       character*80 fname2
       character*4 rwflag
c      'new1','new2','read'
       integer nmap,nbft,ifmt1,lrec1,nmax1,ierr,ifmt2,lrec2,nmax2
       integer nenrgy
       integer recordlen,fhdle,nsym,nbpsy(nsym),ntitle,map(nmap,*)
       character*80  title(ntitle)
       integer nenergy,ietype(nenergy),imtype(nmap)
       integer i,aoints,aoints2,imode
       character*4 :: slabel(8)
       character*4 slabelloc(8)
       character*8 ::  baslab(*)
       character*8 baslabloc(511)
       logical lexist,lopen
       real*8  energy(nenergy)
       integer locntitle,locnsym,locnbft,locninfo,locnenergy,locnmap
       
       include 'sif_data.h'
       logical first
       data first /.true./
       save first

       if (first) then
          status%info=-1
          status%unitno=-1
          status%nbft=-1
          status%nsym=-1
          first=.false.
       endif
c
c      get free handle
c

       do i=1,10
        if (status%unitno(1,i).eq.-1) then
            fhdle=i
            exit
        endif
       enddo

c
c      get unit no
c
       do i=10,99
         inquire(unit=i, opened=lopen) 
         if (.not.lopen) then
              if ( status%unitno(1,fhdle).lt.0)  then
                   status%unitno(1,fhdle)=i
              elseif (status%unitno(2,fhdle).lt.0) then
                   status%unitno(2,fhdle)=i
                   exit
              endif
         endif
       enddo

       write(6,*) 'assigned unitno=',status%unitno(1:2,fhdle)



       if (rwflag(1:3).eq.'new') then 

c
c      total number of basis functions
c
        nbft=0
        do i=1,nsym
         nbft=nbft+nbpsy(i)
         status%nbpsy(i,fhdle)=nbpsy(i)
        enddo
        status%nbft(fhdle)=nbft
        status%nsym(fhdle)=nsym

c
c      irrep labels
c

c      if (.not. present(slabel)) then
       if (iand(imode,1).eq.0) then
        do i=1,nsym
         write(slabelloc(i)(1:4),'(a3,i1)') 'SYM',i
        enddo
       else
        do i=1,nsym
        slabelloc(i)(1:4)=slabel(i)(1:4)
        enddo
       endif
c
c      basis function labels
c
  
c       if (.not. present(baslab)) then
        if (iand(imode,2).eq.0) then
          do i=1,min(9,nbft)
           write(baslabloc(i)(1:8),'(a7,i1)') 'BFN____',i
          enddo
          do i=10,min(99,nbft)
           write(baslabloc(i)(1:8),'(a6,i2)') 'BFN___',i
          enddo
          do i=100,min(999,nbft)
           write(baslabloc(i)(1:8),'(a5,i3)') 'BFN__',i
          enddo
        else
          do i=1,nbft
           baslabloc(i)(1:8)=baslab(i)(1:8)
          enddo
        endif

c
c      create info array
c
c      recordlen=max record length
c      ifmt1: 
c      lrec: actual record length
c      namx: integrals per record
       
       call sifcfg(1,recordlen, nbft,0, ifmt1,lrec1,nmax1,ierr)
       call sifcfg(2,recordlen, nbft,0, ifmt2,lrec2,nmax2,ierr)
       if (rwflag(1:4).eq.'new1') then 
                 status%info(1,fhdle)=1
                 status%unitno(2,fhdle)=status%unitno(1,fhdle)
       elseif (rwflag(1:4).eq.'new2') then
                 status%info(1,fhdle)=2
       else
          call bummer('sif_init: invalid rwflag',0,2)
       endif 
       status%info(2,fhdle)=lrec1
       status%info(3,fhdle)=nmax1
       status%info(4,fhdle)=lrec2
       status%info(5,fhdle)=nmax2
       status%info(6,fhdle)=ifmt1
       status%info(7,fhdle)=ifmt2
       write(6,'(a,7i6)') 'WRITE INFO(*)=',status%info(1:7,fhdle) 

       inquire(file=fname,exist=lexist,opened=lopen)
       if (lexist.or.lopen)
     .     write(6,*) 'file ',fname,' will be overwritten.'
       
       open(unit=status%unitno(1,fhdle),file=fname, form='unformatted',
     .      access='sequential',status='unknown')

c      inquire(file=fname(1:len_trim(fname))//'2',
c    .           exist=lexist,opened=lopen)

        aoints2=status%unitno(2,fhdle)

       call sifo2f( status%unitno(1,fhdle), aoints2, 
     .     fname(1:len_trim(fname))//'2', status%info(1,fhdle), 
     .     status%unitno(2,fhdle), ierr )
         inquire(unit=status%unitno(1,fhdle), opened=lopen) 
         write(6,*) 'sif_init:write units=',status%unitno(1,fhdle),lopen
         inquire(unit=status%unitno(2,fhdle), opened=lopen) 
         write(6,*) 'sif_init:write units=',status%unitno(2,fhdle),lopen

       call sifwh(
     & status%unitno(1,fhdle),  ntitle,  nsym,    nbft,
     & 7,   nenergy,  nmap,    title,
     & nbpsy,   slabelloc,  status%info(1,fhdle),baslabloc,
     & ietype,  energy,  imtype,  map,
     & ierr )

       if (ierr.ne.0) call bummer('sifwh failed',ierr,2)

           

      return
      endif

      if (rwflag(1:4).eq.'read') then 

       inquire(file=fname,exist=lexist,opened=lopen)
       if (.not.lexist)
     .   call bummer(' sif_init: integral file not available',0,2)
       
       open(unit=status%unitno(1,fhdle),file=fname(1:len_trim(fname)), 
     .     form='unformatted', access='sequential',status='old')

       call sifrh1(
     & status%unitno(1,fhdle),  locntitle,  locnsym,    locnbft,
     & locninfo,   locnenergy,  locnmap,    ierr )
        
        if (ierr.ne.0) 
     .   call bummer('sif_init: sifrh1 failed',ierr,2)
        if (locnsym.ne.nsym .and. nsym.ne.-1) 
     .    call bummer('sif_init: inconsistent nsym=',locnsym,2)
        if (locnbft.ne.nbft .and. nsym.ne.-1) 
     .    call bummer('sif_init: inconsistent nbft=',locnbft,2)
       if (locninfo.gt.10) 
     .    call bummer('sif_init: inconsistent ninfo=',locninfo,2)
       if (locntitle.gt.ntitle) 
     .    call bummer('sif_init: ntitle too large, ntitle=',locntitle,2)
       if (locnenergy.gt.nenergy) 
     .    call bummer('sif_init: nenergy too large, nenergy=',
     .                 locnenergy,2)
       if (locnmap.gt.nmap) 
     .    call bummer('sif_init: nmap too large, nmap=',locnmap,2)
      


       nsym=locnsym
       status%nsym(fhdle)=nsym
       status%nbft(fhdle)=locnbft
       ntitle=locntitle
       nenergy=locnenergy
       nmap =locnmap

       call sifrh2(
     & status%unitno(1,fhdle),ntitle,  locnsym, locnbft,
     & locninfo,   locnenergy,  locnmap,    title,
     & nbpsy,   slabel,  status%info(1,fhdle),baslab,
     & ietype,  energy,  imtype,  map,
     & ierr )
       
       status%nbpsy(1:nsym,fhdle)=nbpsy(1:nsym)
       write(6,'(a,7i6)') 'READ INFO(*)=',status%info(1:7,fhdle)

       recordlen=status%info(2,fhdle)
       call sifcfg(1,recordlen, locnbft,0, ifmt1,lrec1,nmax1,ierr)
       call sifcfg(2,recordlen, locnbft,0, ifmt2,lrec2,nmax2,ierr)
       status%info(6,fhdle)=ifmt1
       status%info(7,fhdle)=ifmt2


       aoints2=status%unitno(2,fhdle)
       call sifo2f( status%unitno(1,fhdle), aoints2,
     .     fname(1:len_trim(fname))//'2', status%info(1,fhdle),
     .     status%unitno(2,fhdle), ierr )
         inquire(unit=status%unitno(1,fhdle), opened=lopen)
         write(6,*) 'sif_init:read units=',status%unitno(1,fhdle),lopen
         inquire(unit=status%unitno(2,fhdle), opened=lopen)
         write(6,*) 'sif_init:read units=',status%unitno(2,fhdle),lopen



      return
      endif

      return
      end 

