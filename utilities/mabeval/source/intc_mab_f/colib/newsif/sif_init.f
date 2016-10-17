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

