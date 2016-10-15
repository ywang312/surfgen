       program testsif

        implicit none
        integer  numrem,fhdle2,last2,fhdle,nenergy
       character*8 name
       integer iptl,iptv,iptw,maxint,i,itypea,itypeb,ierr
       integer nsym,nbpsy(8),recordlen
       integer ietype(3),ntitle,nmap,map(2,512)
       integer numint,imtype(2)
       character*80 title(5)
       character*4 slabel(8)
       character*8 bfnlab(512)
       real*8 d1(100),fcore,energy(3),work(100000)
       integer numtot,forbyt,atebyt
       logical lastrec
       include 'sif_data.h'






       nsym=-1
       recordlen=4096
       nmap=2
       ntitle=5
       nenergy=3

       call sif_init('read','aoints',fhdle,nsym,nbpsy,recordlen,
     .               map,nmap,imtype,energy,nenergy,ietype,
     .              title,ntitle,3,slabel,bfnlab)

       write(6,85) 'NSYM',nsym
       write(6,85) 'NBPSY',nbpsy
       write(6,85) 'NMAP',nmap
       write(6,85) 'IMTYPE',imtype(1:nmap)
       write(6,85) 'NENERGY',nenergy
c      write(6,85) 'IETYPE',ietype(1:nenergy)
       do i=1,nenergy
        itypea=ietype(i)/1024
        itypeb=ietype(i)-1024*itypea
        name='--------'
        call siftyp(itypea,itypeb,name,ierr)
        write(6,83) ietype(i),itypea,itypeb,name(1:8)
 83    format(' IETYPE =',i6,' itypea=',i4,' itypeb=',i4,' -> ',a8)
       enddo
       write(6,85) 'NTITLE=',ntitle
       do i=1,ntitle
        write(6,'(a)') title(i)
       enddo
       write(6,'(a,8a4)') 'SLABEL=',slabel(1:nsym)
 
 85    format(1x,a,1x,10i6)

       call sif_init('new2','testao',fhdle2,nsym,nbpsy,recordlen,
     .                map,nmap,imtype,energy,nenergy,ietype,
     .                title,ntitle,3,slabel,bfnlab)


       name(1:8)='S1(*)   '
       fcore=0.0d0
       iptw=1
       do i=1,nsym
         iptw=iptw+atebyt(nbpsy(i)*(nbpsy(i)+1)/2) 
       enddo
       call sif_1e('read',name,work,iptw-1,fhdle,-1,-1,fcore,
     .     work(iptw),100000-iptw)
       call plblks( 'SMAT AO', work, nsym, nbpsy, 'AO',3,6)

       call sif_1e('writ',name,work,iptw-1,fhdle2,2,-1,fcore,
     .        work(iptw),100000-iptw)

       name(1:8)='1/r12   '
c
c  if n2max is different for the two files, we need some
c  extra space 
c
c  subsequently: read/write two records  in one call to sif_2e
c                need space for three records here to have 
c                sufficient space to retain the numrem integrals
c                not yet written.

       iptl=1
       iptv=iptl+
     .      3*4*forbyt(max(status%info(5,fhdle),status%info(5,fhdle2)))
       iptw=iptv+
     .      3*atebyt(max(status%info(5,fhdle),status%info(5,fhdle2)))
       maxint=2*(max(status%info(5,fhdle),status%info(5,fhdle2)))
       call sif_2e('init',name,work,work(iptv),maxint,numtot,
     .      fhdle,work(iptw),100000-iptw,lastrec,0)
       call sif_2e('init',name,work,work(iptv),numtot,numrem,
     .      fhdle2,work(iptw),100000-iptw,lastrec,last2)

       numrem=0
       numtot=0
       
 600   continue
       call sif_2e('read',name,work(iptl+forbyt(4*numrem)),
     .      work(iptv+numrem),maxint,
     .      numint,fhdle,work(iptw),100000-iptw,lastrec,0)

       write(6,*) 'read ',numint,' 2e integrals, lastrec=',lastrec
       numtot=numtot+numint
       last2=0
       if(lastrec) last2=2
       call sif_2e('writ',name,work,work(iptv),numint+numrem,
     .      numrem,fhdle2,work(iptw)
     .     ,100000-iptw,lastrec,last2)

c      call printintegrals(work,work(20001),numint)
         

       if (.not. lastrec) goto 600
 605   continue


        write(6,*) 'total number of 2e integrals:',numtot
       


     
       call sif_finalize(fhdle)
       call sif_finalize(fhdle2)

       stop
       end


        subroutine printintegrals(labels,values,numint)
        implicit none
        integer numint,labels(4,numint)
        real*8 values(*)
        integer i

         do i=1, numint
           write(6,100) labels(1:4,i),values(i)
         enddo
 100     format('(',i4,',',i4,'|',i4,',',i4,')',3x,f14.8)
         return
         end


        
