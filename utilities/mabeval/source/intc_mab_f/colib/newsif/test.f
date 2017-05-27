       program testsif

       character*8 name
       integer i,itypea,itypeb,ierr
       integer nsym,nbpsy(8),recordlen
       integer ietype(3),ntitle,nmap,map(2,512)
       character*80 title(5)
       real*8 d1(100),fcore,energy(3),work(100000)
       character*4 slabel(8)
       character*8 bfnlab(512)





       do itypea=-5,5
        do itypeb=-50,50 
         name='--------'
         call siftyp(itypea,itypeb,name,ierr)
         if (ierr.eq.0) 
     .   write(6,'(a,2x,i4,i4,a,a8,2x,i4)') 
     .      'defined itypea,itypeb=',itypea,itypeb,'=>',name
        enddo
       enddo

       slabel(1)=' a1g'
       slabel(2)=' b1g'
       slabel(3)=' b2g'
       slabel(4)=' b3g'
       slabel(5)=' a2 '
       slabel(6)=' b1 '
       slabel(7)=' b2 '
       slabel(8)=' b  '
       nsym=8
       nbpsy(1)=10
       nbpsy(2)=3
       nbpsy(3)=5
       nbpsy(4)=10
       nbpsy(5)=1
       nbpsy(6)=0
       nbpsy(7)=8
       nbpsy(8)=7
       nbft=0
       do i=1,8
          nbft=nbft+nbpsy(i)
       enddo

       write(6,*) '-------------------------- SAO mode ---------------'
       write(6,'(8(a4,1x)/8(i4,1x))') slabel(1:8),nbpsy(1:8)
       call  sif_baslab(nsym,nbpsy,slabel,bfnlab,'SAO')
       do i=1,nbft
         write(6,'(i8,3x,a8)') i,bfnlab(i)
       enddo  
         
       write(6,*) '-------------------------- MO mode ---------------'
       write(6,'(8(a4,1x)/8(i4,1x))') slabel(1:8),nbpsy(1:8)
       call  sif_baslab(nsym,nbpsy,slabel,bfnlab,'MO ')
       do i=1,nbft
         write(6,'(i8,3x,a8)') i,bfnlab(i)
       enddo  
         
       write(6,*) '-------------------------- SYM mode ---------------'
       write(6,'(8(a4,1x)/8(i4,1x))') slabel(1:8),nbpsy(1:8)
       call  sif_baslab(nsym,nbpsy,slabel,bfnlab,'SYM')
       do i=1,nbft
         write(6,'(i8,3x,a8)') i,bfnlab(i)
       enddo  
         
       


       stop
       end



        
