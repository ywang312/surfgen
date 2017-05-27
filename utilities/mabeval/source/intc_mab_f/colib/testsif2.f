      program testsif
       implicit none
       integer maxreclen,nbft,ifmt,lrec,nmax,ierr
       integer i,usage
       integer basereclen
       integer label(32*34),labelout(32*34) 
       real*8  buffer(32*30)


       write(6,*) ' testing sifcfg - 10bit compression'

       nbft=300
       basereclen=4096
       
       do i=1,8
       call sifcfg(1,basereclen*i,nbft,0,ifmt,lrec,nmax,ierr)
       usage=nmax+(nmax/16)*5+2
       write(6,50) basereclen*i,nbft,ifmt,lrec,nmax,ierr,usage
       if (mod(nmax,16).ne.0) write(6,*) 'failure nmax blocking'
       if (usage.gt.lrec) 
     .    write(6,*) 'failure lrec',usage,'vs.',lrec      

       enddo
 50    format('1e: basereclen=',i5,' nbft=',i4,' ifmt=',i2,
     .           ' lrec=',i6,' nmax=',i6,' ierr=',i2,' usage=',i6)
 

       do i=1,8
       call sifcfg(2,basereclen*i,nbft,0,ifmt,lrec,nmax,ierr)
       usage=nmax+(nmax/8)*5+1
       write(6,51) basereclen*i,nbft,ifmt,lrec,nmax,ierr,usage
       if (mod(nmax,8).ne.0) write(6,*) 'failure nmax blocking'
       if (usage.gt.lrec)         
     .    write(6,*) 'failure lrec',usage,'vs.',lrec      
       enddo
 51    format('2e: basereclen=',i5,' nbft=',i4,' ifmt=',i2,
     .           ' lrec=',i6,' nmax=',i6,' ierr=',i2,' usage=',i6)
 


       write(6,*) ' Testing plab40/ulab40 compression multiple 8 '

       do i=1,32*32
         label(i)=mod(i,1024)
         label(i)=label(i)*label(i)*label(i)
       enddo
       call plab40(buffer,label,32*32)
      
c      call printbin(buffer,4) 

       call ulab40(buffer,labelout,32*32)
     
c     do i=1,10
c        write(*,56) i, label(i),labelout(i)
c     enddo

       ierr=0
       do i=1,32*32
        if (label(i).ne.labelout(i)) then
         write(*,55) i,label(i),labelout(i)
         ierr=ierr+1
        endif
       enddo
 55    format('ERROR i=',i8,' in=',i8, ' out=',i8)
 56    format('COMPARE i=',i8,' in=',i8, ' out=',i8)
       if (ierr.eq.0) write(6,*) 'no errors found'

       write(6,*) ' Testing plab40/ulab40 compression arb. length '

       do i=1,32*31+10
         label(i)=mod(2*i,1024)
       enddo
       call plab40(buffer,label,32*31+10)

c      call printbin(buffer,4) 

       call ulab40(buffer,labelout,32*31+10)

c     do i=1,10
c        write(*,56) i, label(i),labelout(i)
c     enddo

       ierr=0
       do i=1,32*31+10
        if (label(i).ne.labelout(i)) then
         write(*,55) i,label(i),labelout(i)
         ierr=ierr+1
        endif
       enddo
       if (ierr.eq.0) write(6,*) 'no errors found'

       stop     
       end




        subroutine testsif2(basereclen,nbft)
        implicit none
        integer basereclen,nbft,ifmt,l1rec,n1max,ierr
        integer l2rec,n2max,iunit,ntitle,nsymao,norbt
        integer infao(6),nenrgy,nmap,nmbpsy(8),i,ietype(4)
        character*80 title(1)
        character*4 slabel(8)
        character*8 molab(512) 
        integer imtype(5),map(512)
        real*8 energy(8),buffer(basereclen),values(4097)
        integer labels(2,4097),num,labels2(4,4097)
        integer ibvtyp,itypea,itypeb,ioff,isym,j,icount,nrec
        integer iwait,reqnum,nipv

c
c       testing sifew2 sifrd2 
c
c
       call sifcfg(1,basereclen,nbft,0,ifmt,l1rec,n1max,ierr)

        iunit=10
        ntitle=1
        nsymao=2
        infao(1) = 1
        infao(2) = l1rec 
        infao(3) = n1max 
       call sifcfg(2,basereclen,nbft,0,ifmt,l2rec,n2max,ierr)
        infao(4) = l2rec 
        infao(5) = n2max 
        infao(6) = ifmt 
        nenrgy=1
        nmap=0
        title(1)="TESTSIF"
        nmbpsy(1)=60
        nmbpsy(2)=60
         norbt=0
        do i=1,nsymao
          norbt=norbt+nmbpsy(i)
        enddo
 
        do i=1,norbt
         write(molab(i),'(a4,i4)') 'MO__',i
        enddo
        do i=1,nsymao
          write(slabel(i),'(a3,i1)') 'SYM',i
        enddo
        ietype(1)=-1026
        energy(1)=-555.55d0
        imtype(1)=0 
        map(1)=1

        open(unit=iunit, file='moints.test',form='unformatted',
     .       access='sequential')
        call sifwh(iunit,ntitle,nsymao, norbt, 6    ,
     .             nenrgy, nmap, title, nmbpsy, slabel,
     .             infao, molab,ietype,energy,imtype,
     .             map,ierr)  
        itypea=0
        itypeb=1
        ibvtyp=0

        ioff=0
        icount=0
        do isym=1,nsymao
         do i=1,nmbpsy(isym)
           do j=1,i 
            icount=icount+1
            values(icount)=100.0d0+j*1.0d0+i*0.001d0
            labels(1,icount)=i+ioff
            labels(2,icount)=j+ioff
            if (icount.gt.n1max) then
              num=n1max
              call sifew1(iunit,infao,2,num,
     .                    0,itypea,itypeb,ibvtyp,
     .                   values,labels,999.0d0,
     .                   0,buffer,nrec,ierr)
             write(6,'(a,i6,a,i6)') 'X:wrote ',n1max,' integrals',icount
             values(1)=values(icount)
             labels(1,1)=labels(1,icount)
             labels(2,1)=labels(2,icount)
             icount=1
            endif
          enddo
         enddo
        ioff=ioff+nmbpsy(isym)
        enddo
        num=icount
        call sifew1(iunit,infao,2,num,
     .          2,itypea,itypeb,ibvtyp,
     .          values,labels,999.0d0,
     .          0,buffer,nrec,ierr)
        write(6,*) 'wrote ',icount,' integrals'



        itypea=3
        itypeb=0
        ibvtyp=0
        nipv =4

        ioff=0
        icount=0
        do isym=1,nsymao
         do i=1,nmbpsy(isym)
           do j=1,i
            icount=icount+1
            values(icount)=100.0d0+j*1.0d0+i*0.001d0
            labels2(1,icount)=i+ioff
            labels2(2,icount)=j+ioff
            labels2(3,icount)=i+ioff
            labels2(4,icount)=j+ioff
            if (icount.gt.n2max) then
              num=n2max
              call sifew2(iunit,infao,nipv,num,
     .                    0,itypea,itypeb,ibvtyp,
     .                   values,labels2,
     .                   0,buffer,iwait,nrec,reqnum,ierr)
             write(6,'(a,i6,a,i6)') 'Y:wrote ',n2max,' integrals',icount
             values(1)=values(icount)
             labels2(1,1)=labels2(1,icount)
             labels2(2,1)=labels2(2,icount)
             labels2(3,1)=labels2(3,icount)
             labels2(4,1)=labels2(4,icount)
             icount=1
            endif
          enddo
         enddo
        ioff=ioff+nmbpsy(isym)
        enddo
        num=icount
              call sifew2(iunit,infao,nipv,num,
     .                    2,itypea,itypeb,ibvtyp,
     .                   values,labels2,
     .                   0,buffer,iwait,nrec,reqnum,ierr)
        write(6,*) 'wrote ',icount,' integrals'

        close(iunit)
        return
        end



       
     




        subroutine printbin(buffer,n)
        implicit none
        integer*8 buffer(n),n
        integer i,j,joff
        integer word(64)

        joff=0
        do i=1,n
         joff=joff+1
          write(6,*) 'buffer(',joff,')=',buffer(joff)
          word=0
          do j=0,63
            word(64-j)=iand(ishft(buffer(joff),-j),1) 
          enddo
          write(6,'(7(10i1,1x))')word(1:64)
        enddo
        return
        end
 

