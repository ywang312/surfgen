       program testpacking
       use colib

       implicit none

c
c
c      test bit packing routines
c

c
c          32 bit integers <-> 64bit reals
c
       integer,parameter:: maxn=128

       integer  jj,ii,i,iarr(1:maxn)
       integer  ierr
       real*8     darr(1:maxn)
       logical  ll
       call little_endian(ll)
       write(6,*) 'little indian=',ll

        iarr(1)=21   
        iarr(2)=23
        iarr(3)=11264
        iarr(4)=2

        write(6,'(a,4i8)') 'input: ',iarr(1:4) 
        call plab16(darr,iarr,4)
        write(6,'(a,z16)') 'packed=',darr(1)
        iarr=0
        call ulab16(darr,iarr,4)

        write(6,'(a,4i8)') 'output: ',iarr(1:4) 



100    format(a,2(20i4/))

        ierr=0
        do ii=1,maxn
        write(6,*) 'testing sequence of ',ii, ' integers'
        do i=1,maxn
         iarr(i)=i+16500
         darr(i)=0.0d0
        enddo

c       write(6,100) 'Initial 16:',iarr(1:40)
        call plab16(darr,iarr,ii)
        iarr=0
        call ulab16(darr,iarr,ii)
c       write(6,100) 'Final 16  :',iarr(1:40)
         do jj=1,ii
          if (iarr(jj).ne.jj+16500) then 
           write(6,*) 'error got ',iarr(jj),' instead ' ,jj+16500
           ierr=ierr+1
          endif
         enddo
         do jj=ii+1,maxn
          if (iarr(jj).ne.0) then
              write(6,*) 'error2 at ',jj
           ierr=ierr+1
          endif 
         enddo
         enddo


        if (ierr.eq.0) write(6,*) ' plab16/ulab16 NO ERRORS FOUND'
         ierr=0

        do ii=1,maxn
        write(6,*) 'testing sequence of ',ii, ' integers'
        do i=1,maxn
         iarr(i)=i+65536
         darr(i)=0.0d0
        enddo

c       write(6,100) 'Initial 32:',iarr(1:40)
        call plab32(darr,iarr,ii)
        iarr=0
        call ulab32(darr,iarr,ii)
c       write(6,100) 'Final 32  :',iarr(1:40)
         do jj=1,ii
          if (iarr(jj).ne.jj+65536) then  
           write(6,*) 'error at ',jj
           ierr=ierr+1
          endif
         enddo
         do jj=ii+1,maxn
          if (iarr(jj).ne.0) then
              write(6,*) 'error2 at ',jj
           ierr=ierr+1
          endif
         enddo

         enddo

        if (ierr.eq.0) write(6,*) ' plab32/ulab32 NO ERRORS FOUND'

         ierr=0

        do ii=1,maxn
        write(6,*) 'testing sequence of ',ii, ' integers'
        do i=1,maxn
         iarr(i)=i+120
         darr(i)=0.0d0
        enddo

c       write(6,100) 'Initial 8:',iarr(1:40)
        call plab8(darr,iarr,ii)
        iarr=0
        call ulab8(darr,iarr,ii)
c       write(6,100) 'Final 8  :',iarr(1:40)
         do jj=1,ii
          if (iarr(jj).ne.jj+120) then
           write(6,*) 'error got',iarr(jj), 'instead of ',jj+120
           ierr=ierr+1
          endif
         enddo
         do jj=ii+1,maxn
          if (iarr(jj).ne.0) then
              write(6,*) 'error2 at ',jj
           ierr=ierr+1
          endif
         enddo

         enddo

        if (ierr.eq.0) write(6,*) ' plab8/ulab8 NO ERRORS FOUND'


         ierr=0

        do ii=1,maxn
        write(6,*) 'testing sequence of ',ii, ' integers'
        do i=1,maxn
         iarr(i)=mod(i,2)
         darr(i)=0.0d0
        enddo

c       write(6,100) 'Initial 1:',iarr(1:40)
        call plab1(darr,iarr,ii)
        iarr=0
        call ulab1(darr,iarr,ii)
c       write(6,100) 'Final 1  :',iarr(1:40)
         do jj=1,ii
          if (iarr(jj).ne.mod(jj,2)) then
           write(6,*) 'error got',iarr(jj), 'instead of ',mod(jj,2)
           ierr=ierr+1
          endif
         enddo
         do jj=ii+1,maxn
          if (iarr(jj).ne.0) then
              write(6,*) 'error2 at ',jj
           ierr=ierr+1
          endif
         enddo

         enddo

        if (ierr.eq.0) write(6,*) ' plab1/ulab1 NO ERRORS FOUND'




       stop
       end

