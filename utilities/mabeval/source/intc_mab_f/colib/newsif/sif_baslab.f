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

