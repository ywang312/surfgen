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
 
