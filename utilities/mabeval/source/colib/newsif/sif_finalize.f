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

