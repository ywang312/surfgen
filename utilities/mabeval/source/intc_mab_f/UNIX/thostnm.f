      character*64 name
      write(*,*) 'input imax'
      read(*,*) imax
      ierr = hostnm( name(:imax) )
      write(*,*) 'name=', name(:imax)
      stop
      end
