      character*24 cdate, fdate
      external fdate
      cdate = fdate()
      write(*,*) 'cdate=', cdate
      stop
      end
