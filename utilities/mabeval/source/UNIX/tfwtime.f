c     # test of elapsed wall-clock timer, fwtime().
      real*8 second, fwtime
      external fwtime
      do i = 1, 10
         second = fwtime()
         write(*,*) 'second=', second
         write(*,*) 'enter <cr> to continue...'
         read(*,*)
      enddo
      stop
      end
