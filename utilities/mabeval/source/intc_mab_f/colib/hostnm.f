      integer function hostnm( hostn )
c
      implicit       character*(*) hostn
c
mdc*if rs6000 .and. (.not. extname)
c    # interface to the correct fortran external.
      integer  hostnm_
      external hostnm_
      hostnm = hostnm_( hostn )
mdc*else
c
c     # interface to the c_hostnm() c function.
c
      integer  c_hostnm
      external c_hostnm
c
c     # passing character strings directly into c functions
c     # is complicated because the length is concealed
c     # in different ways by different compilers.  this
c     # interface simplifies this by adding an explicit,
c     # and redundant, length argument.
c
      hostn  = ' '
      hostnm = c_hostnm( hostn, len(hostn) )
c
mdc*endif
      return
      end
