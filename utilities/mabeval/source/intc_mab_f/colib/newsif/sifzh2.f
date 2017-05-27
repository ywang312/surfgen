      subroutine sifzh2(
     & aoints,  ntitle,  nsym,    nbft,
     & ninfo,   nenrgy,  nmap,    title,
     & nbpsy,   slabel,  info,    bfnlab,
     & ietype,  energy,  imtype,  map,
     & mapdim,  ierr )
c
c  low-level header2 reading routine.
c
c  *** this routine should not be called directly by user programs. ***
c
c  see sifrh2() for argument description.
c
c  if ( nmap.eq.0 ) then
c     mapdim=1 ;imtype(*) and map(*) are not referenced.
c  else
c     mapdim=nmap ;imtype(*) and map(*) are read.
c  endif
c
c  26-jun-89 written by ron shepard.
c
       implicit logical(a-z)
      integer  aoints, ntitle, nsym,   nbft,   ninfo,
     & nenrgy, nmap,   mapdim, ierr
      character*80     title(ntitle)
      integer  nbpsy(nsym)
      character*4      slabel(nsym)
      integer  info(ninfo)
      character*8      bfnlab(nbft)
      integer  ietype(nenrgy)
      real*8   energy(nenrgy)
      integer  imtype(mapdim)
      integer  map(nbft,mapdim)
c
      if ( nmap.eq.0 ) then
         read(aoints,iostat=ierr)
     &    title,nbpsy,slabel,info,bfnlab,ietype,energy
      else
         read(aoints,iostat=ierr)
     &    title,nbpsy,slabel,info,bfnlab,ietype,energy,imtype,map
      endif
c
      return
      end
