      program testcod34x

      implicit none
      real*8 buffer(10)
      integer intblk,stn121,stn2,stn3,nintkl,
     .    ksym,k1s,k1f,l1,iscr
      integer bits(11),bits4x(11)
      integer n,odata(11),inpdata(11),i,j,k
      character*6 onames(11),inames(11)
      character*6 onames4x(11),inames4x(11)
      data bits /18,18,18,18,18,5,10,10,10,3,0/
      data bits4x /20,20,20,20,4,4,10,10,10,10,20/
      data inames /'intblk','stn121','stn2  ','stn3  ',
     .             'nintkl','ksym  ','k1s   ','k1f   ',
     .             'l1    ','iscr  ','      '/

      data onames /'iblsym','ksym  ','k1s   ','k1f   ',
     .             'l1bl  ','intblk','stn121','stn3  ',
     .             'stn1  ','stn2  ','nintkl'/
      data inames4x /'strtnf','stn121','stn3  ','nintkl',
     .               'lsym  ','ksym  ','l1s   ','l1f   ',
     .               'k1s   ','k1f   ','k1l1  '/
      data onames4x /'intblk','stn121','stn3  ','nintkl',
     .             'lsym  ','ksym  ','l1s   ','l1f   ',
     .             'k1s   ','k1f   ','k1l1  '/

      write(6,*) 'New release has changed fil3w,fil3x header records'
      write(6,*) '** Testing coding/decoding for 3external integrals **'
      write(6,*) '** Using two 64bit words for record header         **'

      write(6,100) (inames(i),bits(i),i=1,11)
 100   format('current bit range for input variables:'/,11(a6,':',i10/))

        do n=3,0,-1
        write(6,*) '-- Testing maximum bit range minus',n
        do i=1,10
         inpdata(i)=2**(max(bits(i)-n,1)-1)-2
         inpdata(i)=max(inpdata(i),0)
         write(6,*) 'inpdata(',i,')=',inpdata(i),'bits=',bits(i)
        enddo
         inpdata(11)=0
         inpdata(6)=min(inpdata(6),8)

        call code3x( buffer, inpdata(1),inpdata(2),
     .    inpdata(3),inpdata(4),inpdata(5),inpdata(6),
     .    inpdata(7),inpdata(8),inpdata(9),inpdata(10))
       call dcod3x(buffer,odata(1),odata(2),odata(3),
     .             odata(4),odata(5),odata(6),odata(7),
     .             odata(8),odata(9),odata(10),odata(11))

 150   format(/,11(2(a6,2x,i8,2x)/))
       write(6, *) '   INPUT DATA                  OUTPUT DATA    '
       write(6,150)(inames(i),inpdata(i),onames(i),odata(i),i=1,11)

       enddo 

      write(6,*) '*****************************************************'

      write(6,*) 'New release has changed fil4w,fil4x header records'
      write(6,*) '** Testing coding/decoding for 4external integrals **'
      write(6,*) '** Using three 64bit words for record header       **'

      write(6,100) (inames4x(i),bits4x(i),i=1,11)


        do n=3,0,-1
        write(6,*) '-- Testing maximum bit range -',n
        do i=1,11
         inpdata(i)=2**(max(bits4x(i)-n,1)-1)-2
         inpdata(i)=max(inpdata(i),1)
        enddo
        do i=7,10
         inpdata(i)=2**(max(bits4x(i)-n,1)-1)-2
         inpdata(i)=max(inpdata(i),1)
        enddo
         inpdata(6)=min(inpdata(6),8)
         inpdata(5)=min(inpdata(5),8)

        call code4x( buffer, inpdata(1),inpdata(2),
     .    inpdata(3),inpdata(4),inpdata(5),inpdata(6),
     .    inpdata(7),inpdata(8),inpdata(9),inpdata(10),inpdata(11))
       call dcod4x(buffer,odata(1),odata(2),odata(3),
     .             odata(4),odata(5),odata(6),odata(7),
     .             odata(8),odata(9),odata(10),odata(11))

       write(6, *) '   INPUT DATA                  OUTPUT DATA    '
       write(6,150)(inames4x(i),inpdata(i),onames4x(i),odata(i),i=1,11)

       enddo



       stop     
       end
