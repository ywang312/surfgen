       type statustype
        integer unitno(2,10)
        integer nbft(10)
        integer nsym(10)
        integer nbpsy(8,10)
        integer info(7,10)
       end type statustype
       type(statustype) :: status

       common /sif_data/ status

