c234567
      program scorefix

      character*200 linebuffer
      character*200 bf1
      character*200 bf2
      character*200 filein
      character*200 fileout
      integer numargs
   
10    format (a200)

      numargs = iargc()
      if (numargs.ne.2) then
        print *,'Usage: scorefix filein fileout'
        goto 999
      end if
      call getarg(1,filein)
      call getarg(2,fileout)

      open (unit=1,file=filein)
      open (unit=2,file=fileout)
      
100   read (1,10,end=900) linebuffer
      if (linebuffer(1:6).ne.'/usr4/') then
        goto 100
      end if
      write (2,10) linebuffer 
      goto 100

900   close(1)
      close(2)
999   end
      
     
