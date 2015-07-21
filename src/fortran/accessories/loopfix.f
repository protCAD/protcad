c234567
      program loopfix

      character*80 linebuffer
	  character*80 bf1
	  character*80 bf2
      character*80 filein
      integer numargs
   
10    format (a80)

      numargs = iargc()
      if (numargs.ne.1) then
        print *,'Usage: loopfix filename'
        goto 999
      end if
      call getarg(1,filein)

      open (unit=1,file=filein)
      open (unit=2,file='out.pdb')
      
100   read (1,10,end=900) linebuffer
      if (linebuffer(14:15).eq.'C ') then
        bf1 = linebuffer
        goto 100
      end if
      if (linebuffer(14:15).eq.'O ') then
        bf2 = linebuffer
        goto 100
      end if
      write (2,10) linebuffer 
      if (linebuffer(14:15).eq.'CA') then
        write (2,10) bf1
        write (2,10) bf2
        bf1 = ''
        bf2 = ''
      end if
      goto 100

900   close(1)
      close(2)
999   end
      
     
