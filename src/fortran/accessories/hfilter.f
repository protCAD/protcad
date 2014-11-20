c234567
      program filterhydrogens

      character*80 linebuffer
	  character*80 bf1
	  character*80 bf2
      character*80 filein
      character*80 fileout
      integer numargs
   
10    format (a80)

      numargs = iargc()
      if (numargs.ne.2) then
        print *,'Usage: hfilter infilename outfilename'
        goto 999
      end if
      call getarg(1,filein)
      call getarg(2,fileout)

      open (unit=1,file=filein)
      open (unit=2,file=fileout)
      
100   read (1,10,end=900) linebuffer
      if (linebuffer(14:14).eq.'H'.or.
     &    linebuffer(13:13).eq.'H') then
        goto 100
      end if
      write (2,10) linebuffer 
      goto 100

900   close(1)
      close(2)
999   end
      
     
