c234567
      program scorefix

      character*100 linebuffer
	  character*100 bf1
	  character*100 bf2
      character*100 filein
      character*100 fileout
      character*10  numstr
      character*100  tempstring
      real tempnum
      real num
      integer numargs
      integer begin
      integer finish
      integer tempend
      logical exist
      integer sumcounter
      integer lowercounter
   
10    format (a100)

      numargs = iargc()
      if (numargs.ne.2) then
        print *,'Usage: scorefix filein number'
        goto 999
      end if
      call getarg(1,filein)
      call getarg(2,numstr)
      call getthereal(num,numstr)

      open (unit=1,file=filein)
      sumcounter = 0
      lowercounter = 0 
100   read (1,10,end=900) linebuffer
      
      call findNextWord(linebuffer,1,begin,finish,exist)
      if (exist) then
      	tempend = finish+1
        tempstring = linebuffer(begin:finish)
        call findNextWord(linebuffer,tempend,begin,finish,exist)
        if (exist) then
          tempend = finish+1
          tempstring = linebuffer(begin:finish)
          call findNextWord(linebuffer,tempend,begin,finish,exist)
          if (exist) then
            tempstring = linebuffer(begin:finish)
          end if
        end if
        call getthereal(tempnum,tempstring)
        if (tempnum.ge.num) then
         lowercounter = lowercounter + 1
        endif
      end if
      sumcounter = sumcounter + 1
      goto 100

900   close(1)
      print *,'total number read = ',sumcounter
      print *,'total number correctly scored = ',lowercounter
999   end
     
