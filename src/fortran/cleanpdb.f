      program cleanpdb

c Author: C. Summa
c-------------------------------------------------------
c Variable declarations


      include 'parameters.h'
      include 'pdb.h'

      character*80 filename
      character*80 outfile
      character*2  directory
      integer      numarguments

      numarguments = iargc()
      if (numarguments.ne.2) then
        print *,'Usage: cleanpdb pdbname outfilename'
        goto 999
      endif
      directory = './'

      call getarg(1,filename)
      call getarg(2,outfile)

      print *,'reading pdbfile : ',filename

      call pdbread_clean (filename,directory)
      call pdbwrite_w_ter (outfile)

999   end
