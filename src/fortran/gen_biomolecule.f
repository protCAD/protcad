      program gen_biomolecule 

c Author: C. Summa
c-------------------------------------------------------
c Variable declarations


      include 'parameters.h'
      include 'pdb.h'

      character*80 filename
      character*2  directory
      integer      numarguments

      numarguments = iargc()
      if (numarguments.ne.1) then
        print *,'Usage: gen_biomolecule pdbname'
        goto 999
      endif
      directory = './'

      call getarg(1,filename)

      print *,'reading pdbfile : ',filename

      call pdbread_w_write (filename,directory)

999   end
