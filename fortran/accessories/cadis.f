      program cadis

c Author: C. Summa
c Last modified: 6/11/98

      include 'parameters.h'
      include 'pdb.h'

      real Distance
      real tempdistance
      real upperlimit
      real lowerlimit
      integer numarguments
      character*4 name
      character*40 pdbname
      character*80 databasedir
c------------------------------------------------------------
c Read in command line arguments

      numarguments = iargc()

      if (numarguments.ne.1) then
        print *,'Usage: cadis pdbname'
        goto 999
      endif
 
      upperlimit = 8.5
      lowerlimit = 0.0

      call getarg(1,pdbname)
      databasedir = 'null'
      call pdbread (pdbname,databasedir)
      name = ' CB ' 
 
        do i=1,numatm
         do j=i+1,numatm
          if (atomname(i).eq.name.and.atomname(j).eq.name) then
            tempdistance =  Distance(coord(1,i),coord(2,i),
     &       coord(3,i),coord(1,j),coord(2,j),coord(3,j)) 
            if (tempdistance.ge.lowerlimit.and.
     &          tempdistance.le.upperlimit.and.
c     &          resnum(j)-resnum(i).ge.30) then
     &           resnum(j)-resnum(i).ge.10) then
             print *,res(i),resnum(i),'  ',res(j),resnum(j),tempdistance
            end if
          end if
         end do
        end do

999     end
