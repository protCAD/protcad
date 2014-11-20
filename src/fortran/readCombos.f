      subroutine ReadCombos (comboTypes,comboTypesComp,contactdir)

      include 'parameters.h'

      integer comboTypes(maxNumCoord),
     +   comboTypesComp(maxNumCoord,570,maxnumatomtypesEnv)
      integer counter
      character*(*)  contactdir
      character*11  filename
      character*120 concat1
      character*120 concat2
      logical exist

10    format (a11)
11    format (15x,i9)
13    format (i13,i12,i12,i12)

      call fileexist('combo0-13list.dat',contactdir,concat1,exist)
      if (exist) then
       open (unit=2,file=concat1)
      else
       print *,'Error opening: ',concat1
       print *,'Stopping'
       stop
      end if

      do counter = 1,14 
       read (2,10) filename
       call fileexist(filename,contactdir,concat2,exist)
       if (exist) then
        open (unit=3,file=concat2)
        read (3,11) comboTypes(counter)
        do i=1,comboTypes(counter)
         read (3,13) (comboTypesComp(counter,i,j),
     &    j=1,maxnumatomtypesEnv)
        end do
       else
        print *,'error reading: ',concat2
       end if
      end do
      close (2)
      comboTypes(15) = 1

      close (15)
      end 
