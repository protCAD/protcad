      subroutine fileexist(fname,dir,concat,exist)

      character*(*) fname
      character*(*) dir
      character*120 tstring
      character*120 concat

      integer fnstart
      integer fnend
      integer dirstart
      integer dirend
      logical wordexist
      logical exist

      call findNextWord(dir,1,dirstart,dirend,wordexist)
      if (.not.wordexist) then
        exist = .false.
        return
      end if

      call findNextWord(fname,1,fnstart,fnend,wordexist)
      if (.not.wordexist) then
       exist = .false.
       return
      end if

      do i=1,120
       tstring(i:i) = ''
      end do

      if (dir.ne.'null') then
        tstring(1:(dirend-dirstart)+1) = dir(dirstart:dirend)
        tstring((dirend-dirstart)+2:
     &        (dirend-dirstart)+2+(fnend-fnstart)+1)=
     &        fname(fnstart:fnend)
        concat = tstring
      else
        tstring(1:(fnend-fnstart)+1) = fname(fnstart:fnend)
        concat = tstring
      end if

      open (unit=1,err = 995,file=concat,status='OLD')
       exist = .true.
      close (1)
      goto 999

995   exist = .false. 
      return

999   return 
      end

c----------------------------------------------------------------
      subroutine writeenergy(energyfile,EwordAlternate,
     & Eword,singleBodyEnergy,comboTypes,fPairwise,
     & PMFENERGY,EPair,contactdir)

      include 'parameters.h'
      character*(*)  energyfile
      character*(*)  contactdir
      character*120  concat
      real Eword(maxnumatomtypes,maxNumCoord,570)
      real EwordAlternate(maxnumatomtypes,maxNumCoord,570)
      real EPair(maxnumatomtypes,maxnumatomtypesEnv)
      real singleBodyEnergy(maxnumatomtypes,binsInHistogram)
      real fPairwise(maxnumatomtypes,maxnumatomtypes)
      real PMFENERGY(maxnumatomtypes,
     &  maxnumatomtypesEnv,pmfMaxNumBins)
      logical exist

      integer  comboTypes(maxNumCoord)

      call fileexist(energyfile,contactdir,concat,exist)

      open (unit=21,file=concat,status='NEW',form='unformatted')

      do i=1,maxnumatomtypes
       do j=1,maxNumCoord
        write (21) (Eword(i,j,k),k=1,comboTypes(j))
        write (21) (EwordAlternate(i,j,k),k=1,comboTypes(j))
       end do
        write (21) (singleBodyEnergy(i,k),k=1,binsInHistogram)
      end do
      do i=1,maxnumatomtypes
        write (21) (fPairwise(i,j),j=1,maxnumatomtypes)
      end do
      do i=1,maxnumatomtypes
       do j=1,maxnumatomtypesEnv
        write (21) EPair(i,j)
       end do
      end do
      do i=1,maxnumatomtypes
       do j=1,maxnumatomtypesEnv
        do k=1,pmfMaxNumBins
         write (21) PMFENERGY(i,j,k)
        end do
       end do
      end do
      close (21)
      end

c---------------------------------------------------
      subroutine readenergy(filename,EwordAlternate,
     & Eword,singleBodyEnergy,comboTypes,fPairwise,
     & PMFENERGY,EPair,fdir)

      include 'parameters.h'

      character*(*)  filename
      character*(*)  fdir
      real Eword(maxnumatomtypes,maxNumCoord,570)
      real EwordAlternate(maxnumatomtypes,maxNumCoord,570)
      real singleBodyEnergy(maxnumatomtypes,binsInHistogram)
      integer  comboTypes(maxNumCoord)
      real fPairwise(maxnumatomtypes,maxnumatomtypes)
      real PMFENERGY(maxnumatomtypes,maxnumatomtypesEnv,pmfMaxNumBins)
      real    EPair(maxnumatomtypes,maxnumatomtypesEnv)
      character*120 concat
      integer i,j,k
      logical exist

      call fileexist(filename,fdir,concat,exist)
      if (.not.exist) then
       print *,'energyfile = ',filename
       print *,'contactdir = ',fdir
       print *,concat
       print *,'Error opening energyfile'
       print *,'Stopping'
       stop
      end if

      open (unit=21,file=concat,status='OLD',form='unformatted')

      print *,'im in readenergy'
      do i=1,maxnumatomtypes
       do j=1,maxNumCoord
        read (21) (Eword(i,j,k),k=1,comboTypes(j))
        read (21) (EwordAlternate(i,j,k),k=1,comboTypes(j))
       end do
        read (21) (singleBodyEnergy(i,k),k=1,binsInHistogram)
      end do
      do i=1,maxnumatomtypes
        read (21) (fPairwise(i,j),j=1,maxnumatomtypes)
      end do
      do i=1,maxnumatomtypes
       do j=1,maxnumatomtypesEnv
        read (21) EPair(i,j)
       end do
      end do
      do i=1,maxnumatomtypes
       do j=1,maxnumatomtypesEnv
        do k=1,pmfMaxNumBins
         read (21) PMFENERGY(i,j,k)
        end do
       end do
      end do

999   close (21)
      end
c----------------------------------------------------------------
      subroutine writeenergy2(energyfile,
     & Eword,singleBodyEnergy,comboTypes,fPairwise,
     & PMFENERGY,EPair,contactdir)

      include 'parameters.h'
      character*(*)  energyfile
      character*(*)  contactdir
      character*120  concat
      real Eword(maxnumatomtypes,maxNumCoord,570)
      real EPair(maxnumatomtypes,maxnumatomtypesEnv)
      real singleBodyEnergy(maxnumatomtypes,binsInHistogram)
      real fPairwise(maxnumatomtypes,maxnumatomtypes)
      real PMFENERGY(maxnumatomtypes,
     &  maxnumatomtypesEnv,pmfMaxNumBins)
      logical exist

      integer  comboTypes(maxNumCoord)

      call fileexist(energyfile,contactdir,concat,exist)

      open (unit=21,file=concat,status='NEW',form='unformatted')

      do i=1,maxnumatomtypes
       do j=1,maxNumCoord
        write (21) (Eword(i,j,k),k=1,comboTypes(j))
       end do
        write (21) (singleBodyEnergy(i,k),k=1,binsInHistogram)
      end do
      do i=1,maxnumatomtypes
        write (21) (fPairwise(i,j),j=1,maxnumatomtypes)
      end do
      do i=1,maxnumatomtypes
       do j=1,maxnumatomtypesEnv
        write (21) EPair(i,j)
       end do
      end do
      do i=1,maxnumatomtypes
       do j=1,maxnumatomtypesEnv
        do k=1,pmfMaxNumBins
         write (21) PMFENERGY(i,j,k)
        end do
       end do
      end do
      close (21)
      end

c---------------------------------------------------
      subroutine readenergy2(filename,
     & Eword,singleBodyEnergy,comboTypes,fPairwise,
     & PMFENERGY,EPair,fdir)

      include 'parameters.h'

      character*(*)  filename
      character*(*)  fdir
      real Eword(maxnumatomtypes,maxNumCoord,570)
      real singleBodyEnergy(maxnumatomtypes,binsInHistogram)
      integer  comboTypes(maxNumCoord)
      real fPairwise(maxnumatomtypes,maxnumatomtypes)
      real PMFENERGY(maxnumatomtypes,maxnumatomtypesEnv,pmfMaxNumBins)
      real    EPair(maxnumatomtypes,maxnumatomtypesEnv)
      character*120 concat
      integer i,j,k
      logical exist

      call fileexist(filename,fdir,concat,exist)
      if (.not.exist) then
       print *,'energyfile = ',filename
       print *,'contactdir = ',fdir
       print *,concat
       print *,'Error opening energyfile'
       print *,'Stopping'
       stop
      end if

      open (unit=21,file=concat,status='OLD',form='unformatted')

      print *,'im in readenergy'
      do i=1,maxnumatomtypes
       do j=1,maxNumCoord
        read (21) (Eword(i,j,k),k=1,comboTypes(j))
       end do
        read (21) (singleBodyEnergy(i,k),k=1,binsInHistogram)
      end do
      do i=1,maxnumatomtypes
        read (21) (fPairwise(i,j),j=1,maxnumatomtypes)
      end do
      do i=1,maxnumatomtypes
       do j=1,maxnumatomtypesEnv
        read (21) EPair(i,j)
       end do
      end do
      do i=1,maxnumatomtypes
       do j=1,maxnumatomtypesEnv
        do k=1,pmfMaxNumBins
         read (21) PMFENERGY(i,j,k)
        end do
       end do
      end do

999   close (21)
      end
c----------------------------------------------------------------
      subroutine writePMF(energyfile,
     & PMF,numatomtypes1,numatomtypes2,numbins)

      include 'parameters.h'
      character*(*)  energyfile
      real PMF(maxnumatomtypes,
     &  maxnumatomtypes,pmfMaxNumBins)
      integer  numbins,numatomtypes1,numatomtypes2

10      format (i2,1x,i2,1x,f7.3,1x,f7.3,1x,f7.3,1x,f7.3,1x,
     &  f7.3,1x,f7.3,1x,f7.3,1x,f7.3,1x,f7.3,1x,f7.3,1x,f7.3,1x,
     &  f7.3,1x,f7.3,1x,f7.3,1x,f7.3,1x,f7.3,1x,f7.3,1x,f7.3,1x,
     &  f7.3,1x,f7.3,1x,f7.3,1x,f7.3,1x,f7.3,1x,f7.3,1x,f7.3,1x,
     &  f7.3,1x,f7.3,1x,f7.3,1x,f7.3,1x,f7.3,1x,f7.3,1x,f7.3,1x,
     &  f7.3,1x,f7.3,1x,f7.3,1x,f7.3,1x,f7.3,1x,f7.3,1x,f7.3,1x,
     &  f7.3,1x,f7.3,1x,f7.3,1x,f7.3,1x,f7.3,1x,f7.3,1x,f7.3,1x,
     &  f7.3,1x,f7.3,1x,f7.3,1x,f7.3,1x,f7.3,1x,f7.3,1x,f7.3,1x,
     &  f7.3,1x,f7.3,1x,f7.3,1x,f7.3,1x,f7.3,1x,f7.3,1x,f7.3)

      open (unit=21,file=energyfile,status='NEW')

      do i=1,numatomtypes1
       if (i.ne.1) then
        do j=1,i-1
          write (21,10) (i-1),(j-1),(PMF(j,i,k),k=1,numbins)
        end do
       end if
       do j=i,numatomtypes2
        write (21,10) (i-1),(j-1),(PMF(i,j,k),k=1,numbins)
       end do
      end do
      close(21)
      end
c----------------------------------------------------------------
      subroutine writePMF_LoRes(energyfile,
     & PMF,numatomtypes1,numatomtypes2,numbins)

      include 'parameters.h'
      character*(*)  energyfile
      real PMF(maxnumatomtypes,
     &  maxnumatomtypesEnv,pmfMaxNumBins)
      integer  numbins,numatomtypes1,numatomtypes2

10      format (i2,1x,i2,1x,f7.3,1x,f7.3,1x,f7.3,1x,f7.3,1x,
     &  f7.3,1x,f7.3,1x,f7.3,1x,f7.3,1x,f7.3,1x,f7.3,1x,f7.3,1x,
     &  f7.3,1x,f7.3,1x,f7.3,1x,f7.3,1x,f7.3,1x,f7.3,1x,f7.3,1x,
     &  f7.3,1x,f7.3,1x,f7.3,1x,f7.3,1x,f7.3,1x,f7.3,1x,f7.3,1x,
     &  f7.3,1x,f7.3,1x,f7.3,1x,f7.3,1x,f7.3,1x,f7.3,1x,f7.3,1x,
     &  f7.3,1x,f7.3,1x,f7.3,1x,f7.3,1x,f7.3,1x,f7.3,1x,f7.3,1x,
     &  f7.3,1x,f7.3,1x,f7.3,1x,f7.3,1x,f7.3,1x,f7.3,1x,f7.3,1x,
     &  f7.3,1x,f7.3,1x,f7.3,1x,f7.3,1x,f7.3,1x,f7.3,1x,f7.3,1x,
     &  f7.3,1x,f7.3,1x,f7.3,1x,f7.3,1x,f7.3,1x,f7.3,1x,f7.3)

      open (unit=21,file=energyfile,status='NEW')

      do i=1,numatomtypes1
       do j=1,numatomtypes2
        write (21,10) (i-1),(j-1),(PMF(i,j,k),k=1,numbins)
       end do
      end do
      close(21)
      end

c----------------------------------------------------------------
      subroutine writeenergy2_formatted(energyfile,
     & Eword,singleBodyEnergy,comboTypes,fPairwise,
     & PMFENERGY,EPair,contactdir)

      include 'parameters.h'
      character*(*)  energyfile
      character*(*)  contactdir
      character*120  concat
      real Eword(maxnumatomtypes,maxNumCoord,570)
      real EPair(maxnumatomtypes,maxnumatomtypesEnv)
      real singleBodyEnergy(maxnumatomtypes,binsInHistogram)
      real fPairwise(maxnumatomtypes,maxnumatomtypes)
      real PMFENERGY(maxnumatomtypes,
     &  maxnumatomtypesEnv,pmfMaxNumBins)
      logical exist

      integer  comboTypes(maxNumCoord)

10    format (700(f12.4))
      call fileexist(energyfile,contactdir,concat,exist)

      open (unit=21,file=concat,status='NEW')

      do i=1,maxnumatomtypes
       do j=1,maxNumCoord
        write (21,10) (Eword(i,j,k),k=1,comboTypes(j))
       end do
	  end do
	  do i=1,maxnumatomtypes
        write (21,10) (singleBodyEnergy(i,k),k=1,binsInHistogram)
      end do
      do i=1,maxnumatomtypes
        write (21,10) (fPairwise(i,j),j=1,maxnumatomtypes)
      end do
      do i=1,maxnumatomtypes
        write (21,10) (EPair(i,j),j=1,maxnumatomtypesEnv)
      end do
      do i=1,maxnumatomtypes
       do j=1,maxnumatomtypesEnv
         write (21,10) (PMFENERGY(i,j,k),k=1,pmfMaxNumBins)
       end do
      end do
      close (21)
      end

c---------------------------------------------------
      subroutine readenergy2_formatted(filename,
     & Eword,singleBodyEnergy,comboTypes,fPairwise,
     & PMFENERGY,EPair,fdir)

      include 'parameters.h'

      character*(*)  filename
      character*(*)  fdir
      real Eword(maxnumatomtypes,maxNumCoord,570)
      real singleBodyEnergy(maxnumatomtypes,binsInHistogram)
      integer  comboTypes(maxNumCoord)
      real fPairwise(maxnumatomtypes,maxnumatomtypes)
      real PMFENERGY(maxnumatomtypes,maxnumatomtypesEnv,pmfMaxNumBins)
      real    EPair(maxnumatomtypes,maxnumatomtypesEnv)
      character*120 concat
      integer i,j,k
      logical exist

10    format (700(f12.4))
      call fileexist(filename,fdir,concat,exist)
      if (.not.exist) then
       print *,'energyfile = ',filename
       print *,'contactdir = ',fdir
       print *,concat
       print *,'Error opening energyfile'
       print *,'Stopping'
       stop
      end if

      open (unit=21,file=concat,status='OLD')

      print *,'im in readenergy'
      do i=1,maxnumatomtypes
       do j=1,maxNumCoord
        read (21,10) (Eword(i,j,k),k=1,comboTypes(j))
       end do
	  end do
	  do i=1,maxnumatomtypes
        read (21,10) (singleBodyEnergy(i,k),k=1,binsInHistogram)
      end do
      do i=1,maxnumatomtypes
        read (21,10) (fPairwise(i,j),j=1,maxnumatomtypes)
      end do
      do i=1,maxnumatomtypes
        read (21,10) (EPair(i,j),j=1,maxnumatomtypesEnv)
      end do
      do i=1,maxnumatomtypes
       do j=1,maxnumatomtypesEnv
         read (21,10) (PMFENERGY(i,j,k),k=1,pmfMaxNumBins)
       end do
      end do

999   close (21)
      end
c----------------------------------------------------------------
      subroutine writeenergy3_formatted(energyfile,
     & Eword,comboTypes,contactdir)

      include 'parameters.h'
      character*(*)  energyfile
      character*(*)  contactdir
      character*120  concat
      real Eword(maxnumatomtypes,maxNumCoord,570)
      logical exist

      integer  comboTypes(maxNumCoord)

10    format (700(f12.4))
c      call fileexist(energyfile,contactdir,concat,exist)

c      open (unit=21,file=concat,status='NEW')
      open (unit=21,file=energyfile,status='NEW')

      do i=1,maxnumatomtypes
       do j=1,maxNumCoord
        write (21,10) (Eword(i,j,k),k=1,comboTypes(j))
       end do
	  end do
      close (21)
      end

c---------------------------------------------------
      subroutine readenergy3_formatted(filename,
     & Eword,comboTypes,fdir)

      include 'parameters.h'

      character*(*)  filename
      character*(*)  fdir
      real Eword(maxnumatomtypes,maxNumCoord,570)
      integer  comboTypes(maxNumCoord)
      character*120 concat
      integer i,j,k
      logical exist

10    format (700(f12.4))
c      call fileexist(filename,fdir,concat,exist)
c      if (.not.exist) then
c       print *,'energyfile = ',filename
c       print *,'contactdir = ',fdir
c       print *,concat
c       print *,'Error opening energyfile'
c       print *,'Stopping'
c       stop
c      end if


c      open (unit=21,file=concat,status='OLD')
      open (unit=21,file=filename,status='OLD')
      print *,'im in readenergy'
      do i=1,maxnumatomtypes
       do j=1,maxNumCoord
	    print *,i,j
        read (21,10) (Eword(i,j,k),k=1,comboTypes(j))
       end do
	  end do
999   close (21)
      end
