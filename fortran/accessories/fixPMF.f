c----------------------------------------------------------------
      program fixPMF

      parameter (maxnumatomtypes=20)
      parameter (pmfMaxNumBins=40)

      character*80  infile
      character*80  outfile
      real PMF(maxnumatomtypes,
     &  maxnumatomtypes,pmfMaxNumBins)
      integer  numbins,numatomtypes1,numatomtypes2
      integer  junk1,junk2

10      format (i2,1x,i2,1x,f7.3,1x,f7.3,1x,f7.3,1x,f7.3,1x,
     &  f7.3,1x,f7.3,1x,f7.3,1x,f7.3,1x,f7.3,1x,f7.3,1x,f7.3,1x,
     &  f7.3,1x,f7.3,1x,f7.3,1x,f7.3,1x,f7.3,1x,f7.3,1x,f7.3,1x,
     &  f7.3,1x,f7.3,1x,f7.3,1x,f7.3,1x,f7.3,1x,f7.3,1x,f7.3,1x,
     &  f7.3,1x,f7.3,1x,f7.3,1x,f7.3,1x,f7.3,1x,f7.3,1x,f7.3,1x,
     &  f7.3,1x,f7.3,1x,f7.3,1x,f7.3,1x,f7.3,1x,f7.3,1x,f7.3,1x,
     &  f7.3,1x,f7.3,1x,f7.3,1x,f7.3,1x,f7.3,1x,f7.3,1x,f7.3,1x,
     &  f7.3,1x,f7.3,1x,f7.3,1x,f7.3,1x,f7.3,1x,f7.3,1x,f7.3,1x,
     &  f7.3,1x,f7.3,1x,f7.3,1x,f7.3,1x,f7.3,1x,f7.3,1x,f7.3)

      numbins = 23
      numatomtypes1 = 20
      numatomtypes2 = 20

      infile = 'PMF.in'
      outfile = 'PMF.out'

      do i=1,numatomtypes1
       do j=1,numatomtypes2
        do k=1,numbins
         PMF(i,j,k) = 0.0
        end do
       end do
      end do

      open (unit=21,file=infile)
      do i=1,numatomtypes1
       do j=i,numatomtypes2
        read (21,10) junk1,junk2,(PMF(i,j,k),k=1,numbins)
        do k=1,numbins
         PMF(j,i,k) = PMF(i,j,k)
        end do
       end do
      end do
      close(21)

      open (unit=22,file=outfile)
      do i=1,numatomtypes1
       do j=1,numatomtypes2
        write(22,10) i-1,j-1,(PMF(i,j,k),k=1,numbins)
       end do
      end do
      close(22)
      end
