
      subroutine PMFmodule (pmfShellCounter,numatomtypes1,
     &   numatomtypes2,numdistances,PMFENERGY,filename)

c-----------------------------------------------------
c Author: C. Summa
c-----------------------------------------------------
c Variable declarations


      include 'parameters.h'

      character*(*) filename
      real PMFENERGY(maxnumatomtypes,maxnumatomtypes,pmfMaxNumBins)
      real pmfShellCounter(maxnumatomtypes,
     &    maxnumatomtypes,pmfMaxNumBins)
      real M(maxnumatomtypes,maxnumatomtypes)
      real RFO(maxnumatomtypes,maxnumatomtypes,pmfMaxNumBins)
      real FOALL
      real FXX(pmfMaxNumBins)
      real RFOALL(pmfMaxNumBins)
      real*8 numatmTotal2Real
      integer  count1,count2,count3,count4,count5,i,j,k
      integer  numatomtypes1
      integer  numatomtypes2
      integer  numdistances
      integer  follow1,follow2,follow3

10    format (i3,2x,i3,2x,f10.4,1x,f10.4,1x,f10.4,1x,f10.4,
     & 1x,f10.4,1x,f10.4,1x,
     & f10.4,1x,f10.4,1x,f10.4,1x,f10.4,1x,f10.4,1x,f10.4,1x,
     & f10.4,1x,f10.4,1x,f10.4,1x,f10.4,1x,f10.4,1x,f10.4,1x,
     & f10.4,1x,f10.4,1x,f10.4,1x,f10.4,1x,f10.4,1x,f10.4,1x,
     & f10.4,1x,f10.4,1x,f10.4,1x,f10.4,1x,f10.4,1x,f10.4,1x,
     & f10.4,1x,f10.4,1x,f10.4,1x,f10.4,1x,f10.4,1x,f10.4,1x,
     & f10.4,1x,f10.4,1x,f10.4,1x,f10.4,1x,f10.4,1x,f10.4,1x,
     & f10.4,1x,f10.4,1x,f10.4,1x,f10.4,1x,f10.4,1x,f10.4)
11    format (f10.4,1x,f10.4,1x,f10.4,1x,f10.4,1x,f10.4,1x,f10.4,1x,
     & f10.4,1x,f10.4,1x,f10.4,1x,f10.4,1x,f10.4,1x,f10.4,1x,
     & f10.4,1x,f10.4,1x,f10.4,1x,f10.4,1x,f10.4,1x,f10.4,1x,
     & f10.4,1x,f10.4,1x,f10.4,1x,f10.4,1x,f10.4,1x,f10.4,1x,
     & f10.4,1x,f10.4,1x,f10.4,1x,f10.4,1x,f10.4,1x,f10.4,1x,
     & f10.4,1x,f10.4,1x,f10.4,1x,f10.4,1x,f10.4,1x,f10.4,1x,
     & f10.4,1x,f10.4,1x,f10.4,1x,f10.4,1x,f10.4,1x,f10.4,1x,
     & f10.4,1x,f10.4,1x,f10.4,1x,f10.4,1x,f10.4,1x,f10.4)
12    format (i3,2x,i3,2x,f10.1,1x,f10.1,1x,f10.1,1x,f10.1,
     & 1x,f10.1,1x,f10.1,1x,
     & f10.1,1x,f10.1,1x,f10.1,1x,f10.1,1x,f10.1,1x,f10.1,1x,
     & f10.1,1x,f10.1,1x,f10.1,1x,f10.1,1x,f10.1,1x,f10.1,1x,
     & f10.1,1x,f10.1,1x,f10.1,1x,f10.1,1x,f10.1,1x,f10.1,1x,
     & f10.1,1x,f10.1,1x,f10.1,1x,f10.1,1x,f10.1,1x,f10.1,1x,
     & f10.1,1x,f10.1,1x,f10.1,1x,f10.1,1x,f10.1,1x,f10.1,1x,
     & f10.1,1x,f10.1,1x,f10.1,1x,f10.1,1x,f10.1,1x,f10.1,1x,
     & f10.1,1x,f10.1,1x,f10.1,1x,f10.1,1x,f10.1,1x,f10.1)
13    format (f10.1,1x,f10.1,1x,f10.1,1x,f10.1,1x,f10.1,1x,f10.1,1x,
     & f10.1,1x,f10.1,1x,f10.1,1x,f10.1,1x,f10.1,1x,f10.1,1x,
     & f10.1,1x,f10.1,1x,f10.1,1x,f10.1,1x,f10.1,1x,f10.1,1x,
     & f10.1,1x,f10.1,1x,f10.1,1x,f10.1,1x,f10.1,1x,f10.1,1x,
     & f10.1,1x,f10.1,1x,f10.1,1x,f10.1,1x,f10.1,1x,f10.1,1x,
     & f10.1,1x,f10.1,1x,f10.1,1x,f10.1,1x,f10.1,1x,f10.1,1x,
     & f10.1,1x,f10.1,1x,f10.1,1x,f10.1,1x,f10.1,1x,f10.1,1x,
     & f10.1,1x,f10.1,1x,f10.1,1x,f10.1,1x,f10.1,1x,f10.1)

14    format (a20)
15    format (a20,f14.2)
16    format (f12.1)
c------------------------------
c First, initialize all relevant arrays:
c------------------------------
      open (unit=2,file=filename)
      do count1 = 1,numatomtypes1
       do count2 = count1,numatomtypes2
        M(count1,count2) = 0.0
        do count3 = 1,numdistances 
         RFO(count1,count2,count3) = 0.0
         PMFENERGY(count1,count2,count3) = 0.0
        end do
       end do
      end do
      FOALL = 0.0
      do count1 = 1,numdistances
       FXX(count1) = 0.0
       RFOALL(count1) = 0.0
      end do


c calculate Mij - the total number of observations for the atomic pair ij
c   over all distances

      do count1 = 1,numatomtypes1
       do count2 = count1,numatomtypes2
        do count3 = 1,numdistances
         M(count1,count2) = M(count1,count2) + 
     &    pmfShellCounter(count1,count2,count3)
        end do
       end do
      end do

      write (2,14) 'Mij'
      do count1 = 1,numatomtypes1
       write (2,13) (M(count1,count2),count2=count1,numatomtypes2)
      end do
      write (2,14) '  '
      write (2,14) '  '

      write (2,14) 'pmfShellCounter'
      do count1 = 1,numatomtypes1
       do count2 = count1,numatomtypes2
        write (2,12) (count1-1),(count2-1),
     &   (pmfShellCounter(count1,count2,count3),
     &   count3=1,numdistances)
       end do
      end do
      write (2,14) '  '
      write (2,14) '  '

c calculate the relative frequency of occurrence of the atomic pair ij
c in class of distance l

      do count1 = 1,numatomtypes1
       do count2 = count1,numatomtypes2
        do count3 = 1,numdistances
         RFO(count1,count2,count3) = 
     &    pmfShellCounter(count1,count2,count3) /
     &    M(count1,count2)
        end do
       end do
      end do

      write (2,14) 'RFO'
      do count1 = 1,numatomtypes1
       do count2 = count1,numatomtypes2
        write (2,10) (count1-1),(count2-1),
     &    (RFO(count1,count2,count3),count3=1,numdistances)
       end do
      end do

c calculate the relative frequency of occurrence of all atomic pairs in 
c class of distance l

      do count1 = 1,numatomtypes1
       do count2 = count1,numatomtypes2
        do count3 = 1,numdistances
         FOALL = FOALL + pmfShellCounter(count1,count2,count3)
         FXX(count3) = FXX(count3) + 
     &    pmfShellCounter(count1,count2,count3)
        end do
       end do
      end do

      do count3 = 1,numdistances
       RFOALL(count3) = FXX(count3)/FOALL
      end do

      write (2,15) 'FOALL = ',FOALL
      do count3 = 1,numdistances
      write (2,16),FXX(count3)
      end do

c now calculate the energy of a given ij pair in a distance of class l

      do count1 = 1,numatomtypes1
       do count2 = count1,numatomtypes2
        do count3 = 1,numdistances
         PMFENERGY(count1,count2,count3) = 
     &    log( 1+( M(count1,count2)*S) ) -
     &    log( 1+( M(count1,count2)*S*
     &             (RFO(count1,count2,count3)/RFOALL(count3))))
        end do
       end do
      end do

      close (2)
999   end

c-----------------------------------------------------

      subroutine PMFmodule_LowRes (pmfShellCounter,numatomtypes1,
     &   numatomtypes2,numdistances,PMFENERGY,filename)

c-----------------------------------------------------
c Author: C. Summa
c-----------------------------------------------------
c Variable declarations


      include 'parameters.h'

      character*(*) filename
      real PMFENERGY(maxnumatomtypes,maxnumatomtypesEnv,pmfMaxNumBins)
      real pmfShellCounter(maxnumatomtypes,
     &    maxnumatomtypesEnv,pmfMaxNumBins)
      real M(maxnumatomtypes,maxnumatomtypesEnv)
      real RFO(maxnumatomtypes,maxnumatomtypesEnv,pmfMaxNumBins)
      real FOALL
      real FXX(pmfMaxNumBins)
      real RFOALL(pmfMaxNumBins)
      real*8 numatmTotal2Real
      integer  count1,count2,count3,count4,count5,i,j,k
      integer  numatomtypes1
      integer  numatomtypes2
      integer  numdistances
      integer  follow1,follow2,follow3

10    format (i3,2x,i3,2x,f10.4,1x,f10.4,1x,f10.4,1x,f10.4,
     & 1x,f10.4,1x,f10.4,1x,
     & f10.4,1x,f10.4,1x,f10.4,1x,f10.4,1x,f10.4,1x,f10.4,1x,
     & f10.4,1x,f10.4,1x,f10.4,1x,f10.4,1x,f10.4,1x,f10.4,1x,
     & f10.4,1x,f10.4,1x,f10.4,1x,f10.4,1x,f10.4,1x,f10.4,1x,
     & f10.4,1x,f10.4,1x,f10.4,1x,f10.4,1x,f10.4,1x,f10.4,1x,
     & f10.4,1x,f10.4,1x,f10.4,1x,f10.4,1x,f10.4,1x,f10.4,1x,
     & f10.4,1x,f10.4,1x,f10.4,1x,f10.4,1x,f10.4,1x,f10.4,1x,
     & f10.4,1x,f10.4,1x,f10.4,1x,f10.4,1x,f10.4,1x,f10.4)
11    format (f10.4,1x,f10.4,1x,f10.4,1x,f10.4,1x,f10.4,1x,f10.4,1x,
     & f10.4,1x,f10.4,1x,f10.4,1x,f10.4,1x,f10.4,1x,f10.4,1x,
     & f10.4,1x,f10.4,1x,f10.4,1x,f10.4,1x,f10.4,1x,f10.4,1x,
     & f10.4,1x,f10.4,1x,f10.4,1x,f10.4,1x,f10.4,1x,f10.4,1x,
     & f10.4,1x,f10.4,1x,f10.4,1x,f10.4,1x,f10.4,1x,f10.4,1x,
     & f10.4,1x,f10.4,1x,f10.4,1x,f10.4,1x,f10.4,1x,f10.4,1x,
     & f10.4,1x,f10.4,1x,f10.4,1x,f10.4,1x,f10.4,1x,f10.4,1x,
     & f10.4,1x,f10.4,1x,f10.4,1x,f10.4,1x,f10.4,1x,f10.4)
12    format (i3,2x,i3,2x,f10.1,1x,f10.1,1x,f10.1,1x,f10.1,
     & 1x,f10.1,1x,f10.1,1x,
     & f10.1,1x,f10.1,1x,f10.1,1x,f10.1,1x,f10.1,1x,f10.1,1x,
     & f10.1,1x,f10.1,1x,f10.1,1x,f10.1,1x,f10.1,1x,f10.1,1x,
     & f10.1,1x,f10.1,1x,f10.1,1x,f10.1,1x,f10.1,1x,f10.1,1x,
     & f10.1,1x,f10.1,1x,f10.1,1x,f10.1,1x,f10.1,1x,f10.1,1x,
     & f10.1,1x,f10.1,1x,f10.1,1x,f10.1,1x,f10.1,1x,f10.1,1x,
     & f10.1,1x,f10.1,1x,f10.1,1x,f10.1,1x,f10.1,1x,f10.1,1x,
     & f10.1,1x,f10.1,1x,f10.1,1x,f10.1,1x,f10.1,1x,f10.1)
13    format (f10.1,1x,f10.1,1x,f10.1,1x,f10.1,1x,f10.1,1x,f10.1,1x,
     & f10.1,1x,f10.1,1x,f10.1,1x,f10.1,1x,f10.1,1x,f10.1,1x,
     & f10.1,1x,f10.1,1x,f10.1,1x,f10.1,1x,f10.1,1x,f10.1,1x,
     & f10.1,1x,f10.1,1x,f10.1,1x,f10.1,1x,f10.1,1x,f10.1,1x,
     & f10.1,1x,f10.1,1x,f10.1,1x,f10.1,1x,f10.1,1x,f10.1,1x,
     & f10.1,1x,f10.1,1x,f10.1,1x,f10.1,1x,f10.1,1x,f10.1,1x,
     & f10.1,1x,f10.1,1x,f10.1,1x,f10.1,1x,f10.1,1x,f10.1,1x,
     & f10.1,1x,f10.1,1x,f10.1,1x,f10.1,1x,f10.1,1x,f10.1)

14    format (a20)
15    format (a20,f14.2)
16    format (f12.1)
c------------------------------
c First, initialize all relevant arrays:
c------------------------------
      open (unit=2,file=filename)
      do count1 = 1,numatomtypes1
       do count2 = 1,numatomtypes2
        M(count1,count2) = 0.0
        do count3 = 1,numdistances 
         RFO(count1,count2,count3) = 0.0
         PMFENERGY(count1,count2,count3) = 0.0
        end do
       end do
      end do
      FOALL = 0.0
      do count1 = 1,numdistances
       FXX(count1) = 0.0
       RFOALL(count1) = 0.0
      end do


c calculate Mij - the total number of observations for the atomic pair ij
c   over all distances

      do count1 = 1,numatomtypes1
       do count2 = 1,numatomtypes2
        do count3 = 1,numdistances
         M(count1,count2) = M(count1,count2) + 
     &    pmfShellCounter(count1,count2,count3)
        end do
       end do
      end do

      write (2,14) 'Mij'
      do count1 = 1,numatomtypes1
       write (2,13) (M(count1,count2),count2=count1,numatomtypes2)
      end do

      write (2,14) 'pmfShellCounter'
      do count1 = 1,numatomtypes1
       do count2 = 1,numatomtypes2
        write (2,12) (count1-1),(count2-1),
     &   (pmfShellCounter(count1,count2,count3),
     &   count3=1,numdistances)
       end do
      end do

c calculate the relative frequency of occurrence of the atomic pair ij
c in class of distance l

      do count1 = 1,numatomtypes1
       do count2 = 1,numatomtypes2
        do count3 = 1,numdistances
         RFO(count1,count2,count3) = 
     &    pmfShellCounter(count1,count2,count3) /
     &    M(count1,count2)
        end do
       end do
      end do

      write (2,14) 'RFO'
      do count1 = 1,numatomtypes1
       do count2 = 1,numatomtypes2
        write (2,10) (count1-1),(count2-1),
     &    (RFO(count1,count2,count3),count3=1,numdistances)
       end do
      end do

c calculate the relative frequency of occurrence of all atomic pairs in 
c class of distance l

      do count1 = 1,numatomtypes1
       do count2 = 1,numatomtypes2
        do count3 = 1,numdistances
         FOALL = FOALL + pmfShellCounter(count1,count2,count3)
         FXX(count3) = FXX(count3) + 
     &    pmfShellCounter(count1,count2,count3)
        end do
       end do
      end do

      do count3 = 1,numdistances
       RFOALL(count3) = FXX(count3)/FOALL
      end do

      write (2,15) 'FOALL = ',FOALL
      do count3 = 1,numdistances
       write (2,16),FXX(count3)
      end do

c now calculate the energy of a given ij pair in a distance of class l

      do count1 = 1,numatomtypes1
       do count2 = 1,numatomtypes2
        do count3 = 1,numdistances
         PMFENERGY(count1,count2,count3) = 
     &    log( 1+( M(count1,count2)*S) ) -
     &    log( 1+( M(count1,count2)*S*
     &             (RFO(count1,count2,count3)/RFOALL(count3))))
        end do
       end do
      end do

      close (2)
999   end
