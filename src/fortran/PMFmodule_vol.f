      subroutine PMFmodule_vol_raw (pmfShellCounter,numatomtypes1,
     &   numatomtypes2,numdistances,PMFENERGY,mindis,maxdis,
     &   stepsize)

c-----------------------------------------------------
c Author: C. Summa
c Modifications: 
c	Version 1.0 - 

c Last modified: 6/30/98
c-------------------------------------------------------
c Variable declarations

      include 'parameters.h'

      real vol(maxnumdistances)
      real mindis,maxdis,stepsize
      real templow,temphigh
      real tempreal
      real PMFENERGY(maxnumatomtypes,maxnumatomtypes,pmfMaxNumBins)
      real pmfShellCounter(maxnumatomtypes,
     &    maxnumatomtypes,pmfMaxNumBins)
      real expectedFreqDist(pmfMaxNumBins)
      real numdistances_real
      real relNumEvents(maxnumatomtypes,maxnumatomtypes)
      real totalNumEvents
      
      integer  count1,count2,count3,count4,count5,i,j,k
      integer  numatomtypes1
      integer  numatomtypes2
      integer  numdistances

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

14	  format (a20)
15	  format (a20,f14.2)
16    format (f12.1)
c------------------------------
c First, initialize all relevant arrays:
c------------------------------
	  open (unit=2,file='pmf.rawdata')

      do count1 = 1,numatomtypes1
       do count2 = count1,numatomtypes2
        relNumEvents(count1,count2) = 0.0
        do count3 = 1,numdistances 
         expectedFreqDist(count3) = 0.0
         PMFENERGY(count1,count2,count3) = 0.0
        end do
       end do
      end do
      totalNumEvents = 0.0

      write (2,14) 'pmfShellCounter'
      do count1 = 1,numatomtypes1
       do count2 = count1,numatomtypes2
        write (2,12) (count1-1),(count2-1),
     &   (pmfShellCounter(count1,count2,count3),
     &   count3=1,numdistances)
       end do
      end do

      templow = mindis
      temphigh = mindis + stepsize
      do count1 = 1,numdistances
        vol(count1) = (1.3333 * 3.14159 * (temphigh**3) ) -
     &                (1.3333 * 3.14159 * (templow**3) )
        templow = temphigh
        temphigh =  temphigh + stepsize
      end do

      write (2,14) 'volume'
      write (2,13) (vol(count1),count1 = 1,numdistances)

      do count1 = 1,numatomtypes1
       do count2 = count1,numatomtypes2
        do count3 = 1,numdistances
         expectedFreqDist(count3) =
     &    expectedFreqDist(count3) +
     &    pmfShellCounter(count1,count2,count3)
        end do
       end do
      end do 

      do count1 = 1,numatomtypes1
       do count2 = count1,numatomtypes2
        do count3 = 1,numdistances
         relNumEvents(count1,count2) =
     &    relNumEvents(count1,count2) +
     &    pmfShellCounter(count1,count2,count3)
         totalNumEvents = totalNumEvents + 
     &    pmfShellCounter(count1,count2,count3)
        end do
       end do
      end do 

c now calculate the energy of a given ij pair in a distance of class l

      do count1 = 1,numatomtypes1
       do count2 = count1,numatomtypes2
        do count3 = 1,numdistances
         if (expectedFreqDist(count3).eq.0.0.or.
     &       relNumEvents(count1,count2).eq.0.0.or.
     &       pmfShellCounter(count1,count2,count3).eq.0.0) then
           do count4 = count3,1,-1
             PMFENERGY(count1,count2,count4) = 20.0
           end do
         else
           PMFENERGY(count1,count2,count3) = 
     &     -log( pmfShellCounter(count1,count2,count3)/
     &          (expectedFreqDist(count3)*
     &           relNumEvents(count1,count2)/totalNumEvents) )
         end if
        end do
       end do
      end do

      close (2)
999   end

c-----------------------------------------------------------------------

      subroutine PMFmodule_vol_norm (pmfShellCounter,numatomtypes1,
     &   numatomtypes2,numdistances,PMFENERGY,mindis,maxdis,
     &   stepsize)

c-----------------------------------------------------
c Author: C. Summa
c Modifications: 
c	Version 1.0 - 

c Last modified: 6/30/98
c-------------------------------------------------------
c Variable declarations

      include 'parameters.h'

      real vol(maxnumdistances)
      real mindis,maxdis,stepsize
      real templow,temphigh
      real tempreal
      real PMFENERGY(maxnumatomtypes,maxnumatomtypes,pmfMaxNumBins)
      real pmfShellCounter(maxnumatomtypes,
     &    maxnumatomtypes,pmfMaxNumBins)
      real rFreq(maxnumatomtypes,maxnumatomtypes,pmfMaxNumBins)
      real avgFreq(maxnumatomtypes,maxnumatomtypes)
      real numdistances_real
      real bigvol
      integer  count1,count2,count3,count4,count5,i,j,k
      integer  numatomtypes1
      integer  numatomtypes2
      integer  numdistances

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

14	  format (a20)
15	  format (a20,f14.2)
16    format (f12.1)
c------------------------------
c First, initialize all relevant arrays:
c------------------------------
	  open (unit=2,file='pmf.rawdata')
      do count1 = 1,numatomtypes1
       do count2 = count1,numatomtypes2
        avgFreq(count1,count2) = 0.0
        do count3 = 1,numdistances 
         rFreq(count1,count2,count3) = 0.0
         PMFENERGY(count1,count2,count3) = 0.0
        end do
       end do
      end do

      write (2,14) 'pmfShellCounter'
      do count1 = 1,numatomtypes1
       do count2 = count1,numatomtypes2
        write (2,12) (count1-1),(count2-1),
     &   (pmfShellCounter(count1,count2,count3),
     &   count3=1,numdistances)
       end do
      end do

c calculate the relative density of atomic pairs ij
c in class of distance l

      templow = mindis
      temphigh = mindis + stepsize
      do count1 = 1,numdistances
        vol(count1) = (1.3333 * 3.14159 * (temphigh**3) ) -
     &                (1.3333 * 3.14159 * (templow**3) )
        templow = temphigh
        temphigh =  temphigh + stepsize
      end do

      write (2,14) 'volume'
      write (2,13) (vol(count1),count1 = 1,numdistances)

      do count1 = 1,numatomtypes1
       do count2 = count1,numatomtypes2
        do count3 = 1,numdistances
         rFreq(count1,count2,count3) =
     &    pmfShellCounter(count1,count2,count3) /
     &    vol(count3)
        end do
       end do
      end do 

      write (2,14) 'rFreq'
      do count1 = 1,numatomtypes1
       do count2 = count1,numatomtypes2
        write (2,10) (count1-1),(count2-1),
     &    (rFreq(count1,count2,count3),count3=1,numdistances)
       end do
      end do

c calculate average density for an ij pair

      do count1 = 1,numatomtypes1
       do count2 = count1,numatomtypes2
        do count3 = 1,numdistances
         avgFreq(count1,count2) = avgFreq(count1,count2) + 
     &    pmfShellCounter(count1,count2,count3)
        end do
       end do
      end do

      bigvol = 1.3333 * 3.14159 * maxdis**3 
      print *,'maxdis = ',maxdis
      print *,'bigvol = ',bigvol
      do count1 = 1,numatomtypes1
       do count2 = count1,numatomtypes2
        avgFreq(count1,count2) = avgFreq(count1,count2)/bigvol
       end do
      end do

      write (2,14) 'AvgFreq'
      do count1 = 1,numatomtypes1
       write (2,13) (avgFreq(count1,count2),count2=count1,numatomtypes2)
      end do

c now calculate the energy of a given ij pair in a distance of class l

      do count1 = 1,numatomtypes1
       do count2 = count1,numatomtypes2
        do count3 = 1,numdistances
         if (rFreq(count1,count2,count3).eq.0.0) then
           do count4 = count3,1,-1
             PMFENERGY(count1,count2,count4) = 10.0
           end do
         else
           PMFENERGY(count1,count2,count3) = 
     &     -log(rFreq(count1,count2,count3)/avgFreq(count1,count2))
         end if
        end do
       end do
      end do

c normalize so that the final value is 0.000
      do count1 = 1,numatomtypes1
       do count2 = count1,numatomtypes2
        tempreal = PMFENERGY(count1,count2,numdistances)
        do count3 = 1,numdistances
          if (PMFENERGY(count1,count2,count3).ne.10.0) then
            PMFENERGY(count1,count2,count3) =
     &      PMFENERGY(count1,count2,count3) - tempreal
          end if 
        end do
       end do
      end do

      close (2)
999   end

c-----------------------------------------------------------------------

      subroutine PMFmodule_vol_norm_lores (pmfShellCounter,
     &   numatomtypes1,
     &   numatomtypes2,numdistances,PMFENERGY,mindis,maxdis,
     &   stepsize)

c-----------------------------------------------------
c Author: C. Summa
c Modifications: 
c	Version 1.0 - 

c Last modified: 6/30/98
c-------------------------------------------------------
c Variable declarations

      include 'parameters.h'

      real vol(maxnumdistances)
      real mindis,maxdis,stepsize
      real templow,temphigh
      real tempreal
      real PMFENERGY(maxnumatomtypes,maxnumatomtypesEnv,pmfMaxNumBins)
      real pmfShellCounter(maxnumatomtypes,
     &    maxnumatomtypesEnv,pmfMaxNumBins)
      real rFreq(maxnumatomtypes,maxnumatomtypesEnv,pmfMaxNumBins)
      real avgFreq(maxnumatomtypes,maxnumatomtypesEnv)
      real numdistances_real
      real bigvol
      integer  count1,count2,count3,count4,count5,i,j,k
      integer  numatomtypes1
      integer  numatomtypes2
      integer  numdistances

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

14	  format (a20)
15	  format (a20,f14.2)
16    format (f12.1)
c------------------------------
c First, initialize all relevant arrays:
c------------------------------
	  open (unit=2,file='pmf.rawdata')
      do count1 = 1,numatomtypes1
       do count2 = 1,numatomtypes2
        avgFreq(count1,count2) = 0.0
        do count3 = 1,numdistances 
         rFreq(count1,count2,count3) = 0.0
         PMFENERGY(count1,count2,count3) = 0.0
        end do
       end do
      end do

      write (2,14) 'pmfShellCounter'
      do count1 = 1,numatomtypes1
       do count2 = 1,numatomtypes2
        write (2,12) (count1-1),(count2-1),
     &   (pmfShellCounter(count1,count2,count3),
     &   count3=1,numdistances)
       end do
      end do

c calculate the relative density of atomic pairs ij
c in class of distance l

      templow = mindis
      temphigh = mindis + stepsize
      do count1 = 1,numdistances
        vol(count1) = (1.3333 * 3.14159 * (temphigh**3) ) -
     &                (1.3333 * 3.14159 * (templow**3) )
        templow = temphigh
        temphigh =  temphigh + stepsize
      end do

      write (2,14) 'volume'
      write (2,13) (vol(count1),count1 = 1,numdistances)

      do count1 = 1,numatomtypes1
       do count2 = 1,numatomtypes2
        do count3 = 1,numdistances
         rFreq(count1,count2,count3) =
     &    pmfShellCounter(count1,count2,count3) /
     &    vol(count3)
        end do
       end do
      end do 

      write (2,14) 'rFreq'
      do count1 = 1,numatomtypes1
       do count2 = 1,numatomtypes2
        write (2,10) (count1-1),(count2-1),
     &    (rFreq(count1,count2,count3),count3=1,numdistances)
       end do
      end do

c calculate average density for an ij pair

      do count1 = 1,numatomtypes1
       do count2 = 1,numatomtypes2
        do count3 = 1,numdistances
         avgFreq(count1,count2) = avgFreq(count1,count2) + 
     &    pmfShellCounter(count1,count2,count3)
        end do
       end do
      end do

      bigvol = 1.3333 * 3.14159 * maxdis**3 
      print *,'maxdis = ',maxdis
      print *,'bigvol = ',bigvol
      do count1 = 1,numatomtypes1
       do count2 = 1,numatomtypes2
        avgFreq(count1,count2) = avgFreq(count1,count2)/bigvol
       end do
      end do

      write (2,14) 'AvgFreq'
      do count1 = 1,numatomtypes1
       write (2,13) (avgFreq(count1,count2),count2=1,numatomtypes2)
      end do

c now calculate the energy of a given ij pair in a distance of class l

      do count1 = 1,numatomtypes1
       do count2 = 1,numatomtypes2
        do count3 = 1,numdistances
         if (rFreq(count1,count2,count3).eq.0.0) then
           do count4 = count3,1,-1
             PMFENERGY(count1,count2,count4) = 10.0
           end do
         else
           PMFENERGY(count1,count2,count3) = 
     &     log( rFreq(count1,count2,count3)/avgFreq(count1,count2) )
         end if
        end do
       end do
      end do

c normalize so that the final value is 0.000
      do count1 = 1,numatomtypes1
       do count2 = 1,numatomtypes2
        tempreal = PMFENERGY(count1,count2,numdistances)
        do count3 = 1,numdistances
          if (PMFENERGY(count1,count2,count3).ne.10.0) then
            PMFENERGY(count1,count2,count3) =
     &      PMFENERGY(count1,count2,count3) - tempreal
          end if 
        end do
       end do
      end do

      close (2)
999   end
