      subroutine CalcSig(numatomtypes,comboTypes,
     &  ChiSqrTable,ChiSqrProb,expected,actual,chisq,sig,sigBool,
     &  sigHis,PerSigGT95)

c This version assigns a value of 1 to SigBool under the following conditions:
c 1) The bin has more than 10 counts in it
c 2) The deviation from the Chi-square distribution is greater than 95%

      include 'parameters.h'

      integer numatomtypes
      integer comboTypes(maxNumCoord)
      integer sigBool(maxnumatomtypes,maxNumCoord,570)
      integer*8  sigHis(14)
      real actual(maxnumatomtypes,maxNumCoord,570)
      real sig(maxnumatomtypes,maxNumCoord,570)
      real expected(maxnumatomtypes,maxNumCoord,570)
      real chisq(maxnumatomtypes,maxNumCoord,570)
      real ChiSqrTable(13,30)
      real ChiSqrProb(13)
      real Poisson
      real SumSigHis
      real SigGT95
      real PerSigGT95
      real TotalSigEvents
      real TotalNonSigEvents
      real PercentSigEvents
      logical SigNotAssigned


c-------------------------------------------------------------------
c  Significance Calculation
c-------------------------------------------------------------------
c Calculate chisq statistic for each atomtype-coordination
c as deviation from random for all cases
c Although chisq is not a good measure in cases where
c the number of events is less than ten, we will deal with
c this using the poisson distribution later in the program

       do count1 = 1,numatomtypes
        do count2 = 1,numCoord
         do count3 = 1,comboTypes(count2)
          if (expected(count1,count2,count3).ne.0.000) then
           chisq(count1,count2,count3) =
     +     ((actual(count1,count2,count3) -
     +     expected(count1,count2,count3))**2)/
     +     expected(count1,count2,count3)
          else
           chisq(count1,count2,count3) = 0.00
          end if
         end do
        end do
       end do

c-------------------------------------------------
c Calculate the Significance of each Bin into which
c the data from the pdb coordinates has been placed.

       do i=1,numCoord
        sigHis(i) = 0.0
       end do
       do count1 = 1,numatomtypes
        do count2 = 1,numCoord
         do count3 = 1,comboTypes(count2)

          if (actual(count1,count2,count3).ge.10) then
           sigBool(count1,count2,count3) = 1
           SigNotAssigned = .true.
           do count4 = 1,13
            if (chisq(count1,count2,count3).ge.
     &          ChiSqrTable(count4,1)) then
             if (SigNotAssigned) then
              sig(count1,count2,count3) = 1.0 - ChiSqrProb(count4)
              SigNotAssigned = .false.
             end if
            end if
            if (count4.eq.13) then
             if (SigNotAssigned) then
              sig(count1,count2,count3) = 0.0
             end if
            end if
           end do
          else
           sig(count1,count2,count3) = 1 -
     &     Poisson(actual(count1,count2,count3),
     &     expected(count1,count2,count3))
           if (sig(count1,count2,count3).ge.0.95) then
            sigBool(count1,count2,count3) = 1
           end if
           SigNotAssigned = .true.
           do count4 = 1,13
            if (SigNotAssigned) then
             if (sig(count1,count2,count3).gt.
     &           (1.0-ChiSqrProb(count4))) then
              SigNotAssigned = .false.
             end if
            end if
           end do
          end if

         end do
        end do
       end do

       TotalSigEvents = 0.0
       TotalNonSigEvents = 0.0
       do count1 = 1,numatomtypes
        do count2 = 1,numCoord
         do count3 = 1,comboTypes(count2)
          if (sigBool(count1,count2,count3).eq.1) then
           TotalSigEvents = TotalSigEvents +  
     &      actual(count1,count2,count3)
          end if
          if (sigBool(count1,count2,count3).eq.0) then
           TotalNonSigEvents = TotalNonSigEvents + 
     &      actual(count1,count2,count3)
          end if
         end do
        end do
       end do

       do count1 = 1,numatomtypes
        TotalNonSigEvents = TotalNonSigEvents +
     &  actual(count1,15,1)
       end do
       PercentSigEvents = TotalSigEvents /
     & (TotalSigEvents + TotalNonSigEvents)

       SumSigHis = 0.0
       SigGT95 = 0.0
       PerSigGT95 = 0.0

       do count1 = 1,14
        SumSigHis = SumSigHis + sigHis(count1)
       end do

       do count1 = 1,4
        SigGT95 = SigGT95 + sigHis(count1)
       end do

       PerSigGT95 = SigGT95 / SumSigHis
       end
c-----------------------------------------------------------------
      subroutine CalcSig_pure(numatomtypes,comboTypes,
     &  ChiSqrTable,ChiSqrProb,expected,actual,chisq,sig,sigBool,
     &  sigHis,PerSigGT95)


c A significant bin is one in which the deviation from the mean is greater
c than 95% significant 

      include 'parameters.h'

      integer numatomtypes
      integer comboTypes(maxNumCoord)
      integer sigBool(maxnumatomtypes,maxNumCoord,570)
      integer*8  sigHis(14)
      real actual(maxnumatomtypes,maxNumCoord,570)
      real sig(maxnumatomtypes,maxNumCoord,570)
      real expected(maxnumatomtypes,maxNumCoord,570)
      real chisq(maxnumatomtypes,maxNumCoord,570)
      real ChiSqrTable(13,30)
      real ChiSqrProb(13)
      real Poisson
      real SumSigHis
      real SigGT95
      real PerSigGT95
      real TotalSigEvents
      real TotalNonSigEvents
      real PercentSigEvents
      logical SigNotAssigned


c-------------------------------------------------------------------
c  Significance Calculation
c-------------------------------------------------------------------
c Go through the following calculations, once with the 'normal'
c event data, and then a second time with the 'altered' event data

c Calculate chisq statistic for each atomtype-coordination
c as deviation from random for all cases
c Although chisq is not a good measure in cases where
c the number of events is less than ten, we will deal with
c this using the poisson distribution later in the program

       do count1 = 1,numatomtypes
        do count2 = 1,numCoord
         do count3 = 1,comboTypes(count2)
          if (expected(count1,count2,count3).ne.0.000) then
           chisq(count1,count2,count3) =
     +     ((actual(count1,count2,count3) -
     +     expected(count1,count2,count3))**2)/
     +     expected(count1,count2,count3)
          else
           chisq(count1,count2,count3) = 0.00
          end if
         end do
        end do
       end do


c-------------------------------------------------
c Calculate the Significance of each Bin into which
c the data from the pdb coordinates has been placed.


       do i=1,numCoord
        sigHis(i) = 0.0
       end do

       do count1 = 1,numatomtypes
        do count2 = 1,numCoord
         do count3 = 1,comboTypes(count2)
          sigBool(count1,count2,count3) = 0.0
c Case where # counts is greater than 10

          if (actual(count1,count2,count3).ge.10) then
           SigNotAssigned = .true.
           do count4 = 1,13
            if (chisq(count1,count2,count3).ge.
     &          ChiSqrTable(count4,1)) then
             if (SigNotAssigned) then
              sig(count1,count2,count3) = 1.0 - ChiSqrProb(count4)
              sigHis(count4) = sigHis(count4)+1
              SigNotAssigned = .false.
             end if
            end if
           end do
           if (SigNotAssigned) then
             sig(count1,count2,count3) = 0.0
             sigHis(14) = sigHis(14) + 1
           end if

           if (.not.SigNotAssigned.and.
     &         sig(count1,count2,count3).ge.0.95) then
             sigBool(count1,count2,count3) = 1
           end if

          else

c          Case where # counts is less than 10
           SigNotAssigned = .true.
           sig(count1,count2,count3) = 1 -
     &     Poisson(actual(count1,count2,count3),
     &     expected(count1,count2,count3))

           if (sig(count1,count2,count3).ge.0.95) then
            sigBool(count1,count2,count3) = 1
           end if

           do count4 = 1,13
            if (SigNotAssigned) then
             if (sig(count1,count2,count3).gt.
     &           (1.0-ChiSqrProb(count4))) then
              sigHis(count4) = sigHis(count4) + 1
              SigNotAssigned = .false.
             end if
            end if
           end do

          end if
         end do
        end do
       end do

       TotalSigEvents = 0.0
       TotalNonSigEvents = 0.0
       do count1 = 1,numatomtypes
        do count2 = 1,numCoord
         do count3 = 1,comboTypes(count2)
          if (sigBool(count1,count2,count3).eq.1) then
           TotalSigEvents = TotalSigEvents +  
     &      actual(count1,count2,count3)
          end if
          if (sigBool(count1,count2,count3).eq.0) then
           TotalNonSigEvents = TotalNonSigEvents + 
     &      actual(count1,count2,count3)
          end if
         end do
        end do
       end do

       do count1 = 1,numatomtypes
        TotalNonSigEvents = TotalNonSigEvents +
     &  actual(count1,15,1)
       end do

       PercentSigEvents = TotalSigEvents /
     & (TotalSigEvents + TotalNonSigEvents)

       SumSigHis = 0.0
       SigGT95 = 0.0
       PerSigGT95 = 0.0

       do count1 = 1,14
        SumSigHis = SumSigHis + sigHis(count1)
       end do

       do count1 = 1,4
        SigGT95 = SigGT95 + sigHis(count1)
       end do

       PerSigGT95 = SigGT95 / SumSigHis
       end
c---------------------------------------	
      function Poisson(x1,mu)
      
      real x1
      real mu
      real Poisson
      double precision Fact

      Poisson = (exp(-mu)*(mu**x1))/Fact(int(x1))
      end
c--------------------------------------
       subroutine ReadChiSqrTable(contactdir,ChiSqrTable,ChiSqrProb)

       character*80 chisqrname
       character*(*) contactdir
       character*120 concat1
      
       real ChiSqrTable(13,30)
       real ChiSqrProb(13)
       logical exist

10     format(4x,f6.2,1x,f6.2,1x,f6.2,1x,f6.2,1x,f6.2,
     &   1x,f6.2,1x,f6.2,1x,f6.2,1x,f6.2,1x,f6.2,1x,f6.2,
     &   1x,f6.2,1x,f6.2)
11     format(4x,f6.3,1x,f6.3,1x,f6.3,1x,f6.3,1x,f6.3,
     &   1x,f6.3,1x,f6.3,1x,f6.3,1x,f6.3,1x,f6.3,1x,f6.3,
     &   1x,f6.3,1x,f6.3)

      chisqrname = "ChiSquared.dat"
      call fileexist(chisqrname,contactdir,concat1,exist)
      if (exist) then
       open (unit=2,file=concat1)
      else
       print *,'Error opening: ',concat1
       print *,'Stopping'
       stop
      end if

      read (2,10) (ChiSqrProb(i),i=1,13)
      do j=1,30
       read (2,11) (ChiSqrTable(i,j),i=1,13)  
      end do
      close (2)

      end
