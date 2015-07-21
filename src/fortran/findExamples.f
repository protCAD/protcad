      program deriveMicroEnv

c Author: C. Summa
c 3/1/98- Passes significance data into bubblesort algortithm. Only those
c 	  clusters which are significant get sorted into TwentyBest/TwentyWorst
c 6/19/99-
c	Tests for existence of database and testset and acts accordingly
c 
c 8/29/00 - cleaned up a lot of junk variables
c 9/21/00 - if a microEnvironment is not seen, instead of assigning
c 	to it the singleBodyEnergy, assign a value of 5.0 (Bad)
c   (this is on line ~1220)
c   - Removed all references to RT
c-------------------------------------------------------
c Variable declarations

	

      include 'parameters.h'
      include 'pdb.h'

c Related to pdbselect database
      character*11 pdbnames(maxNumPDBs)
      character*1  pdbChain(maxNumPDBs)
      integer      thrsh(maxNumPDBs)
      integer      nHlx(maxNumPDBs)
      integer      nBta(maxNumPDBs)
      integer      numpdbnames

      real Distance
      real tempdistance
      real criticalmaxdistance
      real criticalmindistance 
      real fraction(maxnumatomtypes)
      real probCatCoord(maxnumatomtypes,maxNumCoord,
     &       maxnumatomtypesEnv)
      real PA,PB,PC,PD
      real TheoreticalProbWordGCoord(maxnumatomtypes,maxNumCoord,570)
      real TheoreticalProbWord(maxnumatomtypes,maxNumCoord,570)
      real TheoreticalProbCoord(maxnumatomtypes,maxNumCoord)
      real ActualProbWordGCoord(maxnumatomtypes,maxNumCoord,570)
      real ActualProbWord(maxnumatomtypes,maxNumCoord,570)
      real ActualProbCoord(maxnumatomtypes,maxNumCoord)
      real Eword(maxnumatomtypes,maxNumCoord,570)
      real sumEvents(maxnumatomtypes,maxNumCoord)
      real totalSumAtomsCatCoord
      real ExpectedEvents(maxnumatomtypes,maxNumCoord,570)
      real ChiSquared(maxnumatomtypes,maxNumCoord,570)
      real ChiSqrTable(13,30)
      real ChiSqrProb(13)
      real Significance(maxnumatomtypes,maxNumCoord,570)
      real PerSigGT95 
      real SumUnusableEvents
      real SumTotalEvents
      real NumberSignificantlyPopulated
      real NumberNotSignificantlyPopulated
      real sumAtomsCatCoord(maxnumatomtypes,maxNumCoord,
     &       maxnumatomtypesEnv)
      real event(maxnumatomtypes,maxNumCoord,570)
      real alteredEvent(maxnumatomtypes,maxNumCoord,570)
      real TotalSignificantEvents
      real TotalNonSignificantEvents
      real PercentageSignificantEvents
      real rms
      real singleBodyEnergy(maxnumatomtypes,binsInHistogram)
      real normConst
      real CoordHist(maxnumatomtypes,binsInHistogram)
      real CoordHistAverage(binsInHistogram)
      real CoordHistAverageNorm(binsInHistogram)
      integer numBinsInHistogram
      real binsInHistogramReal
      real overallSumEvents
      real version
      real TotalTheoreticalProbWord
      real numatomsInCategory(maxnumatomtypes,3)
      real numatomsInCategoryTotal(maxnumatomtypes,3)
      real tempreal
      real TwentyBestE(maxnumatomtypes,20)
      real TwentyWorstE(maxnumatomtypes,20) 
      real numenvTotal
      real pmfShellCounter(maxnumatomtypes,maxnumatomtypesEnv,
     &       pmfMaxNumBins)
      real pmfShellCounterHiRes(maxnumatomtypes,maxnumatomtypes,
     &       pmfMaxNumBins)
      real pmfMaxdis
      real pmfMindis
      real PMFENERGY(maxnumatomtypes,maxnumatomtypesEnv,pmfMaxNumBins)
      real PMFHiRes(maxnumatomtypes,maxnumatomtypes,pmfMaxNumBins)
      real pmfAvg(maxnumatomtypes,maxnumatomtypesEnv)
      real pmfAvgMaxshellReal
      real*8 numatmTotal2Real
      double precision Fact

      integer iseed

      real ThingsinEnv(maxnumatomtypes,maxNumCoord,maxnumatomtypesEnv)
      real SumThingsinEnv(maxnumatomtypes,maxNumCoord)
      real numInCatCoord(maxnumatomtypes,maxNumCoord)

      real    pair(maxnumatomtypes,maxnumatomtypesEnv)
      real    pairAvg(maxnumatomtypesEnv)
      real    N1(maxnumatomtypes)
      real    p1(maxnumatomtypes,maxnumatomtypesEnv)
      real    EPair(maxnumatomtypes,maxnumatomtypesEnv)

      real A1
      real A2
      real A3

      integer  null
      integer  category(maxat,3)
      integer  i1,i2,i3,i4,i,j,k
      integer  pdbsread
      integer  numarguments
      integer  numskip
      integer  Coordination(maxat)
      integer  tempInteger
      integer  numatomtypes
      integer  numatomtypesEnv
      integer  numatomtypesTop
      integer  comboTypes(maxNumCoord)
      integer  comboTypesComp(maxNumCoord,570,maxnumatomtypesEnv)
      integer  atomData(maxat,7)
      integer  criticalCoordination
      integer  numA,numB,numC,numD
      integer  SignificanceBool(maxnumatomtypes,maxNumCoord,570)
      integer  mapAtomTopology(maxnumatomtypes) 
      integer  follow1,follow2,follow3
      integer  jason
      integer  TwentyBestEMap(maxnumatomtypes,20,2)
      integer  TwentyWorstEMap(maxnumatomtypes,20,2)
      integer  numchainused
      integer  nHlxTot
      integer  nBtaTot
      integer  pmfShell
      integer  maxpmfShell
      integer  numdistances
      integer  pmfAvgMaxshell
      integer*8  numatmTotal2
      integer*8  Pmult
      integer*8  SignificanceHistogram(14)
      character*80  numskipstr,crdiststr,iseedstr
      character*50  selectname
      character*44  databasedir
      character*50  pdbname
      character*45  contactpdbdir
      character*80  misfolddir
      character*80  contactdir
      character*30  energyfile
      character*30  categorylistname
      character*30  pdblistname2
      character*20  reason
      character*11  filename
      character*1   perResidueString
      logical  relevant
      logical  Abool,Bbool,Cbool,Dbool
      logical  firstTimeThrough
      logical  firsttime
      logical  badAtom
      logical  perResidue
      logical  rejectflag(maxNumPDBs)

      character*160 dbconcat
      character*160 mfconcat
      character*160 mflistconcat
      character*160 enconcat
      logical dbexist
      logical mfexist
      logical enexist

      real BestEnergy
      real WorstEnergy
      real numbestnearhetero
      real numworstnearhetero
      real atomEnergy(maxat)
      integer sortedMap(maxat)
      integer bestarray(10)
      integer worstarray(10)
      integer numtoflag
      character*80 pmfFile

11    format (i8,i8,i8,i8,i8,f8.3,f8.3,f8.3,f8.3)
12    format (a50)
14    format (i4,1x,f8.1,f8.1,f8.1,f8.1,f8.1,
     + f8.1,f8.1,f8.1,f8.1,f8.1,f8.1,
     + f8.1,f8.1,f8.1,f8.1,f8.1,f8.1,
     + f8.1,f8.1,f8.1,f8.1,f8.1,f8.1,
     + f8.1,f8.1,f8.1,f8.1,f8.1,f8.1)
15    format (i4,1x,f8.4,f8.4,f8.4,f8.4,f8.4,f8.4,
     + f8.4,f8.4,f8.4,f8.4,f8.4,f8.4,f8.4,f8.4,f8.4,
     + f8.4,f8.4,f8.4,f8.4,f8.4,f8.4,f8.4,f8.4,f8.4,
     + f8.4,f8.4,f8.4,f8.4,f8.4,f8.4,f8.4,f8.4,f8.4)
18    format (i4,1x,f8.0,f8.0,f8.0,f8.0,f8.0,f8.0,
     + f8.0,f8.0,f8.0,f8.0,f8.0,f8.0,f8.0,f8.0,
     + f8.0,f8.0,f8.0,f8.0,f8.0,f8.0,f8.0,f8.0,
     + f8.0,f8.0,f8.0,f8.0,f8.0,f8.0,f8.0,
     + f8.0,f8.0,f8.0,f8.0,f8.0,f8.0,f8.0)
19    format (i2,2x,f7.3,f7.3,f7.3,f7.3)
28    format (5x,f8.1,f8.1,f8.1,f8.1,f8.1,f8.1,
     + f8.1,f8.1,f8.1,f8.1,f8.1,f8.1,f8.1,f8.1,
     + f8.1,f8.1,f8.1,f8.1,f8.1,f8.1,f8.1,f8.1,
     + f8.1,f8.1,f8.1,f8.1,f8.1,f8.1,f8.1,
     + f8.1,f8.1,f8.1,f8.1,f8.1,f8.1,f8.1)
29    format (5x,i8,i8,i8,i8,i8,i8,
     + i8,i8,i8,i8,i8,i8,i8,i8,
     + i8,i8,i8,i8,i8,i8,i8,i8,
     + i8,i8,i8,i8,i8,3x,1a,i2)
33    format (i2,1x,i2,1x,i3,1x,f9.6)
34    format (i4,1x,f8.3,f8.3,f8.3,f8.3,f8.3,
     + f8.3,f8.3,f8.3,f8.3,f8.3,f8.3,
     + f8.3,f8.3,f8.3,f8.3,f8.3,f8.3,
     + f8.3,f8.3,f8.3,f8.3,f8.3,f8.3,
     + f8.3,f8.3,f8.3,f8.3,f8.3,f8.3)
35    format (i4,1x,f10.1,f10.1,f10.1,f10.1,f10.1,
     + f10.1,f10.1,f10.1,f10.1,f10.1,f10.1,
     + f10.1,f10.1,f10.1,f10.1,f10.1,f10.1,
     + f10.1,f10.1,f10.1,f10.1,f10.1,f10.1,
     + f10.1,f10.1,f10.1,f10.1,f10.1,f10.1,
     + f10.1,f10.1,f10.1,f10.1,f10.1,f10.1,
     + f10.1,f10.1,f10.1,f10.1,f10.1,f10.1,
     + f10.1,f10.1,f10.1,f10.1,f10.1,f10.1,
     + f10.1,f10.1,f10.1,f10.1,f10.1,f10.1,
     + f10.1,f10.1,f10.1,f10.1,f10.1,f10.1,
     + f10.1,f10.1,f10.1,f10.1,f10.1,f10.1)
36    format (8x,a4,i8,i8,i8,i8,i8,i8,
     + i8,i8,i8,i8,i8,i8,i8,i8,
     + i8,i8,i8,i8,i8,i8,i8,i8,
     + i8,i8,i8,i8,i8)
38    format (5x,f8.1,f8.1,f8.1,f8.1,f8.1,f8.1,
     + f8.1,f8.1,f8.1,f8.1,f8.1,f8.1,f8.1,f8.1,
     + f8.1,f8.1,f8.1,f8.1,f8.1,f8.1,f8.1,f8.1,
     + f8.1,f8.1,f8.1,f8.1,f8.1,f8.1,f8.1,
     + f8.1,f8.1,f8.1,f8.1,f8.1,f8.1,f8.1)
40    format (a50,2x,i7,2x,i4,2x,f10.3,1x,f10.3,1x,f10.3,
     & 1x,f10.3,1x,f10.3,1x,f8.0,1x,f8.0)
41    format (i4,1x,i8,i8,i8,i8,i8,i8,
     + i8,i8,i8,i8,i8,i8,i8,i8,
     + i8,i8,i8,i8,i8,i8,i8,i8,
     + i8,i8,i8,i8,i8,i8,i8,
     + i8,i8,i8,i8,i8,i8,i8)
42    format (i4,1x,i8,i8,i8,i8,i8,i8,f8.3)
43    format (i4,1x,f10.1)
44    format (f8.3,1x,i4,2x,i4,1x,i4,1x,i4,1x,i4,1x,i4,
     & 1x,f8.0,1x,f8.0)
66    format (i4,i4,i4,i4,i4,i4,i4,i4,i4)
67    format (a15,a15,3x,f8.5)
68    format (f8.3,1x,f8.3,1x,f8.3,1x,f8.3,1x,f8.3,1x,f8.3,1x,
     + f8.3,1x,f8.3,1x,f8.3,1x,f8.3,1x,f8.3,1x,f8.3,1x,
     + f8.3,1x,f8.3,1x,f8.3,1x,f8.3,1x,f8.3,1x,f8.3,1x,
     + f8.3,1x,f8.3,1x,f8.3,1x,f8.3,1x,f8.3,1x,f8.3,1x,
     + f8.3,1x,f8.3,1x,f8.3,1x,f8.3,1x,f8.3,1x,f8.3,1x,
     + f8.3,1x,f8.3,1x,f8.3,1x,f8.3,1x,f8.3,1x,f8.3,1x,
     + f8.3,1x,f8.3,1x,f8.3,1x,f8.3,1x,f8.3,1x,f8.3,1x,
     + f8.3,1x,f8.3,1x,f8.3,1x,f8.3,1x,f8.3,1x,f8.3,1x,
     + f8.3,1x,f8.3,1x,f8.3,1x,f8.3,1x,f8.3,1x,f8.3,1x,
     + f8.3,1x,f8.3,1x,f8.3,1x,f8.3,1x,f8.3,1x,f8.3,1x,
     + f8.3,1x,f8.3,1x,f8.3,1x,f8.3,1x,f8.3,1x,f8.3)
69    format (f6.0,f6.0,f6.0,f6.0,f6.0,f6.0,f6.0,f6.0,f6.0,f6.0,
     + f6.0,f6.0,f6.0,f6.0,f6.0,f6.0,f6.0,f6.0,f6.0,f6.0,
     + f6.0,f6.0,f6.0,f6.0,f6.0,f6.0,f6.0,f6.0,f6.0,f6.0,
     + f6.0,f6.0,f6.0,f6.0,f6.0,f6.0,f6.0,f6.0,f6.0,f6.0,
     + f6.0,f6.0,f6.0,f6.0,f6.0,f6.0,f6.0,f6.0,f6.0,f6.0)
70    format (i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,
     + i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,
     + i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,
     + i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,
     + i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,
     + i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,
     + i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,
     + i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,
     + i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8)
71    format (a11,2x,a6,1x,f7.3,1x,i2,1x,i2,1x,i2,
     &  1x,a3,1x,i5,1x,a1,2x,a4)
72    format (i4,1x,a3,1x,a1,1x,i2,1x,a3,1x,i2,1x,i2,1x,i2,1x,
     + i2,1x,i2,f8.3,1x,i1)
73    format (a150)
74    format (a9,42x,a10,2x,a10,a10,2x,a10,1x,a10,1x,a10)
75    format (a50,2x,i10,1x,f10.3,1x,f10.3,1x,f10.3,
     & 1x,f10.3,1x,f10.3,1x,f10.3)
76    format (a50,2x,f10.3,1x,f10.3,1x,f10.3,1x,f10.3,1x,f10.3,1x)
77    format (i3,1x,i3,1x,f10.7,1x,f10.7,1x,f10.7,1x,f10.7)
78    format (a2,1x,i2,1x,i2,1x,i3,1x,a1,1x,i2,1x,i2,1x,i2,1x,i2,
     & 1x,a1,1x,f8.1,3x,a1,f8.1,2x,f6.2)
79    format (a20,1x,a5,1x,a5,i12,i12)

c------------------------------------------------------------
c Read in command line arguments

      numarguments = iargc()

      if (numarguments.ne.3) then
        print *,'Usage: deriveMicroEnv PDBselectname #residues_to_skip r
     &adius'
        goto 999
      endif

      call getarg(1,selectname)
      call getarg(2,numskipstr)
      call gettheinteger(numskip,numskipstr)
      call getarg(3,crdiststr)
      call getthereal(criticalmaxdistance,crdiststr)

      call getenv('CONTACTDIR',contactdir)
      call getenv('CONTACTPDBDIR',contactpdbdir)
      call getenv('CONTACTPDBDIR2',databasedir)
      call getenv('MISFOLDDIR',misfolddir)

c      criticalmaxdistance = criticalmaxdistance / 10.0
      version = 11.1
 
      perResidueString = 'n'

      energyfile = ''
      energyfile(1:3) = crdiststr(1:3)
      energyfile(4:5) = '_'
      energyfile(5:5) = numskipstr(1:1)
      energyfile(6:6) = '_'
      energyfile(7:14) = '11.1_std'

      print *,'energyfile = ',energyfile

      follow1 = 5
      follow2 = 3
      follow3 = 3
      jason = 6

      if (perResidueString.eq.'y'.or.perResidueString.eq.'Y') then
       perResidue = .true.
      else
       perResidue = .false.
      end if


     
c------------------------------------------------
c print out some of the input data and values
c of relevant variables

      print *,'Contact version: ', version
      print *,'PDBSelect list used to generate potential: ',selectname
      print *,'Number of residues to skip: ',numskip
      print *,'Counting Radius: ',criticalmaxdistance
      print *

c     Test for existence of database for derivation:
c          databasedir

      call fileexist(selectname,databasedir,dbconcat,dbexist)
      if (.not.dbexist) then
       print *
       print *,'No database...running in test only mode'
       print *
      end if

c------------------------------------------------------------
c Initialize all the variables before starting

      pmfMaxdis = 2*criticalmaxdistance
      pmfMindis = 1.8
      criticalmindistance = pmfMindis
      maxpmfShell = anint((pmfMaxdis-pmfMindis)/pmfStep)
      numdistances = maxpmfShell
      pmfAvgMaxshell = anint((criticalmaxdistance-pmfMindis)/pmfStep)
      pmfAvgMaxshellReal = float(pmfAvgMaxshell)
      numatomtypes = 20
      numatomtypesEnv = 4
      

      A1 = 1.00
      A2 = 1.00
      A3 = 1.00


c      print *,'pmfMaxdis = ',pmfMaxdis
c      print *,'pmfMindis = ',pmfMindis
c      print *,'maxpmfShell = ',maxpmfShell
c      print *,'numdistances = ',numdistances
c      print *,'pmfAvgMaxshell = ',pmfAvgMaxshell
c      print *,'pmfAvgMaxshellReal = ',pmfAvgMaxshellReal

      categorylistname = 'FullAtomList9.dat'
     
      numBinsinHistogram = binsInHistogram
      binsInHistogramReal = float(numBinsInHistogram)
      firstTimeThrough = .true.
      null = 0
      overallSumEvents = 0.0
      SumTotalEvents = 0.0
      SumUnusableEvents = 0.0
      do i=1,maxNumPDBs
       nHlx(i) = 0
       nBta(i) = 0
      end do
      nHlxTot = 0
      nBtaTot = 0
      do i=1,maxnumatomtypes
       do j=1,15
        do k=1,570
         event(i,j,k) = 0.0
         alteredEvent(i,j,k) = 0.0
         SignificanceBool(i,j,k) = 0
        end do
        do k=1,4
         sumAtomsCatCoord(i,j,k) = 0.0
        end do
       end do
       numatomsInCategory(i,1) = 0.0
       numatomsInCategoryTotal(i,1) = 0.0
       do k=1,numBinsInHistogram
        CoordHist(i,k) = 0.0
       end do
       do k=1,maxnumatomtypes
	 do l=1,pmfMaxNumBins
	   pmfShellCounterHiRes(i,k,l) = 0.0
	   PMFHiRes(i,k,l) = 0.0
         end do
       end do
       do k=1,maxnumatomtypesEnv
        do l=1,pmfMaxNumBins
         pmfShellCounter(i,k,l) = 0.0
         PMFENERGY(i,k,l) = 0.0
        end do
       end do
      end do

      do k=1,numBinsInHistogram
       CoordHistAverage(k) = 0.0
      end do
 
      do i=1,maxnumatomtypesEnv
       numatomsInCategory(i,2) = 0.0
       numatomsInCategoryTotal(i,2) = 0.0
       do j = 1, maxnumatomtypes
        pair(j,i) = 0.00
        pairAvg(i) = 0.00
        EPair(j,i) = 0.00
       end do
      end do

      do i=1,maxnumatomtypesTop
       numatomsInCategory(i,3) = 0.0
       numatomsInCategoryTotal(i,3) = 0.0
      end do


      numatmTotal2 = 0
      numenvTotal = 0.0
      pdbsread = 0
      numchainused = 0
      totalInteraction = 0
      criticalCoordination = 13

c--------------------------------------------------------

      if (firstTimeThrough) then
       call ReadCombos (comboTypes,comboTypesComp,contactdir)
       comboTypes(maxNumCoord) = 1
       call ReadChiSqrTable (contactdir,ChiSqrTable,ChiSqrProb)
       if (dbexist) then
        call Read_PDBselect(dbconcat,pdbnames,pdbChain,numpdbnames,
     &  nHlx,nBta,thrsh)
       end if
      end if

c-------------------------------------------------------
c Check for existence of energy file corresponding to user-defined
c parameters.  If it exits (i.e. iocode = 0 ) then read in values
c and skip the derivation.

90    if (enexist) then
       print *,'File containing derived energies exists....'
       print *
       print *,'Skipping derivation....'
       print *
       if (.not.mfexist) then
        print *,'No misfolds to test...Nothing to do!!!'
        print *,'Stopping'
        stop
       end if
        call readenergy3_formatted(energyfile,
     &    Eword,comboTypes,contactdir)
        print *,'Energies read from file'
        goto 499
      else
       if (.not.dbexist) then 
        print *,'No energyfile, no database, no way of calculating
     & energies.'
        print *,'Stopping'
        stop
       end if
      end if

c-------------------------------------------------------

      print *
      print *,'Beginning Derivation....'
      do ii = 1, numpdbnames
       filename(1:11) = pdbnames(ii)(1:11)
       rejectflag(ii) = .false.

c-------------------------------------------------------
c Initialize all protein-specific variables before reading
c new pdb file

      numres = 0
      numatm = 0
      do i=1,maxnumatomtypes
       numatomsInCategory(i,1) = 0.0
      end do
      do i=1,maxnumatomtypesEnv
       numatomsInCategory(i,2) = 0.0
      end do
      do i=1,maxnumatomtypesTop
       numatomsInCategory(i,3) = 0.0
      end do
      do i=1,maxat
       multiconf(i) = ''
       unknown(i) = .false.
       do j=1,7
        atomData(i,j) = 0
       end do
       Coordination(i) = 0
       res(i) = ''
       atomname(i) = ''
       category(i,1) = 0
       category(i,2) = 0
       category(i,3) = 0
       do j=1,3
        coord(j,i) = 0.0
       end do
       do j=1,5
        connectivity(j,i) = 0
       end do
      end do
c----------------------------------------------------------
       reject = .false.
       reason = ''
       print *,'reading pdbfile : ',filename

       call pdbread_no_biomt_reduced (filename,databasedir)
       print *,'# of atoms = ',numatm
  
       badAtom = .false. 
       firsttime = .true.

       call pdbcheck (contactdir,filename,reason,firsttime)

       print *,'checked the pdbfile'

       if (reject) then
        print *,'Rejected : ',filename,'  ',reason
        rejectflag(ii) = .true.
        goto 222
       else 
       end if

      numchainused = numchainused + 1 

      call categorizeatoms (numatomtypesEnv,numatomtypesTop,
     &   categorylistname,contactdir,mapAtomTopology,firstTimeThrough,
     &   badAtom,category)
 

      if (badAtom) then
       call pdbcollapse(category)
      end if

      call findRelevantChain (pdbChain(ii))

      badAtom = .false.
  
      call pdbcheck (contactdir,filename,reason,firsttime)

      if (badAtom) then
       print *,'weve still got a bad atom in file ',filename
      end if

      call countcategories (category,numatomsInCategory)

c---------------------------------------------------------
c Calculate a running sum of the total number of atoms in each
c category 

      do j=1,numatomtypes
       numatmTotal2 = numatmTotal2 + numatomsInCategory(j,1)
       numatomsInCategoryTotal(j,1) = numatomsInCategoryTotal(j,1) +
     +     numatomsInCategory(j,1)
      end do
      do j=1,numatomtypesEnv
       numenvTotal = numenvTotal + numatomsInCategory(j,2)
       numatomsInCategoryTotal(j,2) = numatomsInCategoryTotal(j,2) +
     +     numatomsInCategory(j,2)
      end do
       do j=1,numatomtypesTop
       numatomsInCategoryTotal(j,3) = numatomsInCategoryTotal(j,3) +
     +     numatomsInCategory(j,3)
      end do

c-----------------------------------------------------------------
c Now calculate distances, relevancy of pairwise contact 
c (based in part on connectivity and atom type)

       do i1 = chainStart,chainEnd
        Coordination(i1) = 0
        if (.not.unknown(i1)) then
         atomData(i1,1) = category(i1,1)
         do i4=2,7
           atomData(i1,i4) = 0
         end do
         do i2 = 1,numatm
          if (.not.unknown(i2)) then
          call relevantDistance2(i1,i2,numskip,relevant)
           if (relevant) then
            tempdistance = Distance(coord(1,i1),coord(2,i1),
     &       coord(3,i1),coord(1,i2),coord(2,i2),coord(3,i2))
c------ -------------------------------------------
            call Make_Bin_Decision(tempdistance,pmfMaxdis,pmfMindis,
     &         pmfStep,pmfShell)
            if (pmfShell.eq.1) then
              print *,'first_shell!  ',filename,i1,i2
            endif
            if (pmfShell.ne.999) then
             pmfShellCounter(category(i1,1),
     &                       category(i2,2),pmfShell) =
     &       pmfShellCounter(category(i1,1),
     &                       category(i2,2),pmfShell) + 1.0
             pmfShellCounterHiRes(category(i1,1),
     &                       category(i2,1),pmfShell) =
     &       pmfShellCounterHiRes(category(i1,1),
     &                       category(i2,1),pmfShell) + 1.0
            endif
            if (tempdistance.le.criticalmaxdistance.and.
     &          tempdistance.gt.criticalmindistance) then
              Coordination(i1) = Coordination(i1) + 1
              atomData(i1,2) = atomData(i1,2) + 1
              atomData(i1,(category(i2,2)+2)) =
     &          atomData(i1,(category(i2,2)+2)) + 1
            endif
           endif
          endif
         enddo
        else
c        what to do if atom at i1 is unknown         
         do i4 = 1,7
           atomData(i1,i4) = 999
         end do
        endif
       end do


c Now we'll group the atoms based on atomtype(1-20),coordination(0-13),
c and the #s of atoms in each of the four environment variables.

       do i1 = chainStart,chainEnd
        if (.not.unknown(i1)) then
         if (atomData(i1,2).le.criticalCoordination) then
          do i2 = 1,comboTypes(atomData(i1,2)+1)
           if (atomData(i1,3).eq.
     &        comboTypesComp(atomData(i1,2)+1,i2,1)) then
            if (atomData(i1,4).eq.
     &         comboTypesComp(atomData(i1,2)+1,i2,2)) then
             if (atomData(i1,5).eq.
     &          comboTypesComp(atomData(i1,2)+1,i2,3)) then
              if (atomData(i1,6).eq.
     &           comboTypesComp(atomData(i1,2)+1,i2,4)) then
                event(atomData(i1,1),atomData(i1,2)+1,i2) = 
     &           event(atomData(i1,1),atomData(i1,2)+1,i2) + 1.0
                alteredEvent(atomData(i1,1),atomData(i1,2)+1,
     &           i2) = alteredEvent(atomData(i1,1),
     &          atomData(i1,2)+1,i2) + 1.0
              end if
             end if
            end if
           end if
          end do
         else 
          event(atomData(i1,1),15,1) = 
     &     event(atomData(i1,1),15,1) + 1.0
          alteredEvent(atomData(i1,1),15,1) =
     &     alteredEvent(atomData(i1,1),15,1) + 1.0
         end if
        end if
       end do

c- repeat counting for pairwise calculations.....
c- change range of loop variables (i2, specifically) to
c- avoid double counting

      do i1 = chainStart,chainEnd
       if (.not.unknown(i1)) then
        do i2 = 1,numatm
         if (.not.unknown(i2)) then
          call relevantDistance2(i1,i2,numskip,relevant)
          if (relevant) then
           tempdistance = Distance(coord(1,i1),coord(2,i1),
     &       coord(3,i1),coord(1,i2),
     &       coord(2,i2),coord(3,i2))
           if (tempdistance.lt.criticalmaxdistance.and.
     &        tempdistance.gt.criticalmindistance) then
              pair(category(i1,1),category(i2,2)) = 
     &        pair(category(i1,1),category(i2,2)) + 1
           endif
          endif
         endif
        enddo
       endif
      enddo

c- Calculating the number of environment atoms in the sphere
c- of influence surrounding each atomtype/coordination number
c- combination....

      do i1 = chainStart,chainEnd
       if (.not.unknown(i1)) then
       if (atomData(i1,2).ge.15) then
        atomData(i1,2) = 14
       end if
       do i2 = 1,4
        ThingsinEnv(atomData(i1,1),atomData(i1,2)+1,i2) =
     &   ThingsinEnv(atomData(i1,1),atomData(i1,2)+1,i2) +
     &   atomData(i1,i2+2)
       end do
       numInCatCoord(atomData(i1,1),atomData(i1,2)+1) =
     &  numInCatCoord(atomData(i1,1),atomData(i1,2)+1) +
     &  1.00 
       end if
      end do


c---------------------------------------------------------------
c Calculation of Histograms for Coordination Numbers
c based on atom type, and based on topological type

      do i1 = chainStart, chainEnd
      if (.not.unknown(i1)) then
       if (Coordination(i1).ge.0.and.Coordination(i1)
     +   .le.(numBinsInHistogram-2)) then
          CoordHist(category(i1,1),Coordination(i1)+1) = 
     +     CoordHist(category(i1,1),Coordination(i1)+1) + 1
       else
        CoordHist(category(i1,1),numBinsInHistogram) = 
     +   CoordHist(category(i1,1),numBinsInHistogram) + 1
       end if
      end if
      end do

c Go back and see if there another pdb file to read in!

      pdbsread = pdbsread + 1
      firstTimeThrough = .false.

c End of major loop
c-----------------------------------------------------------
222   end do

      call colpdblist (rejectflag,pdbnames,pdbChain,numpdbnames,
     &    nHlx,nBta,thrsh)

      print *,'Number of chains used in derivation = ',
     &  numchainused,numpdbnames

      do i=1,numpdbnames
       print *,pdbnames(i)
       nHlxTot = nHlxTot + nHlx(i)
       nBtaTot = nBtaTot + nBta(i)
      end do

      print *
      print *,'# helical residues = ',nHlxTot
      print *,'# beta-sheet residues  = ',nBtaTot
      print *
      
      numatmTotal2Real = float(numatmTotal2)


c-------------------------
c pair(i,j) is the number of pairwise contacts seen
c in the database between atomtype i and environment
c atom type j.
c
c N1(i) is the sum of all the pairwise contacts between
c an atom type i and any other environment atom type
c
c p1(i,j) is the relative probability of finding a
c pairwise interaction between an atom of type i and an
c enviroment atom j and is defined as
c
c p1(i,j) = pair(i,j)/N1(i)
c--------------------------

      do i1 = 1,numatomtypes
       do i2 = 1,numatomtypesEnv
        N1(i1) = N1(i1) + pair(i1,i2)
       end do
      end do

      do i1 = 1,numatomtypes
       do i2 = 1,numatomtypesEnv
        p1(i1,i2) = pair(i1,i2)/N1(i1)
       end do
      end do


c-----------------------------------------------------
c Calculation of CoordHistaverage

      do i1 = 1,numBinsInHistogram
       do i2 = 1,numatomtypes
        CoordHistAverage(i1) = CoordHistAverage(i1) +
     +    CoordHist(i2,i1)
       end do
      end do

c-----------------------------------------------------
c Now, compare the coordination histograms on an atomtype
c basis to the average coordination histogram.  The
c metric here will be analogous to an RMSD over the
c coordination number bins after the average coordination
c histogram has been normalized to contain the same number
c of counts as the histogram to which we would like to compare
c the average.

      do i1 = 1,numatomtypes
c First, calculate a normalization constant
       rms = 0.0
       normConst = numatomsInCategoryTotal(i1,1) / numatmTotal2Real
       do i2 = 1,maxNumCoord
        CoordHistAverageNorm(i2) = CoordHistAverage(i2) *
     +    normConst
        tempreal = CoordHist(i1,i2)/CoordHistAverageNorm(i2)
        if ((CoordHist(i1,i2).ne.0.0).and.
     +   (CoordHistAverageNorm(i2).ne.0.0)) then
          singleBodyEnergy(i1,i2) =
     +     (-1)* log(CoordHist(i1,i2)/
     +     CoordHistAverageNorm(i2))
        end if
        if (CoordHist(i1,i2).eq.0.0) then
         singleBodyEnergy(i1,i2) = 
     +    (-1) * log(1/CoordHistAverageNorm(i2))
        end if
        if (CoordHistAverageNorm(i2).eq.0.0) then
          singleBodyEnergy(i1,i2) = 0.0
        end if
        rms = rms + (CoordHistAverageNorm(i2) -
     +     CoordHist(i1,i2))**2
       end do
       rms = sqrt(rms/maxNumCoord)
      end do

      open (unit=4,file='CoordHist.dat')
      write (4,12) 'Coordination # Histogram (Atomtype)'
      write (4,12) 'X axis is atomtype, Y axis is Coordination #'
      tempInteger = numBinsInHistogram
      write (4,36) 'Avg.',(j,j=1,numatomtypes)
      do i1 = 1,tempInteger
       write (4,14) i1-1,CoordHistAverage(i1),
     +   (CoordHist(j,i1),j=1,numatomtypes)
      end do
      close(4)

      print 12, 'Coordination # Histogram (Atomtype)'
      print 12, 'X axis is atomtype, Y axis is Coordination #'
      tempInteger = numBinsInHistogram
      print 36, 'Avg.',(j,j=1,numatomtypes)
      do i1 = 1,tempInteger
       print 14, i1-1,CoordHistAverage(i1),
     +   (CoordHist(j,i1),j=1,numatomtypes)
      end do

c Calculation of the actual probability of seeing a particular 
c coordination number about a given atom type

      do i2 = 1,maxNumCoord
        do i1 = 1,numatomtypes
         if (numatomsInCategoryTotal(i1,1).ne.0.0) then
          ActualProbCoord(i1,i2) = 
     &    CoordHist(i1,i2) / 
     &    numatomsInCategoryTotal(i1,1)
         else
          ActualProbCoord(i1,1) = 0.0
         end if
        end do
      end do

c Calculate the number of word bins that are significantly populated
c (i.e. >=10 counts) and the number that are not (i.e. <10 counts)

      NumberSignificantlyPopulated = 0
      NumberNotSignificantlyPopulated = 0
      do i1 = 1,numatomtypes
       do i2 = 1,maxNumCoord
        do i3 = 1,comboTypes(i2)
         SumTotalEvents = SumTotalEvents + event(i1,i2,i3)
         if (event(i1,i2,i3).ge.6.0) then
          NumberSignificantlyPopulated = 
     &     NumberSignificantlyPopulated + 1.0
         else
          NumberNotSignificantlyPopulated =
     &     NumberNotSignificantlyPopulated + 1.0
         end if 
        end do
       end do
      end do

c Those events which fall into a coordination number of 14 or greater fall
c into the category of Unusable Events because we cannot enumerate the
c nondegenerate words with this large a coordination number

      do i1 = 1,numatomtypes
       SumUnusableEvents = SumUnusableEvents + event(i1,15,1)
      end do
  
      SumTotalEvents = SumTotalEvents + SumUnusableEvents


c Calculate sum of Atoms in each atomtype/coordination number
c class

      do i1 = 1,numatomtypes
       do i2 = 1,maxNumCoord
        do i3 = 1,4
         SumThingsinEnv(i1,i2) = SumThingsinEnv(i1,i2) +
     &    ThingsinEnv(i1,i2,i3)
        end do
       end do
      end do

c NEW VERSION of calculation of probCatCoord
c The theoretical probability of drawing an environment atom of 
c type A,B,C,or D given that the atom type has been chosen, as
c well as the coordination number

      do i1 = 1,numatomtypes
       do i2 = 1,maxNumCoord
        do i3 = 1,4
         if (SumThingsinEnv(i1,i2).ne.0) then
          probCatCoord(i1,i2,i3) =
     &     ThingsinEnv(i1,i2,i3)/
     &     SumThingsinEnv(i1,i2)
         else
          probCatCoord(i1,i2,i3) = 0.0
         end if
        end do
       end do
      end do

c Print out values of probCatCoord

      print *
      print *,'-----------------------------'
      print *,'P(envatom|atomtype,coord)'
      print *,'First field: atomtype'
      print *,'Second field: coordination'
      print *,'Fields 3-6: Probabilities assuming pairwise ind.'
      do i1 = 1,numatomtypes
       do i2 = 1,maxNumCoord
        print 77,i1,i2-1,(probCatCoord(i1,i2,i3),i3=1,4)
       end do
      end do
      print *,'-----------------------------'
      print *

c Now calculate the the theoretical probability of a given environment given
c a coordination number:  based on A,B,C,D probabilities calculated above.

      do i1 = 1,numatomtypes
       do i2 = 1,maxNumCoord
        do i3 = 1,comboTypes(i2)
        Abool = .true.
        Bbool = .true.
        Cbool = .true.
        Dbool = .true.
        numA = comboTypesComp(i2,i3,1)
        numB = comboTypesComp(i2,i3,2)
        numC = comboTypesComp(i2,i3,3)
        numD = comboTypesComp(i2,i3,4)
        Pmult = Fact(numA+numB+numC+numD)/
     &     (Fact(numA)*Fact(numB)*Fact(numC)*Fact(numD))

c This code block assigns the probability of picking a particular
c environment atom to probCatCoord
c
c Used when testing misfolds

        PA = probCatCoord(i1,i2,1)**numA
        PB = probCatCoord(i1,i2,2)**numB
        PC = probCatCoord(i1,i2,3)**numC
        PD = probCatCoord(i1,i2,4)**numD  

        if (Pmult.eq.0) then
         TheoreticalProbWordGCoord(i1,i2,i3) = 0.0
        end if 

        if (numA.eq.0) then
         Abool = .false.
        end if
        if (numB.eq.0) then
         Bbool = .false.
        end if
        if (numC.eq.0) then
         Cbool = .false.
        end if
        if (numD.eq.0) then
         Dbool = .false.
        end if

        TheoreticalProbWordGCoord(i1,i2,i3) = 1.0
        if (Abool) then
         TheoreticalProbWordGCoord(i1,i2,i3) = 
     +     TheoreticalProbWordGCoord(i1,i2,i3) * PA
        end if
        if (Bbool) then
         TheoreticalProbWordGCoord(i1,i2,i3) = 
     +     TheoreticalProbWordGCoord(i1,i2,i3) * PB
        end if
        if (Cbool) then
         TheoreticalProbWordGCoord(i1,i2,i3) = 
     +     TheoreticalProbWordGCoord(i1,i2,i3) * PC
        end if
        if (Dbool) then
         TheoreticalProbWordGCoord(i1,i2,i3) = 
     +     TheoreticalProbWordGCoord(i1,i2,i3) * PD
        end if
        if ((.not.Abool).and.(.not.Bbool).and.(.not.Cbool).and.
     +     (.not.Dbool).and.(i2.ne.1)) then
         TheoreticalProbWordGCoord(i1,i2,i3) = 0.0
        end if
        if ((.not.Abool).and.(.not.Bbool).and.(.not.Cbool).and.
     +     (.not.Dbool).and.(i2.eq.1)) then
         TheoreticalProbWordGCoord(i1,i2,i3) = 1.0
        end if

         TheoreticalProbWordGCoord(i1,i2,i3) = 
     +     TheoreticalProbWordGCoord(i1,i2,i3) * Pmult

c -----
c In this version of the program, the theoretical probability of
c seeing a given cluster (or word) given an atomtype and
c coordination is calculated based on the assumption that
c pairwise interaction events are independent
c -----

        end do
       end do
      end do

c Calculate the 'actual' probability of a given microstate from the
c events data
c
c    p_actual = #events in microstate / total # events in coordination
c
c Also, calculate overallSumEvents

      do i1 = 1,numatomtypes
       do i2 = 1,maxNumCoord
         sumEvents(i1,i2) = 0.0
        do i3 = 1,comboTypes(i2)
         sumEvents(i1,i2) = sumEvents(i1,i2) + 
     &     event(i1,i2,i3)
         overallSumEvents = overallSumEvents + 
     &    event(i1,i2,i3)
        end do
       end do
      end do
      do i1 = 1,numatomtypes
       overallSumEvents = overallSumEvents + event(i1,15,1)
      end do

      do i1 = 1,numatomtypes
       do i2 = 1,maxNumCoord
        do i3 = 1,comboTypes(i2)
         if (sumEvents(i1,i2).ne.0.0) then
          ActualProbWordGCoord(i1,i2,i3) =
     &     event(i1,i2,i3) /
     &     sumEvents(i1,i2)
        else
         ActualProbWordGCoord(i1,i2,i3) = 0.0
        end if
        end do
       end do
      end do

      print *,'overallSumEvents = ',overallSumEvents

      do i1 = 1,numatomtypes
       do i2 = 1,maxNumCoord
        TheoreticalProbCoord(i1,i2) = 
     +   CoordHistAverage(i2) / numatmTotal2Real
       end do
      end do

c*************************************************
c  NEW VERSION
c Calculate ActualProbWord,TheoreticalProbWord

      do i1 = 1,numatomtypes
       do i2 = 1,maxNumCoord
        do i3 = 1,comboTypes(i2)
         ActualProbWord(i1,i2,i3) =
     &      ActualProbWordGCoord(i1,i2,i3)
         TheoreticalProbWord(i1,i2,i3) =
     &      TheoreticalProbWordGCoord(i1,i2,i3)
        end do
       end do
      end do

c**************************************************
c NEW VERSION
c Calculate ExpectedEvents assuming random sampling

      do i1 = 1,numatomtypes
       do i2 = 1,maxNumCoord
        do i3 = 1,comboTypes(i2)
         ExpectedEvents(i1,i2,i3) =
     +     TheoreticalProbWord(i1,i2,i3) *
     +     numInCatCoord(i1,i2)
        end do
       end do
      end do

      
      open (unit=4,file="ActualProb.dat")
      do i1 = 1,numatomtypes
       tempreal = 0.0
       do i2 = 1,maxNumCoord
        do i3 = 1,comboTypes(i2)
         write (4,33) i1,i2,i3,ActualProbWord(i1,i2,i3)
         tempreal = tempreal + ActualProbWord(i1,i2,i3)
        end do
       end do
       write (4,33) i1,0,0,tempreal
      end do
      close(4)
      open (unit=4,file="TheoreticalProb.dat")
      do i1 = 1,numatomtypes
       tempreal = 0.0
       do i2 = 1,maxNumCoord
        do i3 = 1,comboTypes(i2)
         write (4,33) i1,i2,i3,TheoreticalProbWord(i1,i2,i3)
         tempreal = tempreal + TheoreticalProbWord(i1,i2,i3)
        end do
       end do
       write (4,33) i1,0,0,tempreal
      end do
      close(4)

c-------------------------------------------------------------------
c  Significance Calculation
c-------------------------------------------------------------------


c Calculate ChiSquared statistic for each atomtype-coordination
c as deviation from random for all cases

      call CalcSig_pure(numatomtypes,comboTypes,ChiSqrTable,
     +  ChiSqrProb,ExpectedEvents,event,ChiSquared,Significance,
     +  SignificanceBool,SignificanceHistogram,PerSigGT95)

      do i1 = 1,numatomtypes
       do i2 = 1,15
        do i3 = 1,comboTypes(i2)
         if (SignificanceBool(i1,i2,i3).eq.1) then
          TotalSignificantEvents = TotalSignificantEvents + 
     +      event(i1,i2,i3)
         else
          TotalNonSignificantEvents = TotalNonSignificantEvents + 
     +      event(i1,i2,i3)
         end if
        end do
       end do
      end do

      PercentageSignificantEvents = (TotalSignificantEvents/
     + (TotalSignificantEvents + TotalNonSignificantEvents))*100

c-------------------------------------------------
       print *,'Significance Data for Energy Calculation'
       print *,'Significance Histogram over full data set'
       print 29,(SignificanceHistogram(i1),i1 = 1,14)

       print *,'PerSigGT95 = ',PerSigGT95
       print *,'SumTotalEvents = ',SumTotalEvents
       print *,'SumUnusableEvents = ',SumUnusableEvents
       print *
       print *,'Ratio Unusable/Total Events = ',
     +   SumUnusableEvents/SumTotalEvents

       print *,'# of Bins significantly populated = ',
     +   NumberSignificantlyPopulated
       print *,'# of Bins not significantly populated = ',
     +   NumberNotSignificantlyPopulated
       print *,'Ratio = ',NumberSignificantlyPopulated/
     +   (NumberSignificantlyPopulated +
     +    NumberNotSignificantlyPopulated)

       print *
       print *
       print *,'TotalSignificantEvents',TotalSignificantEvents
       print *,'TotalNonSignficantEvents',TotalNonSignificantEvents
       print *,'PercentageSignficantEvents',PercentageSignificantEvents

c**********************************************
c if events=0, then add 1 to events for upcoming
c energy calculation

c The value 4.50 is the critical value of the mean
c given an observed value of 1, where the significance is 
c greater than 95% (Assures that our 0 to 1 approximation
c does not alter the significance of the counting statistics

       do i1 = 1,numatomtypes
        do i2 = 1,maxNumCoord
         do i3 = 1,comboTypes(i2)
c! Chris is this valid?
          if (event(i1,i2,i3).eq.0.0) then
            alteredEvent(i1,i2,i3) = 
     +      alteredEvent(i1,i2,i3) + 1.0
          end if
          if (ExpectedEvents(i1,i2,i3).lt.4.5) then
           SignificanceBool(i1,i2,i3) = 0
          end if
         end do
        end do
       end do

c***********************************************
c Energy Calculation

c      do i1 = 1,numatomtypes
c       do i2 = 1,maxNumCoord
c        do i3 = 1,comboTypes(i2)
c         if (SignificanceBool(i1,i2,i3).eq.1) then
c          Eword(i1,i2,i3) = 
c     +     (-1)*log(alteredEvent(i1,i2,i3)/
c     +     ExpectedEvents(i1,i2,i3))
c         else
c          Eword(i1,i2,i3) = singleBodyEnergy(i1,i2)
c         end if
c        end do
c       end do
c      end do 

c***********************************************
c Energy Calculation
c
c      do i1 = 1,numatomtypes
c       do i2 = 1,maxNumCoord
c        do i3 = 1,comboTypes(i2)
c         if (alteredEvent(i1,i2,i3).eq.1.0.and.
c     &       event(i1,i2,i3).eq.0.0) then
c                Eword(i1,i2,i3) = 5.0
c		 else
c          Eword(i1,i2,i3) = 
c     &     (-1)*log(alteredEvent(i1,i2,i3)/
c     &     ExpectedEvents(i1,i2,i3))
c         end if
c		 if (i3.eq.maxNumCoord) then
c           Eword(i1,i2,i3) =
c     &        singleBodyEnergy(i1,i2)
c         end if	
c        end do
c       end do
c      end do 
c***********************************************
c Energy Calculation

      do i1 = 1,numatomtypes
       do i2 = 1,maxNumCoord
        do i3 = 1,comboTypes(i2)
         if (TheoreticalProbWord(i1,i2,i3).eq.0.0.and.
     &      ActualProbWord(i1,i2,i3).eq.0.0) then
          Eword(i1,i2,i3) = singleBodyEnergy(i1,i2)
		 endif
         if (TheoreticalProbWord(i1,i2,i3).ne.0.0.and.
     &      ActualProbWord(i1,i2,i3).eq.0.0) then
c Some very, very bad large number
          Eword(i1,i2,i3) = 10.0
		 endif
         if (TheoreticalProbWord(i1,i2,i3).eq.0.0.and.
     &      ActualProbWord(i1,i2,i3).ne.0.0) then
c This is a total hack!  this should never, ever occur
          Eword(i1,i2,i3) = -10.0
		 endif
         if (TheoreticalProbWord(i1,i2,i3).ne.0.0.and.
     &      ActualProbWord(i1,i2,i3).ne.0.0) then
          Eword(i1,i2,i3) = 
     &     (-1)*log(ActualProbWord(i1,i2,i3)/
     &     TheoreticalProbWord(i1,i2,i3))
		 endif
         if (i2.eq.maxNumCoord) then
           Eword(i1,i2,i3) =
     &        singleBodyEnergy(i1,i2)
         end if	
        end do
       end do
      end do 

c*************************************************
c write out Expected,Actual for all significant events
      open (unit=3,file='raw_counts',status='NEW')
      do i1 = 1,numatomtypes 
       do i2 = 1,maxNumCoord
        do i3 = 1,comboTypes(i2)
         if (SignificanceBool(i1,i2,i3).eq.1) then
          write (3,78) '**',i1,i2-1,i3,'|',
     &     comboTypesComp(i2,i3,1),comboTypesComp(i2,i3,2),
     &     comboTypesComp(i2,i3,3),comboTypesComp(i2,i3,4),
     &     '|',event(i1,i2,i3),'/',ExpectedEvents(i1,i2,i3),
     &     Eword(i1,i2,i3)
         else
          write (3,78) '  ',i1,i2-1,i3,'|',
     &     comboTypesComp(i2,i3,1),comboTypesComp(i2,i3,2),
     &     comboTypesComp(i2,i3,3),comboTypesComp(i2,i3,4),
     &     '|',event(i1,i2,i3),'/',ExpectedEvents(i1,i2,i3),
     &     Eword(i1,i2,i3)
         end if
        end do
       end do
      end do
      close (3)

c skips TwentyBest code
c      goto 299
      
c-------------------------------------------------
c TwentyBest,TwentyWorst

      call BubbleSort(Eword,TwentyBestE,TwentyBestEMap,TwentyWorstE,
     & TwentyWorstEMap,numatomtypes,comboTypes,SignificanceBool)

      print *,'Contents of TwentyBestEMap'
      do i=1,20
       do j=1,20
        print *,TwentyBestEMap(i,j,1),TwentyBestEMap(i,j,2)
       end do
      end do

      do i = 1,20
       print *
       print *,'Twenty Best for atom type ',i
       print *,'--------------------------------------------------------
     &---'
       print *,'  energy   C# cluster N    C    O    S     actual  expec
     &ted'
       print *,'--------------------------------------------------------
     &---'
       do j=1,20
        print 44,TwentyBestE(i,j),
     &   TwentyBestEMap(i,j,1)-1,
     &   TwentyBestEMap(i,j,2),
     &   comboTypesComp(TwentyBestEMap(i,j,1),
     &   TwentyBestEMap(i,j,2),1),
     &   comboTypesComp(TwentyBestEMap(i,j,1),
     &   TwentyBestEMap(i,j,2),2),
     &   comboTypesComp(TwentyBestEMap(i,j,1),
     &   TwentyBestEMap(i,j,2),3),
     &   comboTypesComp(TwentyBestEMap(i,j,1),
     &   TwentyBestEMap(i,j,2),4),
     &   event(i,TwentyBestEMap(i,j,1),
     &   TwentyBestEMap(i,j,2)),
     &   ExpectedEvents(i,
     &   TwentyBestEMap(i,j,1),
     &   TwentyBestEMap(i,j,2))
       end do
       print *
       print *,'Twenty Worst for atom type ',i
       print *,'--------------------------------------------------------
     &---'
       print *,'  energy   C# cluster N    C    O    S     actual  expec
     &ted'
       print *,'--------------------------------------------------------
     &---'
       do j=1,20
        print 44,TwentyWorstE(i,j),
     &   TwentyWorstEMap(i,j,1)-1,
     &   TwentyWorstEMap(i,j,2),
     &   comboTypesComp(TwentyWorstEMap(i,j,1),
     &   TwentyWorstEMap(i,j,2),1),
     &   comboTypesComp(TwentyWorstEMap(i,j,1),
     &   TwentyWorstEMap(i,j,2),2),
     &   comboTypesComp(TwentyWorstEMap(i,j,1),
     &   TwentyWorstEMap(i,j,2),3),
     &   comboTypesComp(TwentyWorstEMap(i,j,1),
     &   TwentyWorstEMap(i,j,2),4),
     &   event(i,TwentyWorstEMap(i,j,1),
     &   TwentyWorstEMap(i,j,2)),
     &   ExpectedEvents(i,
     &   TwentyWorstEMap(i,j,1),
     &   TwentyWorstEMap(i,j,2))
       end do
      end do
299   continue

c skip examples of best and worst
c      goto 335

c------------------------------------------------------- 
c Here, we're going to go through the protein database again to pull out
c examples of the twenty best and twenty worst.  Unfortunately, this involves
c an almost total repetition of the first section of code for reading in the
c PDB files and categorizing the atoms, etc. (A cut and paste job!)
c-------------------------------------------------------

      print *,'Examples of best and worst'

c      do ii=1,pdbsread
c following line is a hack...
      do ii=1,20
       pdbname(1:11) = pdbnames(ii)(1:11)

c      open (unit=4,file = mfconcat)
c300   read (4,12,end = 334) pdbname
c-------------------------------------------------------
c Initialize all protein-specific variables before reading
c new pdb file

      numres = 0
      numatm = 0
      do i=1,maxat
       do j=1,7
        atomData(i,j) = 0
       end do
       atheader(i) = ''
       multiconf(i) = ''
       chain(i) = ''
       Coordination(i) = 0
       res(i) = '   '
       atomname(i) = '    '
       category(i,1) = 0
       category(i,2) = 0
       category(i,3) = 0
       do j=1,3
        coord(j,i) = 0.0
       end do
       do j=1,5
        connectivity(j,i) = 0
       end do
      end do

c----------------------------------------------------------

       reject = .false.
       reason = ''

       call pdbread_no_biomt_reduced (pdbname,misfolddir)
       do i=1,maxat
       	 Bfac(i) = 0.0
       end do

       print *,'reading pdbfile: ',pdbname

       badAtom = .false.

       call pdbcheck (contactdir,pdbname,reason,firsttime)

       if (reject) then
        print *,'Rejected : ',pdbname,'  ',reason
        rejectflag(ii) = .true.
        goto 333
       else
        print *,'Read : ',pdbname
       end if

       call categorizeatoms (numatomtypesEnv,numatomtypesTop,
     &   categorylistname,contactdir,mapAtomTopology,firstTimeThrough,
     &   badAtom,category)


      if (badAtom) then
       call pdbcollapse (category)
      end if

c      call findRelevantChain (pdbChain(ii))
c In lieu of findRelevantChain in this version, manually
c set chainStart and chainEnd
       chainStart = 1
       chainEnd = numatm


c-----------------------------------------------------------------
c Now calculate distances, relevancy of pairwise contact
c (based in part on connectivity and atom type)

      do i1 = chainStart,chainEnd
      if (.not.unknown(i1)) then
       atomData(i1,1) = category(i1,1)
       do i2 = 1,numatm
        if (.not.unknown(i2)) then
        call relevantDistance2(i1,i2,numskip,relevant)
        if (relevant) then
         tempdistance = Distance(coord(1,i1),coord(2,i1),
     &     coord(3,i1),coord(1,i2),
     &     coord(2,i2),coord(3,i2))
         if (tempdistance.lt.criticalmaxdistance.and.
     &       tempdistance.gt.criticalmindistance) then
           Coordination(i1) = Coordination(i1) + 1
           atomData(i1,2) = atomData(i1,2) + 1
           atomData(i1,(category(i2,2)+2)) =
     &       atomData(i1,(category(i2,2)+2)) + 1
         end if
        end if
       end if
       end do
      end if
      end do

c---------
c Does this file contain any of the best (or worst) environments???
c If so, print out the atom number, the atom name, the residue name,
c and the Best/Worst map in question....

      BestEnergy = 20.0
      WorstEnergy = -20.0
      do i1 = chainStart,chainEnd
      if (.not.unknown(i1)) then
       do i2 = 1,20
        i3 = TwentyBestEMap(atomData(i1,1),i2,1)
        i4 = TwentyBestEMap(atomData(i1,1),i2,2)
        if ((i3-1).eq.Coordination(i1)) then
         if (comboTypesComp(i3,i4,1).eq.atomData(i1,3)) then
          if (comboTypesComp(i3,i4,2).eq.atomData(i1,4)) then
           if (comboTypesComp(i3,i4,3).eq.atomData(i1,5)) then
            if (comboTypesComp(i3,i4,4).eq.atomData(i1,6)) then
             tempreal = Eword(atomData(i1,1),
     &                        atomData(i1,2)+1,i4)
             atomEnergy(i1) = tempreal
             if (tempreal.lt.BestEnergy) then
               BestEnergy = tempreal
             end if
            end if
           end if
          end if
         end if
        end if
       end do
       do i2 = 1,20
        i3 = TwentyWorstEMap(atomData(i1,1),i2,1)
        i4 = TwentyWorstEMap(atomData(i1,1),i2,2)
        if ((i3-1).eq.Coordination(i1)) then
         if (comboTypesComp(i3,i4,1).eq.atomData(i1,3)) then
          if (comboTypesComp(i3,i4,2).eq.atomData(i1,4)) then
           if (comboTypesComp(i3,i4,3).eq.atomData(i1,5)) then
            if (comboTypesComp(i3,i4,4).eq.atomData(i1,6)) then
             tempreal = Eword(atomData(i1,1),
     &                        atomData(i1,2)+1,i4)
             atomEnergy(i1) = tempreal
             if (tempreal.gt.WorstEnergy) then
               WorstEnergy = tempreal
             end if
            end if
           end if
          end if
         end if
        end if
       end do
      end if
      end do

      call bsort(atomEnergy,sortedMap,numatm)

      do i1 = chainStart,chainEnd
       Bfac(i1) = 25.0
      end do

      numtoflag = 10

      do i1 = 1,numtoflag
       Bfac(sortedMap(i1)) = 50.0
       Bfac(sortedMap(numatm+1 - i1)) = 1.0
      end do


      print *,'pdbname = ',pdbname
      do i1 = chainStart,chainEnd
      if (.not.unknown(i1)) then
       do i2 = 1,20
        i3 = TwentyBestEMap(atomData(i1,1),i2,1)
        i4 = TwentyBestEMap(atomData(i1,1),i2,2)
        if ((i3-1).eq.Coordination(i1)) then
         if (comboTypesComp(i3,i4,1).eq.atomData(i1,3)) then
          if (comboTypesComp(i3,i4,2).eq.atomData(i1,4)) then
           if (comboTypesComp(i3,i4,3).eq.atomData(i1,5)) then
            if (comboTypesComp(i3,i4,4).eq.atomData(i1,6)) then
             tempreal = Eword(atomData(i1,1),
     &                        atomData(i1,2)+1,i4)
             print 71,pdbname,'Best',tempreal,category(i1,1),i2,
     &        Coordination(i1),res(i1),trueResnum(i1),chain(i1),
     &        atomname(i1)
c             Bfac(i1) =  (tempreal-BestEnergy)/
c     &           (WorstEnergy-BestEnergy) * 50
            end if
           end if
          end if
         end if
        end if
       end do
       do i2 = 1,20
        i3 = TwentyWorstEMap(atomData(i1,1),i2,1)
        i4 = TwentyWorstEMap(atomData(i1,1),i2,2)
        if ((i3-1).eq.Coordination(i1)) then
         if (comboTypesComp(i3,i4,1).eq.atomData(i1,3)) then
          if (comboTypesComp(i3,i4,2).eq.atomData(i1,4)) then
           if (comboTypesComp(i3,i4,3).eq.atomData(i1,5)) then
            if (comboTypesComp(i3,i4,4).eq.atomData(i1,6)) then
             tempreal = Eword(atomData(i1,1),
     &                        atomData(i1,2)+1,i4)
             print 71,pdbname,'Worst',tempreal,category(i1,1),i2,
     &        Coordination(i1),res(i1),trueResnum(i1),chain(i1),
     &        atomname(i1)
c             Bfac(i1) =  (tempreal-BestEnergy)/
c     &           (WorstEnergy-BestEnergy) * 50
            end if
           end if
          end if
         end if
        end if
       end do
      end if
      end do

c      numbestnearhetero = 0.0
c      numworstnearhetero = 0.0
c      do i1 = chainStart,chainEnd
c      if (.not.unknown(i1)) then
c       do i2 = 1,numatm
c        if (hetatm(i2)) then
c        if (bestarray(i1).eq.1) then
c         tempdistance = Distance(coord(1,i1),coord(2,i1),
c     &     coord(3,i1),coord(1,i2),
c     &     coord(2,i2),coord(3,i2))
c         if (tempdistance.lt.criticalmaxdistance.and.
c     &       tempdistance.gt.criticalmindistance) then
c           numbestnearhetero = numbestnearhetero + 1.0
c         end if
c        else
c         if (worstarray(i1).eq.1) then
c         tempdistance = Distance(coord(1,i1),coord(2,i1),
c     &     coord(3,i1),coord(1,i2),
c     &     coord(2,i2),coord(3,i2))
c          if (tempdistance.lt.criticalmaxdistance.and.
c     &       tempdistance.gt.criticalmindistance) then
c           numworstnearhetero = numworstnearhetero + 1.0
c          end if
c         end if
c        end if
c        end if
c       end do
c      end if
c      end do 
c
c      print *,'# of Twenty best near a heteroatom  ',
c     &  numbestnearhetero
c      print *,'# of Twenty worst near a heteroatom  ',
c     &  numworstnearhetero


c      call pdbwrite(pdbname)

333   end do

334   continue
      close (4)

335   continue
c End of TwentyBestCode


      print *,'ActualProbWordGCoord'
      do i2 = 1, numatomtypes
       do i1=1,5
        print 68,(ActualProbWordGCoord(i2,i1,k),
     &   k=1,comboTypes(i1))
       end do
      end do

      print *,'TheoreticalProbWordGCoord'
      do i2=1,numatomtypes
       do i1 = 1,5
        print 68,(TheoreticalProbWordGCoord(i2,i1,k),
     &   k=1,comboTypes(i1))
       end do 
      end do

      print *,'Significance'
      do i2 =1,numatomtypes
       do i1 = 1,5
        print 68,(Significance(i2,i1,k),k=1,comboTypes(i1))
       end do
      end do 

c---------------------------------------------------
c Calculate the fraction of each atom type

      do i=1,numatomtypes
        fraction(i) = numatomsInCategoryTotal(i,1) /
     &   numatmTotal2Real
      end do

c------------------------------------------------
c Print out some info for debugging purposes

      print *,'numatmTotal2 = ',numatmTotal2
      print *
      print *,'total numatoms in catgory'
      print *,'category','total'
      do j=1,numatomtypes
       print 43,j,numatomsInCategoryTotal(j,1)
      end do
      print *
      print *,'Environment','total'
      do j=1,numatomtypesEnv
       print 43,j,numatomsInCategoryTotal(j,2)
      end do
      print *

      call writeenergy3_formatted(energyfile,Eword,comboTypes,
     &  contactdir)

499   continue

999   end
