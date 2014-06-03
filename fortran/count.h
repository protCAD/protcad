
      program contact9_0

c Author: C. Summa
c 3/1/98- Passes significance data into bubblesort algortithm. Only those
c 	  clusters which are significant get sorted into TwentyBest/TwentyWorst
c 6/19/99-
c	Tests for existence of database and testset and acts accordingly
c 
c-------------------------------------------------------
c Variable declarations
c        1         2         3         4         5         6         7
c2345678901234567890123456789012345678901234567890123456789012345678901234567890

	

      include 'parameters.h'
      include 'pdb.h'
      include 'pdbselect.h'

      real Distance
      real tempdistance
      real R(maxnumatomtypes,maxnumatomtypes)
      real criticalmaxdistance
      real criticalmindistance 
      real realPairwiseEvents(maxnumatomtypes,maxnumatomtypes)
      real tempArray1(maxnumatomtypes)
      real tempArray2(maxnumatomtypes)
      real fraction(maxnumatomtypes)
      real fExpected(maxnumatomtypes,maxnumatomtypes)
      real fPairwise(maxnumatomtypes,maxnumatomtypes)
      real distance1
      real distance2
      real pTypeCoord(maxnumatomtypes,maxNumCoord)
      real pCoordEnv(maxnumatomtypes,maxNumCoord,570)
      real proteinPTypeCoord(maxat),proteinPCoordEnv(maxat)
      real probCatCoord(maxnumatomtypes,maxNumCoord,
     &       maxnumatomtypesEnv)
      real PA,PB,PC,PD
      real TheoreticalProbWordGCoord(maxnumatomtypes,maxNumCoord,570)
      real TheoreticalProbWord(maxnumatomtypes,maxNumCoord,570)
      real TheoreticalProbWordAlternate(maxnumatomtypes,maxNumCoord,570)
      real TheoreticalProbCoord(maxnumatomtypes,maxNumCoord)
      real TheoreticalProbCoordAlternate(maxnumatomtypes,maxNumCoord)
      real ActualProbWordGCoord(maxnumatomtypes,maxNumCoord,570)
      real ActualProbWord(maxnumatomtypes,maxNumCoord,570)
      real ActualProbCoord(maxnumatomtypes,maxNumCoord)
      real Eword(maxnumatomtypes,maxNumCoord,570)
      real EwordAlternate(maxnumatomtypes,maxNumCoord,570)
      real energyPerResidue(maxNumRes)
      real energyPerResidueAlternate(maxNumRes)
      real energyPerResiduePairwise(maxNumRes)
      real energyPerResidueHistogram(maxNumPDBs,50)
      real energyPerResidueAlternateHistogram(maxNumPDBs,50)
      real energyPerResidueMean(maxNumPDBs)
      real energyPerResidueAlternateMean(maxNumPDBs)
      real energyPerAtomMean(maxNumPDBs)
      real energyPerAtomAlternateMean(maxNumPDBs)
      real TotalEnergy(maxNumPDBs)
      real TotalEnergyAlternate(maxNumPDBs)
      real TotalEnergyPairwise(maxNumPDBs)
      real TotalEnergySingleBody(maxNumPDBs)
      real ETcomb(maxNumPDBs)
      real realEvent(maxnumatomtypes,maxNumCoord,570)
      real sumEvents(maxnumatomtypes,maxNumCoord)
      real totalSumAtomsCatCoord
      real ExpectedEvents(maxnumatomtypes,maxNumCoord,570)
      real ExpectedEventsAlternate(maxnumatomtypes,maxNumCoord,570)
      real ChiSquared(maxnumatomtypes,maxNumCoord,570)
      real ChiSquaredAlternate(maxnumatomtypes,maxNumCoord,570)
      real ChiSqrTable(13,30)
      real ChiSqrProb(13)
      real Significance(maxnumatomtypes,maxNumCoord,570)
      real SignificanceAlternate(maxnumatomtypes,maxNumCoord,570)
      real SumSignificanceHis 
      real SigGT90 
      real PerSigGT90 
      real PerSigGT90Alternate
      real SumUnusableEvents
      real SumTotalEvents
      real NumberSignificantlyPopulated
      real NumberNotSignificantlyPopulatedAlternate
      real NumberSignificantlyPopulatedAlternate
      real NumberNotSignificantlyPopulated
      real sumAtomsCatCoord(maxnumatomtypes,maxNumCoord,
     &       maxnumatomtypesEnv)
      real event(maxnumatomtypes,maxNumCoord,570)
      real alteredEvent(maxnumatomtypes,maxNumCoord,570)
      real pairwiseEvents(maxnumatomtypes,maxnumatomtypes)
      real totalInteractions
      real Poisson
      real TotalSignificantEvents
      real TotalNonSignificantEvents
      real PercentageSignificantEvents
      real rms
      real singleBodyEnergy(maxnumatomtypes,binsInHistogram)
      real normConst
      real minenergy,maxenergy,energyStepsize
      real CoordHist(maxnumatomtypes,binsInHistogram)
      real CoordHistAverage(binsInHistogram)
      real CoordHistTop(maxnumatomtypesTop,binsInHistogram)
      real CoordHistAverageNorm(binsInHistogram)
      integer numBinsInHistogram
      real binsInHistogramReal
      real overallSumEvents
      real version
      real TotalTheoreticalProbWord
      real numatomsInCategory(maxnumatomtypes,3)
      real numatomsInCategoryTotal(maxnumatomtypes,3)
      real tempreal
      real Undetermined(maxNumPDBs)
      real numatmReal
      real TwentyBestE(maxnumatomtypes,20)
      real TwentyBestAltE(maxnumatomtypes,20)
      real TwentyWorstE(maxnumatomtypes,20) 
      real TwentyWorstAltE(maxnumatomtypes,20)
      real energyPerAtom(maxat)
      real energyPerAtomAlt
      real numenvTotal
      real Emean
      real EmeanA
      real EmeanSB
      real EmeanP
      real EmeanT
      real numpdbsreal
      real sigma
      real sigmaA
      real sigmaSB
      real sigmaP
      real sigmaT
      real Z
      real ZA
      real ZSB
      real ZP
      real ZT
      real pmfShellCounter(maxnumatomtypes,maxnumatomtypesEnv,
     &       pmfMaxNumBins)
      real pmfMaxdis
      real pmfMindis
      real PMFENERGY(maxnumatomtypes,maxnumatomtypesEnv,pmfMaxNumBins)
      real pmfAvg(maxnumatomtypes,maxnumatomtypesEnv)
      real pmfTestenergy(maxNumPDBs,maxat2)
      real pmfTesttotal(maxNumPDBs)
      real numBad(maxNumPDBs)
      real numGood(maxNumPDBs)
      real pmfAvgMaxshellReal
      real randm
      real ran1
      real*8 populationValency
      real*8 numatmTotal2Real
      real*8 probSum(maxat)
      double precision Fact

      real    percentToSample
      integer numToSample
      integer subset(maxat)
      integer iseed

      real    TESubset(maxNumPDBs)
      real    TEASubset(maxNumPDBs)
      real    pmfSubset(maxNumPDBs)
      integer sscount
      logical foundOne

      real ThingsinEnv(maxnumatomtypes,maxNumCoord,maxnumatomtypesEnv)
      real SumThingsinEnv(maxnumatomtypes,maxNumCoord)
      real numInCatCoord(maxnumatomtypes,maxNumCoord)

      real    pair(maxnumatomtypes,maxnumatomtypesEnv)
      real    pairAvg(maxnumatomtypesEnv)
      real    N1(maxnumatomtypes)
      real    p1(maxnumatomtypes,maxnumatomtypesEnv)
      real    TEnergyPair(maxNumPDBs)
      real    EPair(maxnumatomtypes,maxnumatomtypesEnv)

      real A1
      real A2
      real A3

      integer  null
      integer  category(maxat,3)
      integer  i1,i2,i3,i4,i5,i,j,k
      integer  pdbsread
      integer  numarguments
      integer  numskip
      integer  Coordination(maxat)
      integer  tempInteger
      integer  numatomtypes
      integer  numatomtypesEnv
      integer  numatomtypesTop
      integer  numBackboneToCount
      integer  backboneCategory
      integer  backboneCounted
      integer  comboTypes(maxNumCoord)
      integer  comboTypesComp(maxNumCoord,570,maxnumatomtypesEnv)
      integer  atomData(maxat,7)
      integer  criticalCoordination
      integer  probSumHist(10)
      integer  numA,numB,numC,numD
      integer  SignificanceBool(maxnumatomtypes,maxNumCoord,570)
      integer  SignificanceBoolAlternate(maxnumatomtypes,
     &           maxNumCoord,570)
      integer  mapAtomTopology(maxnumatomtypes) 
      integer  whichPdb
      integer  follow1,follow2,follow3
      integer  jason
      integer  TwentyBestEMap(maxnumatomtypes,20,2)
      integer  TwentyBestAltEMap(maxnumatomtypes,20,2)
      integer  TwentyWorstEMap(maxnumatomtypes,20,2)
      integer  TwentyWorstAltEMap(maxnumatomtypes,20,2)
      integer  SigBoolAt(maxat)
      integer  numchainused
      integer  nHlxTot
      integer  nBtaTot
      integer  iocode
      integer  pmfShell
      integer  maxpmfShell
      integer  numdistances
      integer  pmfAvgMaxshell
      integer  StringLength
      integer*8  numatmTotal1
      integer*8  numatmTotal2
      integer*8  Pmult
      integer*8  SignificanceHistogram(14)
      integer*8  SignificanceHistogramAlternate(14)
      character*80  numskipstr,crdiststr,numBackboneStr,iseedstr
      character*80  percentstr
      character*50  selectname
      character*44  databasedir
      character*50  pdbname
      character*45  contactpdbdir
      character*44  contactdir
      character*41  misfolddir
c      character*80  misfolddir
c      character*80  contactdir
      character*30  energyfile
      character*30  categorylistname
      character*30  pdblistname2
      character*20  reason
      character*11  filename
      character*1   perResidueString
      logical  relevant,hetatm
      logical  Abool,Bbool,Cbool,Dbool
      logical  firstTimeThrough
      logical  firsttime
      logical  badAtom
      logical  perResidue
      logical  rejectflag(maxNumPDBs)

      integer templen
      integer tempstart
      integer tempend
      character*120 tstring
      character*120 dbconcat
      character*120 mfconcat
      logical dbexist
      logical mfexist


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
44    format (f8.3,1x,i4,1x,i4,1x,i4,1x,i4,1x,i4,1x,i4,
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
     + f8.3,1x,f8.3,1x,f8.3,1x,f8.3,1x,f8.3,1x,f8.3,1x,
     + f8.3,1x,f8.3,1x,f8.3,1x,f8.3,1x,f8.3,1x,f8.3,1x,
     + f8.3,1x,f8.3,1x,f8.3,1x,f8.3,1x,f8.3,1x,f8.3,1x,
     + f8.3,1x,f8.3,1x,f8.3,1x,f8.3,1x,f8.3,1x,f8.3,1x,
     + f8.3,1x,f8.3,1x,f8.3,1x,f8.3,1x,f8.3,1x,f8.3,1x,
     + f8.3,1x,f8.3,1x,f8.3,1x,f8.3,1x,f8.3,1x,f8.3,1x,
     + f8.3,1x,f8.3,1x,f8.3,1x,f8.3,1x,f8.3,1x,f8.3,1x,
     + f8.3,1x,f8.3,1x,f8.3,1x,f8.3,1x,f8.3,1x,f8.3,1x,
     + f8.3,1x,f8.3,1x,f8.3,1x,f8.3,1x,f8.3,1x,f8.3,1x)
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
     + i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,
     + i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,
     + i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,
     + i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8,1x)
71    format (a11,2x,a6,1x,i2,1x,i2,1x,i2,1x,a3,1x,i5,1x,a1,2x,a4)
72    format (i4,1x,a3,1x,a1,1x,i2,1x,a3,1x,i2,1x,i2,1x,i2,1x,
     + i2,1x,i2,f8.3,1x,i1)
73    format (a150)
74    format (a9,42x,a10,2x,a10,a10,2x,a10,1x,a10,1x,a10)
75    format (a50,2x,f10.3,1x,f10.3,1x,f10.3,
     & 1x,f10.3,1x,f10.3,1x,f10.3)
76    format (a50,2x,f10.3,1x,f10.3,1x,f10.3,1x,f10.3,1x,f10.3,1x)
77    format (i3,1x,i3,1x,f10.7,1x,f10.7,1x,f10.7,1x,f10.7)

c------------------------------------------------------------
c Read in command line arguments

      numarguments = iargc()

      if (numarguments.ne.6) then
        print *,'Usage: con9.0 PDBselectname
     +  #_residues_to_skip critcaldistance 
     + percent_to_sample misfoldfile iseed'
        goto 999
      endif

      call getarg(1,selectname)
      call getarg(2,numskipstr)
      call gettheinteger(numskip,numskipstr)
      call getarg(3,crdiststr)
      call getthereal(criticalmaxdistance,crdiststr)
c      call getarg(4,numBackboneStr)
c      call gettheinteger(numBackboneToCount,numBackboneStr)
      call getarg(4,percentstr)
      call getthereal(percentToSample,percentstr)
      call getarg(5,pdblistname2)
      call getarg(6,iseedstr)
      call gettheinteger(iseed,iseedstr)

      call getenv('CONTACTDIR',contactdir)
      call getenv('CONTACTPDBDIR',contactpdbdir)
      call getenv('CONTACTPDBDIR2',databasedir)
      call getenv('MISFOLDDIR',misfolddir)

      criticalmaxdistance = criticalmaxdistance / 10.0
      version = 9.1
 
      perResidueString = 'n'
      numBackboneStr = '1'
      numBackboneToCount = 1

      energyfile = ''
      energyfile(1:2) = crdiststr(1:2)
      energyfile(3:3) = '.'
      energyfile(4:4) = numskipstr(1:1)
      energyfile(5:5) = '.'
      energyfile(6:6) = numBackboneStr(1:1)
      energyfile(7:13) = '.ene9.0'

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
      print *,'Tested list: ',pdblistname2
      print *,'Number of residues to skip: ',numskip
      print *,'Critical distance : ',criticalmaxdistance
c      print *,'Modulus of backbone residues to count: ',numBackboneToCount
      print *,'Percent sampled in test subset: ',percentToSample
      print *,'iseed = ',iseed
      print *

c     Test for existence of database for derivation:
c  databasedir

      call fileexist(selectname,databasedir,dbconcat,dbexist)
      if (.not.dbexist) then
       print *
       print *,'No database...running in test only mode'
       print *
      end if

c     Test for existence of test set
c  misfolddir

      call fileexist(pdblistname2,misfolddir,mfconcat,mfexist)
      if (.not.mfexist) then
       print *
       print *,'No misfold test set...runnining in derive only mode'
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


      print *,'pmfMaxdis = ',pmfMaxdis
      print *,'pmfMindis = ',pmfMindis
      print *,'maxpmfShell = ',maxpmfShell
      print *,'numdistances = ',numdistances
      print *,'pmfAvgMaxshell = ',pmfAvgMaxshell
      print *,'pmfAvgMaxshellReal = ',pmfAvgMaxshellReal

      categorylistname = 'FullAtomList8.1.dat'
     
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
       do k=1,maxnumatomtypes
        pairwiseEvents(i,k) = 0.0
        fPairwise(i,k) = 0.0
        R(i,k) = 0.0
       end do
       do k=1,numBinsInHistogram
        CoordHist(i,k) = 0.0
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
 
      do i=1,maxnumatomtypesTop
       do k=1,numBinsInHistogram
        CoordHistTop(i,k) = 0.0
       end do
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


      numatmTotal1 = 0
      numatmTotal2 = 0
      numenvTotal = 0.0
      pdbsread = 0
      numchainused = 0
      hetatm = .false.
      totalInteraction = 0
      backboneCategory = 1
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

      open (unit = 21,file=contactdir//energyfile,
     &  status = 'old',iostat = iocode,err = 90)
      close (unit = 21, err = 90)

90    if (iocode.eq.0) then
       print *,'File containing derived energies exists....'
       print *
       print *,'Skipping derivation....'
       print *
       if (.not.mfexist) then
        print *,'No misfolds to test...Nothing to do!!!'
        print *,'Stopping'
        stop
       end if
       if (.not.dbexist) then 
        print *,'No energyfile, no database, no way of calculating
     & energies.'
        print *,'Stopping'
        stop
       end if
       call readenergy(energyfile,EwordAlternate,Eword,singleBodyEnergy,
     &    comboTypes,fPairwise,PMFENERGY,EPair,contactdir)

       print *,'Energies read from file'
       goto 499
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

c       print *,'numatm = ',numatm
       call pdbread_w_write (filename,databasedir)
c       print *,'numatm = ',numatm

c       print *,'read pdbfile' 
  
       badAtom = .false. 
       firsttime = .true.

       call pdbcheck (contactdir,filename,reason,firsttime)

c        do i=1,numatm
c        print *,atheader(i),atomnum(i),atomname(i),
c     &   multiconf(i),res(i),chain(i),
c     &   trueResnum(i),subnum(i)
c        end do

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
 
c      print *,'categorized atoms'

      if (badAtom) then
c       print *,'numatm before collapse',numatm
c       call pdbcollapse(category)
c       print *,'numatm after collapse',numatm
c       do i=1,numatm
c        print *,atheader(i),atomnum(i),atomname(i),
c     &   multiconf(i),res(i),chain(i),
c     &   trueResnum(i),subnum(i)
c        end do   
c       print *,'collapsed pdb'
      end if

c      call pdbwrite (filename,databasedir,cmnt,atheader,
c     &   atomnum,atomname,res,chain,trueResnum,subnum,coord,
c     &   numatm,numres,connectivity,numchains,numcomlines,
c     &   Bfac,rad)


c      print *,'Begin printing unknown atoms'
c      do i =1,numatm
c       if (unknown(i)) then
c        print *,i,atheader(i),atomname(i),res(i)
c       end if
c      end do
c      print *,'End printing unknown atoms'

c      do i=1,numatm
c       print *,i,category(i,1),category(i,2)
c      end do

       call findRelevantChain (pdbChain(ii))

c       print *,'done with find Relevant Chain'
c       print *,'outside:chainStart :',chainStart
c       print *,'outside:chainEnd :',chainEnd

c      do i=chainStart,chainEnd
c       print *,atomname(i),res(i),trueResnum(i),
c     &  category(i,1),unknown(i)
c      end do 

c      print *,'entering second iteration of pdbcheck'

       badAtom = .false.
  
       call pdbcheck (contactdir,filename,reason,firsttime)

c      print *,'done with second check of pdbdata'   

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
       atomData(i1,1) = category(i1,1)
       backboneCounted = 0

       do i2 = 1,numatm

        call Relevant_Distance(i1,i2,numskip,res(i1),
     &     res(i2),atomname(i1),atomname(i2),
     &     chain(i1),chain(i2),resnum,relevant)
        if (relevant) then
         tempdistance = Distance(coord(1,i1),coord(2,i1),
     &     coord(3,i1),coord(1,i2),
     &     coord(2,i2),coord(3,i2))
         call Make_Bin_Decision(tempdistance,pmfMaxdis,pmfMindis,
     &      pmfStep,pmfShell)
c!!!
         if (pmfShell.eq.1) then
          print *,'first_shell!  ',filename,i1,i2
         end if
c!!!
         if (pmfShell.ne.999) then
           pmfShellCounter(category(i1,1),
     &                      category(i2,2),pmfShell) =
     &     pmfShellCounter(category(i1,1),
     &                      category(i2,2),pmfShell) + 1.0
         end if
         if (tempdistance.lt.criticalmaxdistance.and.
     &       tempdistance.gt.criticalmindistance) then
          if (category(i2,2).eq.backboneCategory) then
           if (mod(backboneCounted,numBackboneToCount).eq.0) then
            Coordination(i1) = Coordination(i1) + 1
            atomData(i1,2) = atomData(i1,2) + 1
            atomData(i1,(category(i2,2)+2)) = 
     &        atomData(i1,(category(i2,2)+2)) + 1
           end if
          backboneCounted = backboneCounted + 1
          else
           Coordination(i1) = Coordination(i1) + 1
           atomData(i1,2) = atomData(i1,2) + 1
           atomData(i1,(category(i2,2)+2)) =
     &       atomData(i1,(category(i2,2)+2)) + 1
          end if
         end if
        end if
       end do
      end do


c- repeat counting for pairwise calculations.....
c- change range of loop variables (i2, specifically) to
c- avoid double counting

      do i1 = chainStart,chainEnd
       do i2 = 1,numatm
        call Relevant_Distance(i1,i2,numskip,res(i1),
     &     res(i2),atomname(i1),atomname(i2),
     &     chain(i1),chain(i2),resnum,relevant)
        if (relevant) then
         tempdistance = Distance(coord(1,i1),coord(2,i1),
     &     coord(3,i1),coord(1,i2),
     &     coord(2,i2),coord(3,i2))
         if (tempdistance.lt.criticalmaxdistance.and.
     &       tempdistance.gt.criticalmindistance) then
          pairwiseEvents(category(i1,1),category(i2,1)) =
     &       pairwiseEvents(category(i1,1),category(i2,1)) + 1
          pair(category(i1,1),category(i2,2)) = 
     &       pair(category(i1,1),category(i2,2)) + 1
         end if
        end if
       end do
      end do

      numatmTotal1 = numatmTotal1 + numatm

c- Calculating the number of environment atoms in the sphere
c- of influence surrounding each atomtype/coordination number
c- combination....

      do i1 = chainStart,chainEnd
       if (atomdata(i1,2).ge.15) then
        atomdata(i1,2) = 14
       end if
       do i2 = 1,4
        ThingsinEnv(atomdata(i1,1),atomdata(i1,2)+1,i2) =
     &   ThingsinEnv(atomdata(i1,1),atomdata(i1,2)+1,i2) +
     &   atomdata(i1,i2+2)
       end do
       numInCatCoord(atomdata(i1,1),atomdata(i1,2)+1) =
     &  numInCatCoord(atomdata(i1,1),atomdata(i1,2)+1) +
     &  1.00 
      end do


c Now we'll group the atoms based on atomtype(1-20),coordination(0-13),
c and the #s of atoms in each of the four environment variables.

      do i1 = chainStart,chainEnd
       if (atomData(i1,2).le.criticalCoordination) then
        do i2 = 1,comboTypes(atomData(i1,2)+1)
         if (atomData(i1,3).eq.
     &      comboTypesComp(atomData(i1,2)+1,i2,1)) then
          if (atomData(i1,4).eq.
     &       comboTypesComp(atomData(i1,2)+1,i2,2)) then
           if (atomData(i1,5).eq.
     &        comboTypesComp(atomData(i1,2)+1,i2,3)) then
            if (atomData(i1,6).eq.
     &         comboTypesComp(atomData(i1,2)+1,i2,4)) then
              event(atomData(i1,1),atomData(i1,2)+1,i2) = 
     &         event(atomData(i1,1),atomData(i1,2)+1,i2) + 1
              alteredEvent(atomData(i1,1),atomData(i1,2)+1,
     &         i2) = alteredEvent(atomData(i1,1),
     &         atomData(i1,2)+1,i2) + 1
             continue
            end if
           end if
          end if
         end if
        end do
       else 
c      - Approximation to lower coordination composition 
c        may eventually be inserted here
c       print *,pdbname,i1,atomData(i1,2)
        event(atomData(i1,1),15,1) = 
     +    event(atomData(i1,1),15,1) + 1
        alteredEvent(atomData(i1,1),15,1) =
     +    alteredEvent(atomData(i1,1),15,1) + 1
       end if
      end do

c---------------------------------------------------------------
c Calculation of Histograms for Coordination Numbers
c based on atom type, and based on topological type

      do i1 = chainStart, chainEnd
       if (Coordination(i1).ge.0.and.Coordination(i1)
     +   .le.(numBinsInHistogram-2)) then
          CoordHist(category(i1,1),Coordination(i1)+1) = 
     +     CoordHist(category(i1,1),Coordination(i1)+1) + 1
          CoordHistTop(category(i1,3),Coordination(i1)+1) =
     +     CoordHistTop(category(i1,3),Coordination(i1)+1) + 1
       else
        CoordHist(category(i1,1),numBinsInHistogram) = 
     +   CoordHist(category(i1,1),numBinsInHistogram) + 1
        CoordHistTop(category(i1,3),numBinsInHistogram) =
     +   CoordHistTop(category(i1,3),numBinsInHistogram) + 1
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
      
c!!!!
c      goto 999
c!!!!
cEVERYTHING from here up can be put into a subroutine


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
       print *,'atomtype = ',i1
       rms = 0.0
       normConst = numatomsInCategoryTotal(i1,1) / numatmTotal2Real
c       print *,'normConst = ',normConst
       do i2 = 1,maxNumCoord
        CoordHistAverageNorm(i2) = CoordHistAverage(i2) *
     +    normConst
        tempreal = CoordHist(i1,i2)/CoordHistAverageNorm(i2)
c        print *,i1,i2,CoordHist(i1,i2),
c     +   CoordHistAverageNorm(i2),tempreal,-log(tempreal)
        if ((CoordHist(i1,i2).ne.0.0).and.
     +   (CoordHistAverageNorm(i2).ne.0.0)) then
          singleBodyEnergy(i1,i2) = -RT * 
     +      log(CoordHist(i1,i2)/
     +     CoordHistAverageNorm(i2))
        end if
        if (CoordHist(i1,i2).eq.0.0) then
         singleBodyEnergy(i1,i2) = -RT * 
     +     log(1/CoordHistAverageNorm(i2))
        end if
        if (CoordHistAverageNorm(i2).eq.0.0) then
          singleBodyEnergy(i1,i2) = 0.0
        end if
c        print *,'SingleBodyEnergy =', singleBodyEnergy(i1,i2)
        rms = rms + (CoordHistAverageNorm(i2) -
     +     CoordHist(i1,i2))**2
       end do
       rms = sqrt(rms/maxNumCoord)
c       print *,'rms for atom type ',i1,' is ',rms
      end do

      print *,'Coordination # Histogram (Atomtype)'
      print *,'X axis is atomtype, Y axis is Coordination #'
      tempInteger = numBinsInHistogram
      print 29,(j,j=0,numatomtypes)
      do i1 = 1,tempInteger
       print 18,i1-1,CoordHistAverage(i1),
     +   (CoordHist(j,i1),j=1,numatomtypes)
      end do

c      print *
c      print *
c      print *,'Coordination # Histogram (Topology)'
c      print *
c      tempInteger = numBinsInHistogram
c      print 29,(j,j=0,tempInteger-1)
c      print *
c      do i1 = 1,numatomtypesTop
c       print 18,i1,(CoordHistTop(i1,j),j=1,numBinsInHistogram)
c      end do

c----------------------------------------------------------------
c Calculation of the Theoretical probability of seeing a particular
c coordination number about an atom 'topological type' along with
c the actual probability based on atom type.

      do i2 = 1,maxNumCoord
        do i1 = 1,numatomtypesTop
        if (numatomsInCategoryTotal(i1,3).ne.0.0) then
         TheoreticalProbCoord(i1,i2) = 
     &    CoordHistTop(i1,i2) /  
     &    numatomsInCategoryTotal(i1,3)
        else
          TheoreticalProbCoord(i1,i2) = 0.0
        end if
        end do
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


      print *,'ActualProbcoord(5,1) = ',ActualProbCoord(5,1)

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




c OLD VERSION of calculation of probCatCoord
c Calculate the _actual_ probability of drawing an A,B,C,or D 
c given a random distribution
c Do this on an atomtype, coordination # basis

c      totalSumAtomsCatCoord = 0.0
c      do i1 = 1,numatomtypes
c       do i2 = 1,maxNumCoord
c        do i3 = 1,4
c          probCatCoord(i1,i2,i3) = 
c     +     numatomsInCategoryTotal(i3,2)/
c     +     numenvTotal
c        end do
c       end do
c      end do

c      print *,'#Specific Type / Total# -> probCatCoord(n,x,y)'
c      print *,'      A  ','    B  ','    C  ','    D  '
c      do i=1,maxNumCoord
c       print 19,i,(probCatCoord(follow1,i,j),j=1,4)
c      end do

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
        PA = p1(i1,1)**numA
        PB = p1(i1,2)**numB
        PC = p1(i1,3)**numC
        PD = p1(i1,4)**numD

        if (Pmult.eq.0) then
         TheoreticalProbWordGCoord(i1,i2,i3) = 0.0
        end if 

c        if (i1.eq.follow1.and.i2.eq.follow2.
c     &      and.i3.eq.follow3) then
c         print *
c         print *,'Raw theoretical data for atom types "follow"'
c         print *,'     numA    numB    numC    numD     Pmult    
c     &PA    PB     PC     PD'
c         print 11,numA,numB,numC,numD,Pmult,PA,PB,PC,PD
c         print *
c        end if

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
c seeing a given clustter (or word) given an atomtype and
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
        TheoreticalProbCoordAlternate(i1,i2) = 
     +   CoordHistAverage(i2) / numatmTotal2Real
       end do
      end do

c*************************************************
c OLD VERSION
c Calculate ActualProbWord,TheoreticalProbWord,TheoreticalProbWordAlternate

c      do i1 = 1,numatomtypes
c       do i2 = 1,maxNumCoord
c        do i3 = 1,comboTypes(i2)
c         ActualProbWord(i1,i2,i3) =
c     +      ActualProbWordGCoord(i1,i2,i3) *
c     +      ActualProbCoord(i1,i2)
c         TheoreticalProbWord(i1,i2,i3) =
c     +      TheoreticalProbWordGCoord(i1,i2,i3) *
c     +      TheoreticalProbCoord(mapAtomTopology(i1),i2)
c         TheoreticalProbWordAlternate(i1,i2,i3) =
c     +      TheoreticalProbWordGCoord(i1,i2,i3) *
c     +      TheoreticalProbCoordAlternate(i1,i2)
c        end do
c       end do
c      end do

c*************************************************
c NEW VERSION
c Calculate ActualProbWord,TheoreticalProbWord,TheoreticalProbWordAlternate

      do i1 = 1,numatomtypes
       do i2 = 1,maxNumCoord
        do i3 = 1,comboTypes(i2)
         ActualProbWord(i1,i2,i3) =
     +      ActualProbWordGCoord(i1,i2,i3)
         TheoreticalProbWord(i1,i2,i3) =
     +      TheoreticalProbWordGCoord(i1,i2,i3)
         TheoreticalProbWordAlternate(i1,i2,i3) =
     +      TheoreticalProbWordGCoord(i1,i2,i3)
        end do
       end do
      end do


c**************************************************
c OLD VERSION
c Calculate ExpectedEvents assuming random sampling

c      do i1 = 1,numatomtypes
c       do i2 = 1,maxNumCoord
c        do i3 = 1,comboTypes(i2)
c         ExpectedEvents(i1,i2,i3) =
c     +     TheoreticalProbWord(i1,i2,i3) *
c     +     numatomsInCategoryTotal(i1,1)
c         ExpectedEventsAlternate(i1,i2,i3) =
c     +     TheoreticalProbWordAlternate(i1,i2,i3) *
c     +     numatomsInCategoryTotal(i1,1)
c         TotalTheoreticalProbWord = TotalTheoreticalProbWord +
c     +     TheoreticalProbWord(i1,i2,i3)
c        end do
c       end do
c      end do


c**************************************************
c NEW VERSION
c Calculate ExpectedEvents assuming random sampling

      do i1 = 1,numatomtypes
       do i2 = 1,maxNumCoord
        do i3 = 1,comboTypes(i2)
         ExpectedEvents(i1,i2,i3) =
     +     TheoreticalProbWord(i1,i2,i3) *
     +     numInCatCoord(i1,i2)
         ExpectedEventsAlternate(i1,i2,i3) =
     +     TheoreticalProbWordAlternate(i1,i2,i3) *
     +     numInCatCoord(i1,i2)
c         TotalTheoreticalProbWord = TotalTheoreticalProbWord +
c     +     TheoreticalProbWord(i1,i2,i3)
        end do
       end do
      end do



c      print *,'TotalTheoreticalProbWord = ',TotalTheoreticalProbWord

c-------------------------------------------------------------------
c  Significance Calculation
c-------------------------------------------------------------------


c Calculate ChiSquared statistic for each atomtype-coordination
c as deviation from random for all cases

      call CalculateSignificance(numatomtypes,comboTypes,ChiSqrTable,
     +  ChiSqrProb,ExpectedEvents,event,ChiSquared,Significance,
     +  SignificanceBool,SignificanceHistogram,PerSigGT90)

      call CalculateSignificance(numatomtypes,comboTypes,ChiSqrTable,
     +  ChiSqrProb,ExpectedEventsAlternate,event,ChiSquaredAlternate,
     +  SignificanceAlternate,SignificanceBoolAlternate,
     +  SignificanceHistogramAlternate,PerSigGT90Alternate)

      do i1 = 1,numatomtypes
       do i2 = 1,15
        do i3 = 1,combotypes(i2)
         if (SignificanceBool(i1,i2,i3).eq.1) then
          TotalSignificantEvents = TotalSignificantEvents + 
     +      event(i1,i2,i3)
         else
          TotalNonSignificantEvents = TotalNonSignificantEvents + 
     +      event(i1,i2,i3)
         end if
         if (SignificanceBoolAlternate(i1,i2,i3).eq.1) then
          TotalSignificantEventsAlternate = 
     &      TotalSignificantEventsAlternate +
     &      event(i1,i2,i3)
         else
          TotalNonSignificantEventsAlternate = 
     &     TotalNonSignificantEventsAlternate +
     &     event(i1,i2,i3)  
         end if
        end do
       end do
      end do

      PercentageSignificantEvents = (TotalSignificantEvents/
     + (TotalSignificantEvents + TotalNonSignificantEvents))*100
      PercentageSignificantEventsAlternate = 
     & (TotalSignificantEventsAlternate/
     & (TotalSignificantEventsAlternate + 
     & TotalNonSignificantEventsAlternate))*100



c Print out the chisq data

c      print *,'Chi Squared (ExpectedEvents) data for atomtype ',follow1
c      print *
c      print 29,(i1,i1=1,numatomtypes)
c      do i2 = 1,jason
c       print 35,i2,(ChiSquared(follow1,i2,i3),
c     +   i3 = 1,comboTypes(i2))
c      end do
c      print *,'Chi Squared (ExpectedEventsAlternate) data for atomtype '
c     & ,follow1
c      print *
c      do i2 = 1,jason
c       print 35,i2,(ChiSquaredAlternate(follow1,i2,i3),
c     +   i3 = 1,comboTypes(i2))
c      end do

c------------------------------------------------
c  Print Significance data for atom type 18

c       print *
c       print *,'Significance of each data point for atom type ',follow1
c       do i2 = 1,jason
c        print 68,(Significance(follow1,i2,i3),
c     +   i3 = 1,comboTypes(i2))
c       end do
c       do i2 = 1,jason
c        print 68,(SignificanceAlternate(follow1,i2,i3),
c     +   i3 = 1,comboTypes(i2))
c       end do

c-------------------------------------------------
c Print SignificanceBool data for atom type follow1

c       print *
c       print *,'SignificanceBool of each data point for atom type ',
c     &  follow1
c       do i2 = 1,jason
c	   print 70,(SignificanceBool(follow1,i2,i3),
c     +   i3 = 1,comboTypes(i2))
c       end do
c       do i2 = 1,jason
c        print 70,(SignificanceBoolAlternate(follow1,i2,i3),
c     +   i3 = 1,comboTypes(i2))
c       end do

c-------------------------------------------------
       print *,'Significance Data for Original Energy Calculation'
       print *,'Significance Histogram over full data set'
       print 29,(SignificanceHistogram(i1),i1 = 1,14)

       print *,'PerSigGT90 = ',PerSigGT90
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

       print *,'Significance Data for Alternate Energy Calculation'
       print *,'Significance Histogram for Alternate Energy'
       print 29,(SignificanceHistogramAlternate(i1),i1 = 1,14)

       print *
       print *
       print *,'TotalSignificantEventsAlternate',
     &  TotalSignificantEventsAlternate
       print *,'TotalNonSignficantEventsAlternate',
     +    TotalNonSignificantEventsAlternate
       print *,'PercentageSignficantEventsAlternate',
     +    PercentageSignificantEventsAlternate


c------------------------------------------------
c Print out some sample 'event' data

c      print *,'Event data for atom type ',follow1
c      do i=1,jason
c       print 68,(event(follow1,i,k),k=1,comboTypes(i))
c       print *
c      end do


c**********************************************
c if events=0, then add 1 to events for upcoming
c energy calculation

c The value 3.56 is the critical value of the mean
c given a Poisson value of 1, where the significance is 
c greater than 90% (Assures that our 0 to 1 approximation
c does not alter the significance of the counting statistics

       do i1 = 1,numatomtypes
        do i2 = 1,maxNumCoord
         do i3 = 1,comboTypes(i2)
          if (event(i1,i2,i3).eq.0) then
            alteredEvent(i1,i2,i3) = 
     +      alteredEvent(i1,i2,i3) + 1
          end if
          if (ExpectedEvents(i1,i2,i3).lt.3.56) then
           SignificanceBool(i1,i2,i3) = 0
          end if
          if (ExpectedEventsAlternate(i1,i2,i3)
     &     .lt.3.56) then
           SignificanceBoolAlternate(i1,i2,i3) = 0
          end if
         end do
        end do
       end do

c***********************************************
c Energy Calculation


      do i1 = 1,numatomtypes
       do i2 = 1,maxNumCoord
        do i3 = 1,comboTypes(i2)
         if (SignificanceBool(i1,i2,i3).eq.1) then
          Eword(i1,i2,i3) = -RT * 
     +     log(alteredEvent(i1,i2,i3)/
     +     ExpectedEvents(i1,i2,i3))
         else
          Eword(i1,i2,i3) = singleBodyEnergy(i1,i2)
         end if
         if (SignificanceBoolAlternate(i1,i2,i3).eq.1) then
          EwordAlternate(i1,i2,i3) = - RT * 
     +     log(alteredEvent(i1,i2,i3)/
     +     ExpectedEventsAlternate(i1,i2,i3))
         else
          EwordAlternate(i1,i2,i3) = 
     +     singleBodyEnergy(i1,i2)
         end if
        end do
       end do
      end do 

c-------------------------------------------------
c TwentyBest,TwentyWorst

c!!!  goto 320

      call BubbleSort(Eword,TwentyBestE,TwentyBestEMap,TwentyWorstE,
     & TwentyWorstEMap,numatomtypes,comboTypes,SignificanceBool)

      call BubbleSort(EwordAlternate,TwentyBestAltE,TwentyBestAltEMap,
     & TwentyWorstAltE,TwentyWorstAltEMap,numatomtypes,comboTypes,
     & SignificanceBoolAlternate)

      print *,'Contents of TwentyBestAltEMap'
      do i=1,20
       do j=1,20
        print *,TwentyBestAltEMap(i,j,1),TwentyBestAltEMap(i,j,2)
       end do
      end do


      do i = 1,20
       print *
       print *,'Twenty Best for atom type ',i
       do j=1,20
        print 44,TwentyBestAltE(i,j),
     &   TwentyBestAltEMap(i,j,1)-1,
     &   TwentyBestAltEMap(i,j,2),
     &   comboTypesComp(TwentyBestAltEMap(i,j,1),
     &   TwentyBestAltEMap(i,j,2),1),
     &   comboTypesComp(TwentyBestAltEMap(i,j,1),
     &   TwentyBestAltEMap(i,j,2),2),
     &   comboTypesComp(TwentyBestAltEMap(i,j,1),
     &   TwentyBestAltEMap(i,j,2),3),
     &   comboTypesComp(TwentyBestAltEMap(i,j,1),
     &   TwentyBestAltEMap(i,j,2),4),
     &   event(i,TwentyBestAltEMap(i,j,1),
     &   TwentyBestAltEMap(i,j,2)),
     &   ExpectedEventsAlternate(i,
     &   TwentyBestAltEMap(i,j,1),
     &   TwentyBestAltEMap(i,j,2))
       end do
       print *
       print *,'Twenty Worst for atom type ',i
       do j=1,20
        print 44,TwentyWorstAltE(i,j),
     &   TwentyWorstAltEMap(i,j,1)-1,
     &   TwentyWorstAltEMap(i,j,2),
     &   comboTypesComp(TwentyWorstAltEMap(i,j,1),
     &   TwentyWorstAltEMap(i,j,2),1),
     &   comboTypesComp(TwentyWorstAltEMap(i,j,1),
     &   TwentyWorstAltEMap(i,j,2),2),
     &   comboTypesComp(TwentyWorstAltEMap(i,j,1),
     &   TwentyWorstAltEMap(i,j,2),3),
     &   comboTypesComp(TwentyWorstAltEMap(i,j,1),
     &   TwentyWorstAltEMap(i,j,2),4),
     &   event(j,TwentyWorstAltEMap(i,j,1),
     &   TwentyWorstAltEMap(i,j,2)),
     &   ExpectedEventsAlternate(i,
     &   TwentyWorstAltEMap(i,j,1),
     &   TwentyWorstAltEMap(i,j,2))
       end do
      end do


c!!!!
c320    goto 334
      
c------------------------------------------------------- 
c Here, we're going to go through the protein database again to pull out
c examples of the twenty best and twenty worst.  Unfortunately, this involves
c an almost total repetition of the first section of code for reading in the
c PDB files and categorizing the atoms, etc. (A cut and paste job!)

c-------------------------------------------------------

      print *,'Examples of 20 best and 20 worst'

      do ii=1,pdbsread 
       filename(1:11) = pdbnames(ii)(1:11)


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

       call pdbread (filename,databasedir)


       print *,'reading pdbfile: ',filename

       badAtom = .false.

       call pdbcheck (contactdir,filename,reason,firsttime)

c       print *,'checked pdbfile'
 
       if (reject) then
        print *,'Rejected : ',filename,'  ',reason
        rejectflag(ii) = .true.
        goto 333
       else
        print *,'Read : ',filename
       end if

       call categorizeatoms (numatomtypesEnv,numatomtypesTop,
     &   categorylistname,contactdir,mapAtomTopology,firstTimeThrough,
     &   badAtom,category)

c      print *,'exit categorizeatoms'

      if (badAtom) then

       call pdbcollapse (category)

      end if

       call findRelevantChain (pdbChain(ii))


c-----------------------------------------------------------------
c Now calculate distances, relevancy of pairwise contact
c (based in part on connectivity and atom type)

      do i1 = chainStart,chainEnd
       atomData(i1,1) = category(i1,1)
       backboneCounted = 0

       do i2 = 1,numatm

        call Relevant_Distance(i1,i2,numskip,res(i1),
     &     res(i2),atomname(i1),atomname(i2),
     &     chain(i1),chain(i2),resnum,relevant)
        if (relevant) then
         tempdistance = Distance(coord(1,i1),coord(2,i1),
     &     coord(3,i1),coord(1,i2),
     &     coord(2,i2),coord(3,i2))

         if (tempdistance.lt.criticalmaxdistance.and.
     &       tempdistance.gt.criticalmindistance) then
          if (category(i2,2).eq.backboneCategory) then
           if (mod(backboneCounted,numBackboneToCount).eq.0) then
            Coordination(i1) = Coordination(i1) + 1
            atomData(i1,2) = atomData(i1,2) + 1
            atomData(i1,(category(i2,2)+2)) =
     &        atomData(i1,(category(i2,2)+2)) + 1
           end if
          backboneCounted = backboneCounted + 1
          else
           Coordination(i1) = Coordination(i1) + 1
           atomData(i1,2) = atomData(i1,2) + 1
           atomData(i1,(category(i2,2)+2)) =
     &       atomData(i1,(category(i2,2)+2)) + 1
          end if
         end if
        end if
       end do
      end do

c---------
c Does this file contain any of the best (or worst) environments???
c If so, print out the atom number, the atom name, the residue name,
c and the Best/Worst map in question....


      print *,'pdbname = ',filename
      do i1 = chainStart,chainEnd
       do i2 = 1,20
        i3 = TwentyBestAltEMap(atomData(i1,1),i2,1)
        i4 = TwentyBestAltEMap(atomData(i1,1),i2,2)
        if ((i3-1).eq.Coordination(i1)) then
         if (comboTypesComp(i3,i4,1).eq.atomData(i1,3)) then
          if (comboTypesComp(i3,i4,2).eq.atomData(i1,4)) then
           if (comboTypesComp(i3,i4,3).eq.atomData(i1,5)) then
            if (comboTypesComp(i3,i4,4).eq.atomData(i1,6)) then
             print 71,filename,'Best',category(i1,1),i2,
     &        Coordination(i1),res(i1),trueResnum(i1),chain(i1),
     &        atomname(i1)
            end if
           end if
          end if
         end if
        end if
       end do
       do i2 = 1,20
        i3 = TwentyWorstAltEMap(atomData(i1,1),i2,1)
        i4 = TwentyWorstAltEMap(atomData(i1,1),i2,2)
        if ((i3-1).eq.Coordination(i1)) then
         if (comboTypesComp(i3,i4,1).eq.atomData(i1,3)) then
          if (comboTypesComp(i3,i4,2).eq.atomData(i1,4)) then
           if (comboTypesComp(i3,i4,3).eq.atomData(i1,5)) then
            if (comboTypesComp(i3,i4,4).eq.atomData(i1,6)) then
             print 71,filename,'Worst',category(i1,1),i2,
     &        Coordination(i1),res(i1),trueResnum(i1),chain(i1),
     &        atomname(i1)
            end if
           end if
          end if
         end if
        end if
       end do
      end do


333   end do

334   continue

c--------------------------------------------------
c Add up half of contact matrix to get total number
c of contacts counted

      do i=1,numatomtypes
       do j=i,numatomtypes
        totalInteractions = totalInteractions + 
     +     pairwiseEvents(i,j)
       end do
      end do
 
      print *,'totalInteractions',totalInteractions
      print *,'numatmTotal2',numatmTotal2


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
        print 68,(Significance(i2,i1,k),k=1,combotypes(i1))
       end do
      end do 

c      print *,'ActualProbCoord for atom type ',follow1
c      do i=1,jason
c       print *,ActualProbCoord(follow1,i)
c       print *
c      end do

c      print *,'TheoreticalProbCoord for atom type ',follow1
c      do i=1,jason
c       print 68,TheoreticalProbCoord(follow1,i)
c       print *
c      end do

c      print *,'TheoreticalProbCoordAlternate for atomtype ',follow1
c      do i=1,jason
c       print 68,TheoreticalProbCoordAlternate(follow1,i)
c       print *
c      end do



c      print *,'ActualProbWord for atom type ',follow1
c      do i=1,jason
c       print 68,(ActualProbWord(follow1,i,k),k=1,comboTypes(i))
c       print *
c      end do

c      print *,'TheoreticalProbWord for atom type ',follow1
c      do i=1,jason
c       print 68,(TheoreticalProbWord(follow1,i,k),k=1,comboTypes(i))
c       print *
c      end do

c---------------------------------------------------
c Print out some 'energy' data

c      print *,'Energy data for atom type ',follow1
c      print *,'****************************'
c      do i=1,jason
c       print 68,(Eword(follow1,i,k),k=1,comboTypes(i))
c      end do
c      print *
c      print *,' Alternate Energy'
c      do i=1,jason
c       print 68,(EwordAlternate(follow1,i,k),k=1,comboTypes(i))
c      end do

c---------------------------------------------------
c Calculate the fraction of each atom type

      do i=1,numatomtypes
        fraction(i) = numatomsInCategoryTotal(i,1) /
     +   numatmTotal2Real
      end do

c--------------------------------------------------
c Calculate f_expected(i)

      do i=1,numatomtypes
       do j=1,numatomtypes
        fExpected(i,j) = fraction(i) * fraction(j) *
     &     (totalInteractions/2)
       end do
      end do

c--------------------------------------------------
c Convert the pairwise events array to real values for input
c into correlate subroutine
c 0 events approximation is accomplished here

      do i=1, numatomtypes
       do j=1, numatomtypes
        if (pairwiseEvents(i,j).eq.0) then
         print *,'0 to 1 contact approximation invoked for 
     &contact type ',i,j
         pairwiseEvents(i,j) = 1
        end if
       end do
      end do

c--------------------------------------------------
c Calcuate fPairwise(i,j)

      do i=1,numatomtypes
       do j=1,numatomtypes
        fPairwise(i,j) = -RT * 
     &    log(pairwiseEvents(i,j)/fExpected(i,j))
       end do
      end do

c--------------------------------------------------
c OK, now calculate correlation coefficients for each pair of
c atom types

      do i=1, numatomtypes
       do j=1, numatomtypes
        do k=1,numatomtypes
         tempArray1(k) = fPairwise(i,k)
        end do
        do k=1,numatomtypes
         tempArray2(k) = fPairwise(j,k)
        end do
        call correlate(tempArray1,tempArray2,numatomtypes,
     +   R(i,j))
        do k=1,numatomtypes
         tempArray1(k) = 0.0
        end do
        do k=1,numatomtypes
         tempArray2(k) = 0.0
        end do
       end do
      end do 

c--------------------------------------------------
c Print out some ExpectedEvent data

c      print *,'Expected Event data for atom type ',follow1
c      do i=1,jason
c       print 68,(ExpectedEvents(follow1,i,k),k=1,comboTypes(i))
c       print *
c      end do


c      print *
c      print *
c      print *,'singleBodyEnergy (Atomtype)'
c      print *
c      tempInteger = numBinsInHistogram
c      print 29,(j,j=0,tempInteger-1)
c      print *
c      print *,'numatomtypes = ',numatomtypes
c      print *,'numBinsInHistogram = ',numBinsInHistogram
c      do i1 = 1,numatomtypes
c       print 15,i1,(singleBodyEnergy(i1,j),
c     &   j=1,numBinsInHistogram)
c      end do

c------------------------------------------------
c Print out some info for debugging purposes

c!!!!
c      goto 510

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
      print *,'Topology','total'
      do j=1,numatomtypesTop
       print 43,j,numatomsInCategoryTotal(j,3)
      end do
      print *
      print *
      print *,'Pairwise data'
      print 29,(i,i=1,numatomtypes)
      do i=1,numatomtypes
       print 18,i,(pairwiseEvents(i,j),j=1,numatomtypes)
      end do

      print *
      print *,'Expected contacts'
      print 29,(i,i=1,numatomtypes)
      do i=1,numatomtypes
       print 14,i,(fExpected(i,j),j=1,numatomtypes)
      end do
      print *

      print *
      print *,'Pairwise function ln(actual/expected)'
      print 29,(i,i=1,numatomtypes)
      do i=1,numatomtypes
      print 34,i,(fPairwise(i,j),j=1,numatomtypes)
      end do
      print *

      print *
      print *,'Correlation between atom types (based on pairwise data)'
      print 29,(i,i=1,numatomtypes)
      do i=1,numatomtypes
       print 34,i,(R(i,j),j=1,numatomtypes)
      end do

      call PMFmodule(pmfShellCounter,numatomtypes,
     &   numatomtypesEnv,numdistances,PMFENERGY)

      print *,'Back from PMFmodule'

      do i1 = 1,numatomtypes
       do i2 = 1,numatomtypesEnv
        print *,i1,i2,(PMFENERGY(i1,i2,i3),
     &   i3 = 1,numdistances)
       end do
      end do


      call writeenergy(energyfile,EwordAlternate,Eword,singleBodyEnergy,
     &  comboTypes,fPairwise,PMFENERGY,EPair,contactdir)

c      print *,'Im back from writeenergy'
499   continue

      do i1 = 1,numatomtypes
       do i2 = 1,numatomtypesEnv
        pmfAvg(i1,i2) = 0.0
        do i3 = 1,pmfAvgMaxshell
         pmfAvg(i1,i2) = pmfAvg(i1,i2) +
     &    PMFENERGY(i1,i2,i3)         
        end do
       end do
      end do

      print *
      print *,'calculation of pmfAvg'
      print *

c      print *,'before division'
c      do i1 = 1,numatomtypes
c       do i2 = 1,numatomtypesEnv
c        print *,i1,i2,pmfAvg(i1,i2)
c       end do
c      end do


      do i1 = 1,numatomtypes
       do i2 = 1,numatomtypesEnv
        pmfAvg(i1,i2) = pmfAvg(i2,i2)/
     &     pmfAvgMaxshellReal
        print *,i1,i2,pmfAvg(i1,i2)
       end do
      end do

      if (.not.mfexist) then
510   goto 999 
      end if
c---------------------------------------------------
      open (unit=4,file = misfolddir//pdblistname2)

      print *,misfolddir//pdblistname2
      print *
      print *
      do i=1,2000
       numGood(i) = 0.0
       numBad(i) = 0.0
       TotalEnergy(i) = 0.0
       TotalEnergyPairwise(i) = 0.0
       TotalEnergyAlternate(i) = 0.0
       TotalEnergySingleBody(i) = 0.0
       TEnergyPair(i) = 0.0
       ETcomb(i) = 0.0
       pmfTesttotal(i) = 0.0
       TESubset(i) = 0.0
       TEASubset(i) = 0.0
       pmfSubset(i) = 0.0
       do j=1,10
        pmfTestenergy(i,j) = 0.0
       end do
      end do
      whichPdb = 1

c      print 73,'FILE_NAME                                            
c     +Numatm Numres Energy
c     +     AltEnergy   Pairwise   SingleBody '

       print 74,'FILE_NAME','Energy','AltEnergy','PMFAdd',
     &  'Esubset','AESubset','PMFSub'

      firstTimeThrough = .true.

500   read (4,12,end = 998) pdbname 

      print *,'pdbname = ',pdbname 
c-------------------------------------------------------
c Initialize all protein-specific variables before reading
c new pdb file

      Undetermined(whichPdb) = 0.0
      numres = 0
      numatm = 0
      do i=1,maxNumRes
       energyPerResidue(i) = 0.0
       energyPerResidueAlternate(i) = 0.0
       energyPerResiduePairwise(i) = 0.0
      end do
      do i=1,maxat
       multiconf(i) = ''
       SigBoolAt(i) = 1
       do j=1,7
        atomData(i,j) = 0
       end do
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
         call pdbread (pdbname,misfolddir)
       print *,'pdbread is over'
       if (reject) then
        goto 500
       end if


      if (firstTimeThrough) then
       numatmToSample = anint(numatm*(percentToSample/100))
       call random_subset(numatm,numatmToSample,iseed,subset)
       print *,'random_subset is over'
       firstTimeThrough = .false.
      end if

      badAtom = .false.
      call categorizeatoms (numatomtypesEnv,numatomtypesTop,
     &   categorylistname,contactdir,mapAtomTopology,firstTimeThrough,
     &   badAtom,category)
      print *,'categorize_atoms is over'
      if (badAtom) then
       print *,'Skipping file ',filename,'because of bad atom'
       goto 500
      end if 

c-----------------------------------------------------------------
c Now calculate distances, relevancy of pairwise contact
c (based in part on connectivity and atom type)

      sscount = 1
      do i1 = 1,numatm
       if (i1.eq.subset(sscount)) then
        foundOne = .true.
        sscount = sscount + 1
       else
        foundOne = .false.
       endif
       atomData(i1,1) = category(i1,1)
       backboneCounted = 0
       do i2 = 1,numatm
        call Relevant_Distance(i1,i2,numskip,res(i1),
     &     res(i2),atomname(i1),atomname(i2),
     &     chain(i1),chain(i2),resnum,relevant)
         tempdistance = Distance(coord(1,i1),coord(2,i1),
     &     coord(3,i1),coord(1,i2),
     &     coord(2,i2),coord(3,i2))
        if (relevant) then
         if (tempdistance.lt.criticalmaxdistance.and.
     &       tempdistance.gt.criticalmindistance) then
          call Make_Bin_Decision(tempdistance,pmfMaxdis,pmfMindis,
     &       pmfStep,pmfShell)
          pmfTestenergy(whichPdb,i1) = 
     &      pmfTestenergy(whichPdb,i1) + 
     &      (PMFENERGY(category(i1,1),
     &       category(i1,2),pmfShell) -
     &      pmfAvg(category(i1,1),category(i1,2))) 
          pmfTesttotal(whichPdb) = pmfTesttotal(whichPdb) +
     &      (PMFENERGY(category(i1,1),
     &       category(i1,2),pmfShell) -
     &      pmfAvg(category(i1,1),category(i1,2)))
          if (foundOne) then
           pmfSubset(whichPdb) = pmfSubset(whichPdb) +
     &      (PMFENERGY(category(i1,1),
     &       category(i1,2),pmfShell) -
     &      pmfAvg(category(i1,1),category(i1,2)))
          end if


c Note! Pairwise energy calculation does NOT take into account the skipping
c of backbone residues! This is an artificial construction to improve statistics.
          energyPerResiduePairwise(resnum(i1)) = 
     &     energyPerResiduePairwise(resnum(i1)) +
     &     (fPairwise(category(i1,1),category(i2,1)) / 2) 
          TotalEnergyPairwise(whichPdb) = 
     &     TotalEnergyPairwise(whichPdb) +
     &     fPairwise(category(i1,1),category(i2,1)) / 2
          TEnergyPair(whichPdb) = TEnergyPair(whichPdb) +
     &     EPair(category(i1,1),category(i2,2))
          if (category(i2,2).eq.backboneCategory) then
           if (mod(backboneCounted,numBackboneToCount).eq.0) then
            Coordination(i1) = Coordination(i1) + 1
            atomData(i1,2) = atomData(i1,2) + 1
            atomData(i1,(category(i2,2)+2)) =
     &        atomData(i1,(category(i2,2)+2)) + 1
           end if
          backboneCounted = backboneCounted + 1
          else
           Coordination(i1) = Coordination(i1) + 1
           atomData(i1,2) = atomData(i1,2) + 1
           atomData(i1,(category(i2,2)+2)) =
     +       atomData(i1,(category(i2,2)+2)) + 1
          end if
          pairwiseEvents(category(i1,1),category(i2,1)) =
     +     pairwiseEvents(category(i1,1),category(i2,1)) + 1
         end if
        end if
       end do
      end do

    

c Now we'll group the atoms based on atomtype(1-20),coordination(0-13),
c and the #s of atoms in each of the four environment variables.


c      print *,'Resnum Res Chain Name  category coord    #A      
c     & #B      #C      #D     EwordAlt  Sig'

      sscount = 1
      do i1 = 1,numatm
       if (i1.eq.subset(sscount)) then
        foundOne = .true.
        sscount = sscount + 1
       else
        foundOne = .false.
       endif
       if (atomData(i1,2).le.criticalCoordination) then
        TotalEnergySingleBody(whichPdb) = 
     &  TotalEnergySingleBody(whichPdb) +
     &  singleBodyEnergy(category(i1,1),atomData(i1,2)+1)
        do i2 = 1,comboTypes(atomData(i1,2)+1)
         if (atomData(i1,3).eq.comboTypesComp(atomData
     &   (i1,2)+1,i2,1)) then
          if (atomData(i1,4).eq.comboTypesComp(atomData
     &    (i1,2)+1,i2,2)) then
           if (atomData(i1,5).eq.comboTypesComp(atomData
     &     (i1,2)+1,i2,3)) then
            if (atomData(i1,6).eq.comboTypesComp(atomData
     &      (i1,2)+1,i2,4)) then
              if (SignificanceBool(atomData(i1,1),atomData
     &         (i1,2)+1,i2).eq.0) then
               Undetermined(whichPDB) = Undetermined(whichPDB) + 1
              end if
c!!!!!
c             print 72,resnum(i1),res(i1),chain(i1),
c     &        atomData(i1,1),atomname(i1),
c     &        atomData(i1,2),atomData(i1,3),atomData(i1,4),
c     &        atomData(i1,5),atomData(i1,6),
c     &        EwordAlternate(atomData(i1,1),atomData(i1,2)+1,i2),
c     &        SignificanceBool(atomData(i1,1),atomData(i1,2)+1,i2)

              TotalEnergy(whichPdb) = TotalEnergy(whichPdb) +
     +          Eword(atomData(i1,1),atomData(i1,2)+1,i2)
              TotalEnergyAlternate(whichPdb) = 
     &         TotalEnergyAlternate(whichPdb) +
     &          EwordAlternate(atomData(i1,1),
     &          atomData(i1,2)+1,i2)
              if (foundOne) then
               TESubset(whichPdb) = TESubset(whichPdb) +
     &          Eword(atomData(i1,1),atomData(i1,2)+1,i2)
               TEASubset(whichPdb) = TEASubset(whichPdb) +
     &          EwordAlternate(atomData(i1,1),
     &          atomData(i1,2)+1,i2)
              end if

              if (EwordAlternate(atomData(i1,1),
     &         atomData(i1,2)+1,i2).lt.0.00) then
                numGood(whichPdb) = numGood(whichPdb) + 1.00
              else
                numBad(whichPdb) = numBad(whichPdb) + 1.00
              end if
              energyPerResidue(resnum(i1)) = 
     &         energyPerResidue(resnum(i1)) +
     &         Eword(atomData(i1,1),atomData(i1,2)+1,i2)
              energyPerResidueAlternate(resnum(i1)) = 
     &         energyPerResidueAlternate(resnum(i1)) +
     &         EwordAlternate(atomData(i1,1),
     &         atomData(i1,2)+1,i2)
            end if
           end if
          end if
         end if
        end do
       else
c     - Approximation to lower coordination composition may eventually be inserted here
        if (SignificanceBool(atomData(i1,1),15,1).eq.0) then
         Undetermined(whichPDB) = Undetermined(whichPDB) + 1
        end if
c!!!!!
c             print 72,resnum(i1),res(i1),chain(i1),
c     &        atomData(i1,1),atomname(i1),
c     &        atomData(i1,2),atomData(i1,3),atomData(i1,4),
c     &        atomData(i1,5),atomData(i1,6),
c     &        EwordAlternate(atomData(i1,1),15,1),
c     &        SignificanceBool(atomData(i1,1),15,1)
c!!!!!

        TotalEnergy(whichPdb) = TotalEnergy(whichPdb) + 
     +    Eword(atomData(i1,1),15,1)
        TotalEnergyAlternate(whichPdb) = 
     &  TotalEnergyAlternate(whichPdb) +
     &    EwordAlternate(atomData(i1,1),15,1)
        if (foundOne) then
         TESubset(whichPdb) = TESubset(whichPdb) +
     &     Eword(atomData(i1,1),15,1)
         TEASubset(whichPdb) = TEASubset(whichPdb) +
     &     EwordAlternate(atomData(i1,1),15,1)
        end if
        if (EwordAlternate(atomData(i1,1),15,1).lt.0.00) then
         numGood(whichPdb) = numGood(whichPdb) + 1.00
        else
         numBad(whichPdb) = numBad(whichPdb) + 1.00
        end if
        energyPerResidue(resnum(i1)) = 
     &   energyPerResidue(resnum(i1)) +
     &    Eword(atomData(i1,1),15,1)
        energyPerResidueAlternate(resnum(i1)) = 
     &    energyPerResidueAlternate(resnum(i1)) +
     &    EwordAlternate(atomData(i1,1),15,1)
       end if
      end do


      ETcomb(whichPdb) = A1 * TotalEnergySingleBody(whichPdb) +
     &  A2 * TotalEnergyPairwise(whichPdb) +
     &  A3 * TotalEnergyAlternate(whichPdb)
 

c997   print 40,pdbname,numatm,numres,TotalEnergy(whichPdb),
c     & TotalEnergyAlternate(whichPdb),TotalEnergyPairwise(whichPdb),
c     & TotalEnergySingleBody(whichPdb),pmfTesttotal(whichPdb),
c     & numGood(whichPdb),numBad(whichPdb)


997   print 76,pdbname,TotalEnergySingleBody(whichPdb),
     & TEnergyPair(whichPdb),TotalEnergyAlternate(whichPdb),
     & ETcomb(whichPdb)


c997   print 75,pdbname,TotalEnergy(whichPdb),
c     & TotalEnergyAlternate(whichPdb),
c     & pmfTesttotal(whichPdb),
c     & TESubset(whichPdb),
c     & TEASubset(whichPdb),
c     & pmfSubset(whichPdb)

      numatmReal = float(numatm)
      energyPerAtomMean(whichPdb) = TotalEnergy(whichPdb) / numatmReal
      energyPerAtomAlternateMean(whichPdb) = 
     &  TotalEnergyAlternate(whichPdb) / numatmReal
      energyPerResidueMean(whichPdb) = 0.0
      energyPerResidueAlternateMean(whichPdb) = 0.0

      do i1 = 1,numres
       if (perResidue) then
        print *,i1,energyPerResidue(i1),
     &    energyPerResidueAlternate(i1)
       end if
        energyPerResidueMean(whichPdb) = 
     &   energyPerResidueMean(whichPdb) +
     &   energyPerResidue(i1)
        energyPerResidueAlternateMean(whichPdb) =
     &   energyPerResidueAlternateMean(whichPdb) +
     &   energyPerResidueAlternate(i1)
      end do
      energyPerResidueMean(whichPdb) = energyPerResidueMean(whichPdb) /
     &  numres
      energyPerResidueAlternateMean(whichPdb) = 
     & energyPerResidueAlternateMean(whichPdb) / numres

      if (energyPerResidueMean(whichPdb).eq.0.0) then
       print *,'Zero mean',filename,numres
       do i2 = 1,numres
        print *,energyPerResidue(i2)
       end do
      end if
      minenergy = -15.0
      maxenergy = 15.0
      energyStepsize = 1.0
c      call BinEnergyPerResidue(numres,energyPerResidue,
c     &  energyPerResidueHistogram,
c     &  whichPdb,minenergy,maxenergy,energyStepsize)
c      call BinEnergyPerResidue(numres,energyPerResidueAlternate,
c     &  energyPerResidueAlternateHistogram,whichPdb,
c     &  minenergy,maxenergy,energyStepsize)

      whichPdb = whichPdb + 1
      goto 500

998   close (4)
      whichPdb = whichPdb -1


c      print *,'EnergyPerResidueHistogram'
c      do i=1,whichPdb
c        print 69,(energyPerResidueHistogram(i,j),j=1,30)
c      end do
c      print *
c      print *,'EnergyPerResidueAlternateHistogram'
c      do i=1,whichPdb
c        print 69,(energyPerResidueHistogram(i,j),j=1,30)
c      end do


c      print *
c      print *,'Normal_Mean_perRes','      Alternate_Mean_perRes',
c     + '    Normal_Mean_perAtom','  Alternate_Mean_perAtom',
c     + '   #_Undetermined_Atoms'
c      do i=1,whichPdb
c        print *,energyPerResidueMean(i),energyPerResidueAlternateMean(i),
c     +   energyPerAtomMean(i),energyPerAtomAlternateMean(i),Undetermined(i)
c      end do
      numpdbsreal = float(whichPdb-1)


c------------------------------------------------
c Z score calculation for misfolds
 
      do i=2,whichPdb
       Emean = Emean + TotalEnergy(i)
       EmeanA = EmeanA + TotalEnergyAlternate(i)
       EmeanSB = EmeanSB + TotalEnergySingleBody(i)
c       EmeanP = EmeanP + TotalEnergyPairwise(i)
       EmeanP = EmeanP + TEnergyPair(i)
       EmeanT = EmeanT + ETcomb(i)
      end do

      Emean = Emean /  numpdbsreal
      EmeanA = EmeanA / numpdbsreal
      EmeanSB = EmeanSB / numpdbsreal
      EmeanP = EmeanP / numpdbsreal
      EmeanT = EmeanT / numpdbsreal

      do i=1,whichPdb
       sigma = sigma + (TotalEnergy(i) - Emean)**2
       sigmaA = sigmaA + (TotalEnergyAlternate(i) - EmeanA)**2
       sigmaSB = sigmaSB + (TotalEnergySingleBody(i) - EmeanSB)**2
c       sigmaP = sigmaP + (TotalEnergyPairwise(i) - EmeanP)**2
       sigmaP = sigmaP + (TEnergyPair(i) - EmeanP)**2
       sigmaT = sigmaT + (ETcomb(i) - EmeanT)**2
      end do

       sigma = sqrt(sigma/numpdbsreal) 
       sigmaA = sqrt(sigmaA/numpdbsreal)
       sigmaSB = sqrt(sigmaSB/numpdbsreal)
       sigmaP = sqrt(sigmaP/numpdbsreal)
       sigmaT = sqrt(sigmaT/numpdbsreal)

       Z = (Emean - TotalEnergy(1))/sigma
       ZA = (EmeanA - TotalEnergyAlternate(1))/sigmaA
       ZSB = (EmeanSB - TotalEnergySingleBody(1))/sigmaSB
c       ZP = (EmeanP - TotalEnergyPairwise(1))/sigmaP
       ZP = (EmeanP - TEnergyPair(1))/sigmaP
       ZT = (EmeanT - ETcomb(1))/simgaT
   
       print *
       print *,'numpdbsreal = ',numpdbsreal
       print *

c       print *,'Emean = ',Emean
c       print *,'Sigma = ',sigma
c       print *,'Z = ',Z
c       print *

c       print *,'EmeanA = ',Emean
c       print *,'SigmaA = ',sigmaA
c       print *,'ZA = ',ZA
c       print *

       print *,'EmeanSB = ',EmeanSB
       print *,'SigmaSB = ',sigmaSB
       print *,'ZSB = ',ZSB
       print *

       print *,'EmeanP = ',EmeanP
       print *,'SigmaP = ',sigmaP
       print *,'ZP = ',ZP

       print *,'EmeanA = ',Emean
       print *,'SigmaA = ',sigmaA
       print *,'ZA = ',ZA
       print *




999   end
c--------------------------------------------------------

      subroutine Relevant_Distance(i1,i2,
     & numskip,residuename1,residuename2,atomname1,atomname2,
     & chain1,chain2,resnum,relevant)

      include 'parameters.h'
      integer i1,i2,resnum(maxat),numskip
      character*4 atomname1,atomname2
      character*3 residuename1,residuename2
      character*1 chain1,chain2
      logical relevant,sameChain

c------------
c  Is the distance between two atoms (designated by the variables
c  i1 and i2) relevant to our counting statistics?
c  Tests for relevance:
c    Current version:  if i2 is  within numskip residues from the
c       residue containing the atom of interest (denoted by i1)
c       then the variable relevant receives the value .false. 
c       Otherwise, the value is .true.
c-------------

      sameChain = .true.
      relevant = .true. 

      if (chain1.ne.chain2) then
       sameChain = .false.
      end if

c--Case where resnum(i2) is less than resnum(i1)

      if (resnum(i2).lt.resnum(i1)) then
       if (resnum(i2).gt.resnum(i1)-numskip) then
         relevant = .false.
       end if      

       if (resnum(i1)-resnum(i2).eq.1.and.atomname1.eq.
     +  ' N  '.and.atomname2.eq.' C  ') then
          relevant = .false.
       end if
      end if

c--Case where resnum(i2) is greater than resnum(i1)

      if (resnum(i2).ge.resnum(i1)) then
       if (resnum(i2).lt.resnum(i1)+numskip) then
          relevant = .false.
       end if
       if (resnum(i2)-resnum(i1).eq.1.and.atomname1.eq.
     +   ' C  '.and.atomname2.eq.' N  ') then
         relevant = .false.
       end if
      end if

      if (numskip.eq.0) then
       relevant = .true.
      end if

      if (.not.sameChain) then
       relevant = .true.
      end if

      end

c----------------------------------------------------------------------

      subroutine correlate(x,y,n,r)

c x and y are corresponding arrays of type real (maximum of 100 units)
c n is an integer value designating the actual length of the x and y arrays
c r is the returned correlation coefficient

      include 'parameters.h'
 
      integer i,n
      real x(maxnumatomtypes)
      real y(maxnumatomtypes)
      real r
      real xdiff(100)
      real ydiff(100)
      real*8 xbar,ybar,numerator,xdenom,ydenom,denominator

c     initialize varibles

      do i=1,maxnumatomtypes
       xdiff(i) = 0.0
       ydiff(i) = 0.0
      end do
      r = 0.0
      xbar = 0.0
      ybar = 0.0
      xdenom = 0.0
      ydenom = 0.0
      numerator = 0.0
      denominator = 0.0

      do i=1,n
       xbar = xbar + x(i)
       ybar = ybar + y(i)
      end do
      xbar = xbar / n
      ybar = ybar / n

      do i=1,n
       xdiff(i) = x(i) - xbar
       ydiff(i) = y(i) - ybar
      end do

      do i=1,n
       numerator = numerator + (xdiff(i) * ydiff(i))
       xdenom = xdenom + (xdiff(i)**2)
       ydenom = ydenom + (ydiff(i)**2)
      end do
      denominator = sqrt(xdenom * ydenom)
      r = numerator/denominator
      end

c------------------------------------------------------------------
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


c*************************************
      function Fact(foobar)

      double precision Fact
      integer i,foobar

      Fact = 1
      if (foobar.ge.2) then
       do i=2,foobar
        Fact = Fact*i
       end do
      end if
      end

c*************************************
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

c       do j=1,30
c        print 11,(ChiSqrTable(i,j),i=1,13)
c       end do

      end



c
c **************************************************************
c
      subroutine gettheinteger(istr,str)

      character*80 str
        character ch
        integer numvar,istr,ich,power
        integer begin(10),end(10)

        call parseline(str,begin,end,numvar)
        istr = 0
        if(numvar.gt.2)then
           write(6,*)'Error parsing  argument in gettheinteger()'
      else
           do i = end(1),begin(1),-1
                ch = str(i:i)
                ich = ichar(ch) - 48
                power = 10**(end(1)-i)
                istr = istr + power*ich
         enddo
      endif

        return
        end
c
c *************************************************************
c
        subroutine parseline(longstring,begin,end,numvar)

        character*80 longstring
        character*1 char
        integer mxwordsinstring
        parameter (mxwordsinstring=10)
        integer begin(mxwordsinstring),end(mxwordsinstring)
        integer begincount,endcount,numvar

        begincount = 0
        endcount = 0
        do i = 1,80
           char = longstring(i:i)
           if(begincount .eq. endcount)then
c     -     I am looking for the beginning of an entry.
                if(char.ne.' ')then
                   begincount = begincount + 1
                   begin(begincount) = i
            endif
         elseif (begincount .eq. (endcount + 1))then
c     -     I am looking for the end of an entry.
                if(char.eq.' ')then
                   endcount = endcount + 1
                   end(endcount) = i - 1
            else
                   if(i.eq.80)then
                        endcount = endcount + 1
                        end(endcount) = i
               endif
            endif
         else
                write(6,*)'Error in subroutine parseline().'
                goto 20
         endif
      enddo

        numvar = begincount

        return
 20   end
c
c ****************************************************************
c
      subroutine getthereal(rstr,str)

      integer mxcmdlinevar
        parameter (mxcmdlinevar=10)
        integer begin(mxcmdlinevar),end(mxcmdlinevar)
        integer ich,power,decimalptindx,tempint,numvar
        character*80 str
        character ch,chch
        real parity
        real rstr,leftofdecimal,rightofdecimal

        call parseline(str,begin,end,numvar)

c     Begin by checking if the argument is negative.
      i = begin(1)
      chch = str(i:i)
        if(chch.eq.'-')then
           parity = -1.0
           begin(1) = begin(1) + 1
      else
           parity = 1.0
      endif
c
        decimalptindx = end(1) + 1
        if(numvar.gt.2)then
           write(6,*)'Error parsing argument in getthereal()'
      else
           do i = begin(1),end(1)
                ch = str(i:i)
                ich = ichar(ch)
                if(ich.eq.46)then   ! 46 is ascii code for the decimal point
                   decimalptindx = i
            endif
         enddo
         tempint = 0
           if(decimalptindx.gt.begin(1))then
                ii = decimalptindx - 1
              do i = ii,begin(1),-1
                   power = 10**(ii - i)
                   ch = str(i:i)
                   ich = ichar(ch) - 48
                   tempint = tempint + ich*power
            enddo
              leftofdecimal = float(tempint)
         else
                leftofdecimal = 0.0
         endif
           tempint = 0
           ii = decimalptindx + 1
           if(decimalptindx.lt.end(1))then
            do i = end(1),ii,-1
                   power = 10**(end(1)-i)
                   ch = str(i:i)
                   ich = ichar(ch) - 48
                   tempint = tempint + ich*power
            enddo
                rightofdecimal = float(tempint)/float(10*power)
         else
                rightofdecimal = 0.0
         endif
           rstr = leftofdecimal + rightofdecimal
      endif
        rstr = parity*rstr

      return
        end

c--------------------------------------------------------
      subroutine Make_Bin_Decision(tempdistance,maxdistance,
     &    mindistance,distanceStepsize,whereToBin)

      real tempdistance,maxdistance,mindistance,distanceStepsize,
     +     templowerlimit,tempupperlimit
      integer whereToBin

      templowerlimit = mindistance

      if (tempdistance.le.mindistance) then
       goto 500
      end if

      whereToBin = 1
      do tempupperlimit = mindistance+distanceStepsize,maxdistance,
     +   distanceStepsize
       if (tempdistance.gt.templowerlimit.and.tempdistance.le.
     +     tempupperlimit) then
         return
       end if
       whereToBin = whereToBin + 1
       templowerlimit = tempupperlimit
      end do

c------------------------
c If subroutine reaches this point, interaction falls
c outside the bounds of the specified max and min distances.
c A value of 999 acts a flag to tell the calling routine that
c this is a special case.

500   whereToBin = 999

      end

c----------------------------------------------------------------------------
      subroutine BinEnergyPerResidue(numres,energy,energyHistogram,
     +       whichPdb,minenergy,maxenergy,energyStepsize)

      include 'parameters.h'

      integer numres,i,whichPdb,whereToBin
      real minenergy,maxenergy,energyStepsize
      real templowerlimit,tempupperlimit
      real energyHistogram(200,50)
      real energy(maxNumRes)

      do i=1, numres
        templowerlimit = minenergy

        whereToBin = 1

        do tempupperlimit = minenergy+energyStepsize,maxenergy,
     +   energyStepsize

         if (energy(i).gt.templowerlimit.and.energy(i).le.
     +     tempupperlimit) then
         energyHistogram(whichPdb,whereToBin) = 
     +      energyHistogram(whichPdb,whereToBin) + 1
         end if

        whereToBin = whereToBin + 1
        templowerlimit = tempupperlimit
       end do  

      end do
      end
