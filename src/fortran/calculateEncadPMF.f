      program calculateEncadPMF

c Author: C. Summa
c-------------------------------------------------------
c Variable declarations

	

      include 'parameters.h'
      include 'pdb.h'
      include 'encad.h'

      doublePrecision Fact
      real Distance
      real tempdistance
      real pmfStepSize
      real R(maxnumatomtypes,maxnumatomtypes)
      real realPairwiseEvents(maxnumatomtypes,maxnumatomtypes)
      real SumTotalEvents
      real pairwiseEvents(maxnumatomtypes,maxnumatomtypes)
      real totalInteractions
      real minenergy,maxenergy,energyStepsize
      real overallSumEvents
      real version
      real tempreal
      real Undetermined(maxNumPDBs)
      real numatmReal
      real numpdbsreal
      real pmfShellCounter(maxnumatomtypes,maxnumatomtypesEnv,
     &       pmfMaxNumBins)
      real pmfShellCounterHiRes(maxnumatomtypes,maxnumatomtypes,
     &       pmfMaxNumBins)
      real pmfMaxdis
      real pmfMindis
      real PMFENERGY(maxnumatomtypes,maxnumatomtypesEnv,pmfMaxNumBins)
      real PMFHiRes(maxnumatomtypes,maxnumatomtypes,pmfMaxNumBins)
      real pmfTestenergy(maxNumPDBs,maxat2)
      real pmfTesttotal(maxNumPDBs)
      real randm
      real ran1
      real*8 numatmTotal2Real

      real    TESubset(maxNumPDBs)
      real    TEASubset(maxNumPDBs)
      real    pmfSubset(maxNumPDBs)
      integer sscount
      logical foundOne

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
      integer  tempInteger
      integer  numatomtypes
      integer  numatomtypesEnv
      integer  numatomtypesTop
      integer  atomData(maxat,7)
      integer  probSumHist(10)
      integer  numA,numB,numC,numD
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
      integer  numdistances
      integer  StringLength
      integer*8  numatmTotal2
      integer*8  Pmult
      integer*8  SignificanceHistogram(14)
      character*80  numskipstr,crdiststr,iseedstr
      character*80  percentstr
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

      integer templen
      integer tempstart
      integer tempend
      character*120 tstring
      character*120 dbconcat
      character*120 mfconcat
      character*120 enconcat
      logical dbexist
      logical mfexist
      logical enexist

      real func1
      real treal

      real BestEnergy
      real WorstEnergy
      real numbestnearhetero
      real numworstnearhetero
      real atomEnergy(maxat)
      integer sortedMap(maxat)
      integer bestarray(10)
      integer worstarray(10)
      integer numtoflag

      character*80 lowDistStr
      character*80 highDistStr
      character*80 stepStr

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
     & 1x,a1,1x,f8.1,a1,f8.1,2x,f6.2)
79    format (f5.2,1x,f5.2,1x,f5.2,1x,f5.2,1x,f5.2,1x,f5.2,1x,
     & f5.2,1x,f5.2,1x,f5.2,1x,f5.2,1x,f5.2,1x,f5.2,1x,
     & f5.2,1x,f5.2,1x,f5.2,1x,f5.2,1x,f5.2,1x,f5.2,1x,
     & f5.2,1x,f5.2,1x,f5.2,1x,f5.2,1x,f5.2,1x,f5.2,1x,
     & f5.2,1x,f5.2,1x,f5.2,1x,f5.2,1x,f5.2,1x,f5.2,1x,
     & f5.2,1x,f5.2,1x,f5.2,1x,f5.2,1x,f5.2,1x,f5.2,1x,
     & f5.2,1x,f5.2,1x,f5.2,1x,f5.2,1x,f5.2,1x,f5.2,1x,
     & f5.2,1x,f5.2,1x,f5.2,1x,f5.2,1x,f5.2,1x,f5.2,1x,
     & f5.2,1x,f5.2,1x,f5.2,1x,f5.2,1x,f5.2,1x,f5.2,1x,
     & f5.2,1x,f5.2,1x,f5.2,1x,f5.2,1x,f5.2,1x,f5.2,1x)

c------------------------------------------------------------
c Read in command line arguments

      numarguments = iargc()

      if (numarguments.ne.5) then
        print *,'Usage: calculateHiResPMF listname #res_skip
     &lowDist highDist step' 
        goto 999
      endif

      call getarg(1,selectname)
      call getarg(2,numskipstr)
      call gettheinteger(numskip,numskipstr)
      call getarg(3,lowDistStr)
      call getthereal(pmfMindis,lowDistStr)
      call getarg(4,highDistStr)
      call getthereal(pmfMaxdis,highDistStr)
      call getarg(5,stepStr)
      call getthereal(pmfStepSize,stepStr)

      call getenv('CONTACTDIR',contactdir)
      call getenv('CONTACTPDBDIR',contactpdbdir)
      call getenv('CONTACTPDBDIR2',databasedir)
      call getenv('MISFOLDDIR',misfolddir)

c------------------------------------------------
c print out some of the input data and values
c of relevant variables

      print *,'PDBSelect list used to generate potential: ',selectname
      print *,'Number of residues to skip: ',numskip
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


      numdistances = ((pmfMaxdis - pmfMindis) /pmfStepSize)
c      numdistances = ((pmfMaxdis - pmfMindis) /pmfStepSize) + 1
      numatomtypes = 26

      A1 = 1.00
      A2 = 1.00
      A3 = 1.00


      print *,'pmfMaxdis = ',pmfMaxdis
      print *,'pmfMindis = ',pmfMindis
      print *,'pmfStepSize = ',pmfStepSize
      print *,'numdistances = ',numdistances

      categorylistname = 'FullAtomList9.dat'
     
      firstTimeThrough = .true.
      null = 0
      overallSumEvents = 0.0
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

 
      do i=1,maxnumatomtypesEnv
       do j = 1, maxnumatomtypes
        pair(j,i) = 0.00
       end do
      end do

      numatmTotal2 = 0
      numenvTotal = 0.0
      pdbsread = 0
      numchainused = 0
      totalInteraction = 0

c--------------------------------------------------------

      if (firstTimeThrough) then
       if (dbexist) then
        call Read_pdblist(dbconcat)
       end if
      end if

c-------------------------------------------------------
c Check for existence of energy file corresponding to user-defined
c parameters.  If it exits (i.e. iocode = 0 ) then read in values
c and skip the derivation.

c-------------------------------------------------------

      print *
      print *,'Beginning Derivation....'
      do ii = 1, numencadpdbnames
       filename = encadpdbnames(ii)
       rejectflag(ii) = .false.

c-------------------------------------------------------
c Initialize all protein-specific variables before reading
c new pdb file

      numres = 0
      numatm = 0
      do i=1,maxat
       multiconf(i) = ''
       unknown(i) = .false.
       do j=1,7
        atomData(i,j) = 0
       end do
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

       call pdbread_clean_nobiomt_with_h (filename,databasedir)
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

      call categorizeEncadatoms (firstTimeThrough, badAtom,
     &   categorylistname,contactdir,category)
 

      if (badAtom) then
       call pdbcollapse(category)
      end if

      call findRelevantChain (pdbChain(ii))

      badAtom = .false.
  
      call pdbcheck (contactdir,filename,reason,firsttime)

      if (badAtom) then
       print *,'weve still got a bad atom in file ',filename
      end if

c-----------------------------------------------------------------
c Now calculate distances, relevancy of pairwise contact 
c (based in part on connectivity and atom type)


       chainstart = 1
       chainend = numatm
       do i1 = chainstart,chainend
        if (.not.unknown(i1)) then
         do i2 = 1,numatm
          if (.not.unknown(i2)) then
           call relevantDistance2(i1,i2,numskip,relevant)
           if (relevant) then
            tempdistance = Distance(coord(1,i1),coord(2,i1),
     &          coord(3,i1),coord(1,i2),coord(2,i2),coord(3,i2))
            call Make_Bin_Decision(tempdistance,pmfMaxdis,pmfMindis,
     &          pmfStepSize,pmfShell)
            if (pmfShell.ne.999) then
              call pmfAddIt(pmfShellCounterHiRes,pmfShell,
     &          category(i1,1),category(i2,1),i1,i2,
     &          chainStart,chainEnd)
            end if
	   end if
          end if
         end do
        end if
       end do

c Go back and see if there another pdb file to read in!

      pdbsread = pdbsread + 1
      firstTimeThrough = .false.

c End of major loop
c-----------------------------------------------------------
222   end do

      close (4)

      call PMFmodule(pmfShellCounterHiRes,numatomtypes,
     &   numatomtypes,numdistances,PMFHiRes,
     &   'PMF_rawdata.dat')

      print *,'Back from PMFmodule'

      call writePMF('PMF_hires_symmetric.dat',PMFHiRes,numatomtypes,
     &   numatomtypes,numdistances)

      open (unit=5,file='PMFDistanceBinUpperLimits.dat')
        write (5,79) (pmfMinDis + (i*pmfStepSize),i=1,numdistances)
      close (5)

499   continue
999   end

c------------------------------------------------------------------
      subroutine Make_Bin_Decision(tempdistance,maxdistance,
     &    mindistance,distanceStepsize,whereToBin)

      real tempdistance,maxdistance,mindistance,distanceStepsize,
     +     templowerlimit,tempupperlimit
      integer whereToBin
      integer counter

      templowerlimit = mindistance
      tempupperlimit = mindistance + distanceStepsize
      counter = 1
200   if (tempdistance.ge.templowerlimit.and.
     &    tempdistance.lt.tempupperlimit) then
        whereToBin = counter
        return
      end if

      templowerlimit = tempupperlimit
      tempupperlimit = tempupperlimit + distanceStepsize
      counter = counter + 1
      if (tempupperlimit.le.maxdistance + 0.01) then
        goto 200
      end if

c If subroutine reaches this point, interaction falls
c outside the bounds of the specified max and min distances.
c A value of 999 acts a flag to tell the calling routine that
c this is a special case.

500   whereToBin = 999
      return

      end

c--------------------------------------------------------------------------
      subroutine pmfAddIt(pmfShellCounterHiRes,pmfShell,
     &         category1,category2,num1,num2,
     &         chainStart,chainEnd)

      include 'parameters.h'
	  integer pmfShell
	  integer category1,category2,chainStart,chainEnd
      integer tempcat1,tempcat2
	  integer num1,num2
	  integer temp
      real pmfShellCounterHiRes(maxnumatomtypes,maxnumatomtypes,
     &       pmfMaxNumBins)

      if (category1.gt.category2) then
		tempcat1 = category2
		tempcat2 = category1
      else
        tempcat2 = category2
        tempcat1 = category1
      end if

	  if (num2.lt.chainStart.and.num2.gt.chainEnd) then
	   if (num1.lt.num2) then
	    pmfShellCounterHiRes(category1,category2,pmfShell) =
     &     pmfShellCounterHiRes(category1,category2,pmfShell) + 1.0
        return
	   end if
	  else
	   pmfShellCounterHiRes(category1,category2,pmfShell) =
     &   pmfShellCounterHiRes(category1,category2,pmfShell) + 1.0
       return
	  end if
	  end
