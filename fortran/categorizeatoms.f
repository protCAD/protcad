
      subroutine categorizeatoms (numatomtypesEnv,numatomtypesTop,
     &   categorylistname,contactdir,mapAtomTopology,firstTimeThrough,
     &   badAtom,category)

c----------------------------------------------
c  Each atom is assigned a category which describes the physico-chemical
c  properties of that atom
c
c 9-17 The three variable categorylist
c      has been changed to a two dimensional array to allow
c      each atom to be assigned two different categories:
c         categorylist(line,1) = 1 through 16  call these 'atom categories'
c         categorylist(line,2) = 1 through 5   call these 
c                    'environment categories'
c 12-23   categorylist(line,3) = 1 through 8 call these 'topology categories'
c
c      in the same vein, category(numatm) has been changed to allow
c      each atom to be descibed with two different environment designations.
c      The 2D array reflects this.
c      Numatoms_in_category will still reflect (i.e. count) both the numbers
c      of atoms binned via the first category scheme (atom categories) and
c      the second category scheme

      include 'parameters.h'
      include 'pdb.h'

      integer category(maxat,3)
      integer categorylist(maxdatalines,3)
      integer tempnumatomtypes1
      integer tempnumatomtypes2
      integer tempnumatomtypes3
      integer numatomtypes
      integer numatomtypesEnv
      integer numatomtypesTop
      integer mapAtomTopology(maxnumatomtypes)
      integer numunknown
   
      character*(*) contactdir
      character*(*) categorylistname
      character*120 concat
      character*7  junk
      character*4  atomtype(maxdatalines)
      character*3  restype(maxdatalines)
      logical notCategorized
      logical firstTimeThrough
      logical badAtom
      logical exist


13    format (a4,1x,a3,3x,i4,3x,i4,3x,i4)
14    format (a10)
c-----------------------------------------------
c  Initialize category
 
      do i=1,numatm
       category(i,1) = 0
       category(i,2) = 0
       category(i,3) = 0
      end do
      numunknown  = 0
      tempnumatomtypes1 = 1
      tempnumatomtypes2 = 1

      call fileexist(categorylistname,contactdir,concat,exist)
      if (.not.exist) then
       print *,'Error opening Atom/Category map file'
       print *,'Stopping'
       stop
      end if

      open (unit=1,file=concat,status='OLD')

      read (1,14) junk

      i=1
100   read (1,13,end = 200) atomtype(i),restype(i),categorylist(i,1),
     +   categorylist(i,2),categorylist(i,3)
      if (categorylist(i,1).gt.tempnumatomtypes1) then
       tempnumatomtypes1 = categorylist(i,1)
      end if
      if (categorylist(i,2).gt.tempnumatomtypes2) then
       tempnumatomtypes2 = categorylist(i,2)
      end if
      if (categorylist(i,3).gt.tempnumatomtypes3) then
       tempnumatomtypes3 = categorylist(i,3)
      end if

      i=i+1
      goto 100

200   close (1)

      numdatalines = i-1
      numatomtypes = tempnumatomtypes1
      numatomtypesEnv = tempnumatomtypes2
      numatomtypesTop = tempnumatomtypes3

      do i=1,numdatalines
       mapAtomTopology(categorylist(i,1)) = categorylist(i,3)
      end do

c now that I've read in all the atomtypes that I recognize,
c I'll attempt to label each atom in my pdbfile with an
c atomtype...

c If I fail with any atom, label that atom with an unknown flag

300   do i=1,numatm
       notCategorized = .true.
       unknown(i) = .true.
       do j=1,numdatalines
        if (atomname(i).eq.atomtype(j).and.res(i).eq.restype(j)) then
         notCategorized = .false.
         unknown(i) = .false.
         do k=1,3
          category(i,k) = categorylist(j,k)
         end do
        end if
       end do


       if (notCategorized) then
c         print *,'Not Categorized: ',i,' ',res(i),atomname(i)
         badAtom = .true.
       end if
      end do
      close (1)

c      do i=1,numatm
c       print *,i,'category:',category(i,1),category(i,2),category(i,3)
c      end do

      end

c------------------------------------------
      subroutine categorizeEncadAtoms (firstTimeThrough,badAtom,
     &   categorylistname,contactdir,category)

c----------------------------------------------
c  Each atom is assigned a category which describes the physico-chemical
c  properties of that atom
c
c 9-17 The three variable categorylist
c      has been changed to a two dimensional array to allow
c      each atom to be assigned two different categories:
c         categorylist(line,1) = 1 through 16  call these 'atom categories'
c         categorylist(line,2) = 1 through 5   call these 
c                    'environment categories'
c 12-23   categorylist(line,3) = 1 through 8 call these 'topology categories'
c
c      in the same vein, category(numatm) has been changed to allow
c      each atom to be descibed with two different environment designations.
c      The 2D array reflects this.
c      Numatoms_in_category will still reflect (i.e. count) both the numbers
c      of atoms binned via the first category scheme (atom categories) and
c      the second category scheme

      include 'parameters.h'
      include 'pdb.h'
      include 'encad.h'

      integer category(maxat,3)
      integer tempnumatomtypes1
      integer tempnumatomtypes2
      integer tempnumatomtypes3
      integer numatomtypes
      integer numatomtypesEnv
      integer numatomtypesTop
      integer mapAtomTopology(maxnumatomtypes)
      integer numunknown
   
      character*(*) contactdir
      character*(*) categorylistname
      character*120 concat
      character*7  junk
      character*3  convertFromSingleLetterCode
      logical compareIgnoringWhitespace
      logical notCategorized
      logical firstTimeThrough
      logical badAtom
      logical exist


13    format (a4,1x,a3,3x,i4,3x,i4,3x,i4)
14    format (a10)
c-----------------------------------------------
c  Initialize category
 
      do i=1,numatm
       category(i,1) = 0
       category(i,2) = 0
       category(i,3) = 0
      end do
      numunknown  = 0
      tempnumatomtypes1 = 1
      tempnumatomtypes2 = 1

      if (firstTimeThrough) then
        call fileexist(categorylistname,contactdir,concat,exist)
        if (.not.exist) then
         print *,'Error opening Atom/Category map file'
         print *,'Stopping'
         stop
        end if
	
	call read_encad_namelist()
	call generateCategorymap()
      end if

c now that I've read in all the atomtypes that I recognize,
c I'll attempt to label each atom in my pdbfile with an
c atomtype...

c If I fail with any atom, label that atom with an unknown flag

300   do i=1,numatm
       notCategorized = .true.
       unknown(i) = .true.
       do j=1,encadnumentries
        if (compareIgnoringWhitespace(atomname(i),encadatomname(j))
     &      .and.res(i).eq.convertFromSingleLetterCode(encadresname(j)))
     &       then
          notCategorized = .false.
          unknown(i) = .false.
	  category(i,1) = encadatomtypeint(j)
        end if
       end do

       if (notCategorized) then
         print *,'Not Categorized: ',i,' ',res(i),atomname(i)
         badAtom = .true.
       end if
      end do
      close (1)

c      do i=1,numatm
c       print *,i,'category:',category(i,1),category(i,2),category(i,3)
c      end do

      end
