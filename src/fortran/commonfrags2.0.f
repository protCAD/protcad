      program commonfrags

c Author: C. Summa
c-------------------------------------------------------
c Variable declarations

      include 'parameters.h'
      parameter (greatbignum = 1000000)
      include 'pdb.h'

c data from list
      integer i,k
      integer numfrags
      character*30 tempstring
      character*11 thepdbname(greatbignum)
      character*3 theres(greatbignum)
      integer theresnum(greatbignum)
      character*1 thechain(greatbignum)
      character*4 theatomname(greatbignum)
      character*6 thetype(greatbignum)
      integer thetimesseen(greatbignum)
      character*1 hit(greatbignum,10)
      real maxdis(10)
      integer skip(10)
      integer numfilescompared
      real criticaldistance
      integer skipnum

c environment variables
      character*80  listname
      character*80  contactdir
      character*80  contactpdbdir
      character*80  databasedir
      character*80  misfolddir

c searching for cluster in pdb file
      integer i1,i2
      integer bigcount
      integer clusternum
      logical relevant
      integer refAtomIndex
      character*11  pdbname
      real criticalmindistance
      real criticalmaxdistance
      real Distance
      integer maCounter
      integer markedAtoms(20)
      integer resNumsofMarkedAtoms(20)
      integer fragments(20,2)
      integer numfragments
      character*80 fragfilename
      character*5 stringcounter
      character*5 strincrement

11    format (a11)
12    format (i5, f5.1, i5)
13    format (a30)
14    format (1x,f3.1,1x,i1)
15    format (a11,3x,a3,2x,i4,2x,a1,1x,a4,2x,a6,2x,i3,
     &  2x,a1,a1,a1,a1,a1,a1,a1,a1,a1,a1,a1)

c------------------------------------------------------------
c Read in command line arguments

      numarguments = iargc()

      if (numarguments.ne.1) then
        print *,'Usage: commonfrags listname'
        goto 999
      endif

      call getarg(1,listname)
 
      criticalmindistance = 1.8

      open (unit=1,file=listname)
      i = 1
      j = 1
50    read (1,13,end=100) tempstring
      if (tempstring(1:3).ne.'pdb') then
        if (tempstring(2:2).ne.'A') then
         backspace(1)
         read (1,14) maxdis(j),skip(j)
         j = j + 1
        end if
        goto 50
      endif

      backspace (1) 

      read (1,15,end=100) thepdbname(i),theres(i),theresnum(i),
     &   thechain(i),theatomname(i),
     &   thetype(i),thetimesseen(i),(hit(i,k),k=1,10)
      i = i + 1
      goto 50


100   continue
      numfrags = i-1
      numfilescompared = j-1

c for debugging
c      do j=1,numfilescompared
c       print 14,maxdis(j),skip(j)
c      end do
c      do j=1,numfrags
c       print 15,thepdbname(j),theres(j),theresnum(j),
c     &   thechain(j),theatomname(j),
c     &   thetype(j),thetimesseen(j),(hit(j,k),k=1,10)
c      end do

      call getenv('CONTACTPDBDIR2',databasedir)

      stringcounter = '00001'
      do bigcount = 1,numfrags
c        if (thetimesseen(bigcount).ge.3) then
          call initialize_pdb_data()
          clusternum = bigcount
          call pdbread_no_biomt_reduced (thepdbname(clusternum),
     &       databasedir)
          call pick_distance_and_skip (hit,clusternum,
     &       maxdis,skip,criticalmaxdistance,numskip)
          print *,thepdbname(clusternum),criticalmaxdistance,numskip

c Now find reference atom around which cluster resides

          do i1 = 1,numatm
            if (thechain(clusternum).eq.chain(i1)) then
              if (theresnum(clusternum).eq.resnum(i1)) then
                if (theatomname(clusternum).eq.atomname(i1)) then
                  refAtomIndex=i1
                  exit
                end if
              end if
            end if
          end do
      
          do i1 = 1,20
           markedAtoms(i1) = 0
          end do

          maCounter = 0
          do i1 = 1,numatm
            call relevantDistance2(refAtomIndex,i1,numskip,relevant)
            if (relevant) then
              tempdistance = Distance(coord(1,refAtomIndex),
     &          coord(2,refAtomIndex),coord(3,refAtomIndex),
     &          coord(1,i1),coord(2,i1),coord(3,i1))
              if (tempdistance.le.criticalmaxdistance.and.
     &            tempdistance.gt.criticalmindistance) then 
                maCounter = maCounter + 1
                markedAtoms(maCounter) = i1
              endif
            endif
          end do

c label marked atoms with Bfactors

          do i1 = 1,numatm
            Bfac(i1) = 0.5
          end do
          do i1 = 1,maCounter
           Bfac(markedAtoms(maCounter)) = 1.0
          end do
          Bfac(refAtomIndex) = 0.0

          do i1 = 1,maCounter
           resNumsofMarkedAtoms(i1) = resnum(markedAtoms(i1))
          end do

          call find_lengths(resNumsofMarkedAtoms,maCounter,
     &       fragments,numfragments)

          do i1 = 1,numfragments
            print *,(fragments(i1,i2),i2=1,2)
          end do

          fragfilename = ''
          fragfilename(1:11) = thepdbname(clusternum)(1:11)
          fragfilename(12:12) = '_'
          fragfilename(13:17) = stringcounter
          
          call pdbwrite_fragments(fragfilename,fragments,
     &       numfragments,thechain(clusternum))
c        end if
       stringcounter = strincrement(stringcounter)
      end do

999   end

c---------------------------------------------------
      subroutine pick_distance_and_skip (hit,clusternum,
     &   maxdis,skip,criticalDistance,numskip)

      parameter (greatbignum = 1000000)
      real criticalDistance
      integer numskip
      integer clusternum
      character*1 hit(greatbignum,10)
      real maxdis(10)
      integer skip(10)

c     pick largest one - they should all be OK

      criticalDistance = 0.0
      numskip = 0
      do i=1,10
       if (hit(clusternum,i).eq.'x') then
         criticalDistance = maxdis(i)
         numskip = skip(i)
       end if
      end do
      end

c---------------------------------------------------

      subroutine find_lengths(resNumsofMarkedAtoms,
     &       numMarkedAtoms,fragments,numfragments)

      integer numfragments
      integer numMarkedAtoms
      integer fragments(20,2)
      integer resNumsofMarkedAtoms(20)
      logical firstfrag
      logical found

      do i=1,20
        do j=1,2
          fragments(i,j) = 0
        end do
      end do
      numfragments = 0
      firstfrag = .true.
      
      do i=1,numMarkedAtoms
        if (firstfrag) then
          numfragments = 1
          firstfrag = .false.
          fragments(numfragments,1) = resNumsofMarkedAtoms(i)-5
          fragments(numfragments,2) = resNumsofMarkedAtoms(i)+5
        end if
        found = .false.
        do j=1,numfragments
          if (resNumsofMarkedAtoms(i).ge.fragments(j,1).and.
     &        resNumsofMarkedAtoms(i).le.fragments(j,2)) then
            if (resNumsofMarkedAtoms(i).ge.
     &          resNumsofMarkedAtoms(j)) then
              fragments(j,2) = resNumsofMarkedAtoms(i) + 5
            else
              fragments(j,1) = resNumsofMarkedAtoms(i) - 5
            end if
            found = .true.
          end if
        end do
        if (.not.found) then
          numfragments = numfragments + 1
          fragments(numfragments,1) = resNumsofMarkedAtoms(i)-5
          fragments(numfragments,2) = resNumsofMarkedAtoms(i)+5
        end if
      end do
      end

c--------------------------------------------------------
      function strincrement(stringcounter)
     
      character*5 strincrement
      character*5 stringcounter

c  char(48) = '0'
c  char(57) = '9'

      if (ichar(stringcounter(5:5)).ne.57) then
        stringcounter(5:5) = char(ichar(stringcounter(5:5)) + 1)
        strincrement = stringcounter
        return
      end if

      stringcounter(5:5) = '0'

      if (ichar(stringcounter(4:4)).ne.57) then
        stringcounter(4:4) = char(ichar(stringcounter(4:4)) + 1)
        strincrement = stringcounter
        return
      end if

      stringcounter(4:4) = '0'

      if (ichar(stringcounter(3:3)).ne.57) then
        stringcounter(3:3) = char(ichar(stringcounter(3:3)) + 1)
        strincrement = stringcounter
        return
      end if

      stringcounter(3:3) = '0'

      if (ichar(stringcounter(2:2)).ne.57) then
        stringcounter(2:2) = char(ichar(stringcounter(2:2)) + 1)
        strincrement = stringcounter
        return
      end if

      stringcounter(2:2) = '0'

      if (ichar(stringcounter(1:1)).ne.57) then
        stringcounter(1:1) = char(ichar(stringcounter(1:1)) + 1)
        strincrement = stringcounter
        return
      end if

      print *,'Error - stringcounter greater than 9999!'
      stop
      end
      
