c------------------------------------------------------------
      subroutine pdbread (fname,dir)

      parameter (mxBMT = 50)

c----------------------------
c      character*80 cmnt(500)
c----------------------------

      include 'parameters.h'
      include 'pdb.h'
      character*80 pdbline1
      character*80 biomtline(3)

      character*(*) fname 
      integer fnstart
      integer fnend
      character*50 tempfname
      character*120 concat
      character*(*) dir

      integer i1,i2
      integer j1,j2,j3
      integer i,j,k,l,m
      integer t
      integer eOlig
      integer resinc
      integer srnum
      real    TFMMat(mxBMT,3,3)
      real    TLTMat(mxBMT,3)
      character BMTchain(mxBMT,20)
      character tempchain(mxBMT)
      logical BIOMTstart
      logical BIOMT
      logical art
      integer nTrans
      integer nChain(mxBMT)
      integer tempNchain
 
      integer atcnt

      character*6 tempah
      character*4 tempanam
      character*3 tempres
      character*1 tempsub
      character*1 tempc
      character*1 tempmc
      integer     tempanum
      integer     temptru
      real        tempx
      real        tempy
      real        tempz
      real        temprad
      real        tempB  

      integer frnum,nrnum,lrnum
      character*1 lsnum
      integer lanum

      integer nat
      integer origNumChains
      integer newatplace
      real    oldCoord(3)
      real    newCoord(3)
      character IncChar
      real tempMatA(3,3)
      real tempMatB(3)
      logical exist

      integer StringLength
      integer linelen
 

      resinc = 0
      srnum = 0
      do j1 = 1,mxBMT
       do j2 = 1,3
        do j3 = 1,3
         TFMMat(j1,j2,j3) = 0.0
        end do
        TLTMat(j1,j2) = 0.0
       end do
       nChain(j1) = 0
      end do
      nTrans = 0
      atcnt = 0
      art = .false.
      numchains = 1
      BIOMT = .false.

10    format (a80)

      call fileexist(fname,dir,concat,exist)
      if (.not.exist) then
       goto 995
      end if

c      print *,'concat = ',concat

      i1 = 1
      open (unit=1,err = 995,file=concat,status='OLD')

100   read (1,10,end=900) pdbline1
c      print *,pdbline1

      call catEntry(pdbline1,t)

c      print *,'type = ',t

c OK, now if it falls into categories 1 through 9,
c IGNORE IT and add it to the cmnt array

      if ((t.ge.1.and.t.le.9).or.
     &    (t.eq.28).or.
     &    (t.eq.29).or.
     &    (t.eq.30).or.
     &    (t.eq.31).or.
     &    (t.eq.32).or.
     &    (t.eq.33).or.
     &    (t.eq.34).or.
     &    (t.eq.35).or.
     &    (t.eq.36).or.
     &    (t.eq.37).or.
     &    (t.eq.38).or.
     &    (t.eq.999)) then
       cmnt(i1) = pdbline1
       i1 = i1 + 1
       goto 100
      end if

      if (t.eq.39) then
c Ignore it, it's an 'ANISOU' line
       goto 100
      end if

c If t = 10, then it's a remark... 
c There are currenly two pieces of data that we want
c to retrieve from remarks:
c
c	words such as DIMER, TETRAMER, etc.
c	BIOMT data

      if (t.eq.10) then
       if (pdbline1(8:10).eq.'350') then
        if (pdbline1(12:16).eq.'APPLY') then
         call startBIOMT(pdbline1,BIOMTstart,tempNchain,
     &      tempchain,art)
         cmnt(i1) = pdbline1
         i1 = i1 + 1
        else
         BIOMTstart = .false.  
        end if
        if (BIOMTstart) then
c Now, immediately read in the next 3 lines and send this
c data to a special parser surbroutine
         nTrans = nTrans + 1
         nChain(nTrans) = tempNchain
         do i2 = 1,nChain(nTrans)
          BMTchain(nTrans,i2) = tempchain(i2)
         end do
         read (1,10) biomtline(1)
         cmnt(i1) = biomtline(1)
         i1 = i1 + 1
         read (1,10) biomtline(2)
         cmnt(i1) = biomtline(2)
         i1 = i1 + 1
         read (1,10) biomtline(3)
         cmnt(i1) = biomtline(3)
         i1 = i1 + 1
         call parseBIOMT(nTrans,biomtline,TFMMat,TLTMat,art)
         BIOMT = .true.
        else
         if (BIOMT.and.(StringLength(pdbline1).gt.15)) then
          biomtline(1) = pdbline1
          cmnt(i1) = biomtline(1)  
          i1 = i1 + 1
          read (1,10) biomtline(2)
          cmnt(i1) = biomtline(2)
          i1 = i1 + 1
          read (1,10) biomtline(3)
          cmnt(i1) = biomtline(3)
          i1 = i1 + 1
          linelen = StringLength(biomtline(1))
          if (linelen.gt.15.and.
     &        biomtline(1)(14:18).eq.'BIOMT') then
           nTrans = nTrans + 1
           nChain(nTrans) = nChain(nTrans-1)
           do i2 = 1,nChain(nTrans)
            BMTchain(nTrans,i2) = BMTchain(nTrans-1,i2)
           end do
           call parseBIOMT(nTrans,biomtline,TFMMat,TLTMat,art)
          end if
          BIOMT = .true.
         end if
        end if
       end if
       if (pdbline1(8:10).eq.'300') then
        call oligomercheck(pdbline1,eOlig)
        cmnt(i1) = pdbline1
        i1 = i1 + 1
       end if
       goto 100
      end if

c Again, if its one of the following types, we want to simply
c store it without manipulation
      if (t.gt.10.and.t.le.21) then
       cmnt(i1) = pdbline1
       i1 = i1 + 1
       goto 100
      end if
 
      if (t.eq.22.or.t.eq.24) then

c It's an ATOM or HETATM record - parse the line and 
c extract the data

       call parseATOM(pdbline1,tempah,tempanum,tempanam,
     &  tempmc,tempres,tempc,temptru,tempsub,tempx,tempy,
     &  tempz,temprad,tempB)

c OK, now run some tests to decide if we're going to pay
c attention to this particular atom entry

c Ignore Water
c       if (tempres.eq.'H0H'.or.
c     &     tempres.eq.'D0D'.or.
c     &     tempres.eq.'S04'.or.
c     &     tempres.eq.'HOH'.or.
c     &     tempres.eq.'DOD'.or.
c     &     tempres.eq.'SO4') then
c        goto 100
c       end if  

c Ignore Multiple entries of a given atom
c       if (atcnt.gt.1.and.
c     &    tempres.eq.resnum(atcnt).and.
c     &    tempanam.eq.atomname(atcnt).and.
c     &    tempsub.ne.subnum(atcnt)) then
c        goto 100
c       end if 
 
c Ignore Hydrogen
       if (tempanam(2:2).eq.'H') then
        goto 100
       end if
      
       if (atcnt.eq.maxat-1) then
         print *,'Too Many Atoms for pdbreader!'
         print *,'It''s probably an NMR structure'
         stop
       end if

       if (atcnt.eq.0) then
        lrnum = temptru
        frnum = temptru
        lsnum = tempsub
        nrnum = lrnum + 1
       end if
       atcnt = atcnt + 1
       atheader(atcnt) = tempah
       atomnum(atcnt)  = tempanum
       atomname(atcnt) = tempanam
       multiconf(atcnt) = tempmc
       res(atcnt) = tempres
       if (t.eq.24) then
        hetatm(atcnt) = .true.
       end if
       if (tempc.eq.' '.and.art) then
        chain(atcnt) = 'A'
       else
        chain(atcnt) = tempc
       end if
       trueResnum(atcnt) = temptru
       subnum(atcnt) = tempsub
       coord(1,atcnt) = tempx
       coord(2,atcnt) = tempy
       coord(3,atcnt) = tempz
       rad(atcnt) = temprad
       Bfac(atcnt) = tempB

c       print *,'Internal Representation'
c       print *,atheader(atcnt),atomnum(atcnt),atomname(atcnt),
c     &   multiconf(atcnt),res(atcnt),chain(atcnt),
c     &   trueResnum(atcnt),subnum(atcnt)

      end if

      if (t.eq.23) then
c it's a terminator line - do nothing!
       goto 100
      end if

      if (t.eq.25) then
c ignore it (it's a master line)
       goto 100
      end if

      if (t.eq.26) then
c ignore it (it's connectivity data)
       goto 100
      end if

      if (t.eq.27) then
c we're at the end of the pdbfile 
       goto 100
      end if

      if (atcnt.gt.0) then
       if (atcnt.gt.1) then
        if (chain(atcnt).ne.chain(atcnt-1)) then
          numchains = numchains + 1
        end if
       end if 

      if (trueResnum(atcnt).eq.lrnum.and.
     &    subnum(atcnt).eq.lsnum) then
       goto 300
      end if

      if (trueResnum(atcnt).eq.lrnum.and.
     &    subnum(atcnt).ne.lsnum) then
       numres = numres + 1
       resinc = resinc + 1
       goto 300
      end if

      if (trueResnum(atcnt).eq.nrnum.and.
     &    subnum(atcnt).eq.' ') then
       numres = numres + 1
       nrnum = nrnum + 1
       goto 300
      endif

      if (trueResnum(atcnt).gt.nrnum.and.
     &    subnum(atcnt).eq.' ') then
       resinc = resinc +
     &  (trueResnum(atcnt)-nrnum-1)
       numres = numres + 1
       nrnum = nrnum + resinc
       goto 300 
      endif

      if (trueResnum(atcnt).lt.lrnum) then
       resinc = 0
       numres = numres + 1
       frnum = trueResnum(atcnt)
       nrnum = trueResnum(atcnt) +1
       lrnum = trueResnum(atcnt)
       goto 300
      end if
      end if

300   if (atcnt.gt.0) then
      lrnum = trueResnum(atcnt)
      lsnum = subnum(atcnt)
      lanum = atomnum(atcnt)
      resnum(atcnt) = trueResnum(atcnt) +
     &  resinc - srnum
      subnum(atcnt) = ' '
      end if

c      print *,trueResnum(atcnt),resinc,srnum
      goto 100

900   close (1)
      numatm = atcnt 
      numcomlines = i1-1

c      print *,'nTrans = ',nTrans
      if (numatm*nTrans.gt.maxat) then
       print *,'transformations will produce too many atoms!
     &   increase maxat!'
       reject = .true.
       goto 999
      end if   

      origNumchains = numchains

      if (BIOMT) then
       print *,fname,' is a multimer'
       do i=1,nTrans
        do j=1,nChain(i)
         nat = 0
         lanum = atomnum(numatm)
         lrnum = resnum(numatm)
         do k=1,numatm
          if (chain(k).eq.BMTchain(i,j)) then
           nat = nat + 1
           newatplace = numatm + nat

           do l=1,3
            oldCoord(l) = coord(l,k)
            tempMatB(l) = TLTMat(i,l)
            do m=1,3
             tempMatA(l,m) = TFMMat(i,l,m)
            end do
           end do


           call transform(tempMatA,tempMatB,oldCoord,newCoord)

           coord(1,newatplace) = newCoord(1)
           coord(2,newatplace) = newCoord(2)
           coord(3,newatplace) = newCoord(3)

           atheader(newatplace) = atheader(k)
           atomname(newatplace) = atomname(k)
           multiconf(newatplace) = multiconf(k)
           res(newatplace) = res(k)
           subnum(newatplace) = subnum(k)

           chain(newatplace) = incChar(chain(k),numchains)

           atomnum(newatplace) = lanum + nat
           trueResnum(newatplace) = trueResnum(k)

           Bfac(newatplace) = Bfac(k)
           rad(newatplace) = rad(k)
           hetatm(newatplace) = hetatm(k)
                                         
          end if
         end do
         numchains = numchains + 1
         numatm = numatm + nat
        end do
       end do
      end if
      do i=1,numatm
       if (art.and.chain(i).eq.'A') then
        chain(i) = ' '
       end if
      end do
      if (art) then
       print *,'Artificial designation and removal of chain specifier'
      end if
      
      call findNextWord(fname,1,fnstart,fnend,exist)
      if (exist) then
       tempfname= 'mod.'//fname(fnstart:fnend)
      end if
      goto 997

995   print *,'ERROR opening pdbfile!'
      stop


997   continue
999   end

c------------------------------------------------------------
      subroutine pdbread_no_biomt (fname,dir)

      parameter (mxBMT = 50)

c----------------------------
c      character*80 cmnt(500)
c----------------------------

      include 'parameters.h'
      include 'pdb.h'
      character*80 pdbline1
      character*80 biomtline(3)

      character*(*) fname 
      integer fnstart
      integer fnend
      character*50 tempfname
      character*120 concat
      character*(*) dir

      integer i1,i2
      integer j1,j2,j3
      integer i,j,k,l,m
      integer t
      integer eOlig
      integer resinc
      integer srnum
      real    TFMMat(mxBMT,3,3)
      real    TLTMat(mxBMT,3)
      character BMTchain(mxBMT,20)
      character tempchain(mxBMT)
      logical BIOMTstart
      logical BIOMT
      logical art
      integer nTrans
      integer nChain(mxBMT)
      integer tempNchain
 
      integer atcnt

      character*6 tempah
      character*4 tempanam
      character*3 tempres
      character*1 tempsub
      character*1 tempc
      character*1 tempmc
      integer     tempanum
      integer     temptru
      real        tempx
      real        tempy
      real        tempz
      real        temprad
      real        tempB  

      integer frnum,nrnum,lrnum
      character*1 lsnum
      integer lanum

      integer nat
      integer origNumChains
      integer newatplace
      real    oldCoord(3)
      real    newCoord(3)
      character IncChar
      real tempMatA(3,3)
      real tempMatB(3)
      logical exist
      character*9 tempfname1
      character*11 tempfname2

      integer StringLength
      integer linelen

      resinc = 0
      srnum = 0
      do j1 = 1,mxBMT
       do j2 = 1,3
        do j3 = 1,3
         TFMMat(j1,j2,j3) = 0.0
        end do
        TLTMat(j1,j2) = 0.0
       end do
       nChain(j1) = 0
      end do
      nTrans = 0
      atcnt = 0
      art = .false.
      numchains = 1
      BIOMT = .false.

10    format (a80)

      call fileexist(fname,dir,concat,exist)
c      print *,'fname = ',fname
c      print *,'exist = ',exist
      if (.not.exist) then
       tempfname1=''
       tempfname1(1:4) = fname(4:7)
       tempfname1(5:9) = '.mmol'       
       call fileexist(tempfname1,dir,concat,exist)
c       print *,'tempfname1 = ',tempfname1
c       print *,'exist = ',exist
       if (.not.exist) then
        tempfname2 = ''
        tempfname2(1:4) = fname(4:7)
        tempfname2(5:5) = '_'
        tempfname2(6:6) = '1'
        tempfname2(7:11) = '.mmol'
        call fileexist(tempfname2,dir,concat,exist)
c        print *,'tempfname2 = ',tempfname2
c        print *,'exist = ',exist
        if (.not.exist) then
          goto 995
        end if
       end if
      end if

c      print *,'concat = ',concat

      i1 = 1
      open (unit=1,err = 995,file=concat,status='OLD')

100   read (1,10,end=900) pdbline1
c      print *,pdbline1

      call catEntry(pdbline1,t)

c      print *,'type = ',t

c OK, now if it falls into categories 1 through 9,
c IGNORE IT and add it to the cmnt array

      if ((t.ge.1.and.t.le.9).or.
     &    (t.eq.28).or.
     &    (t.eq.29).or.
     &    (t.eq.30).or.
     &    (t.eq.31).or.
     &    (t.eq.32).or.
     &    (t.eq.33).or.
     &    (t.eq.34).or.
     &    (t.eq.35).or.
     &    (t.eq.36).or.
     &    (t.eq.37).or.
     &    (t.eq.38).or.
     &    (t.eq.999)) then
       cmnt(i1) = pdbline1
       i1 = i1 + 1
       goto 100
      end if

      if (t.eq.39) then
c Ignore it, it's an 'ANISOU' line
       goto 100
      end if

c If t = 10, then it's a remark... 
c There are currenly two pieces of data that we want
c to retrieve from remarks:
c
c	words such as DIMER, TETRAMER, etc.
c	BIOMT data

      if (t.eq.10) then
       if (pdbline1(8:10).eq.'350') then
        if (pdbline1(12:16).eq.'APPLY') then
         call startBIOMT(pdbline1,BIOMTstart,tempNchain,
     &      tempchain,art)
         cmnt(i1) = pdbline1
         i1 = i1 + 1
        else
         BIOMTstart = .false.  
        end if
        if (BIOMTstart) then
c Now, immediately read in the next 3 lines and send this
c data to a special parser surbroutine
         nTrans = nTrans + 1
         nChain(nTrans) = tempNchain
         do i2 = 1,nChain(nTrans)
          BMTchain(nTrans,i2) = tempchain(i2)
         end do
         read (1,10) biomtline(1)
         cmnt(i1) = biomtline(1)
         i1 = i1 + 1
         read (1,10) biomtline(2)
         cmnt(i1) = biomtline(2)
         i1 = i1 + 1
         read (1,10) biomtline(3)
         cmnt(i1) = biomtline(3)
         i1 = i1 + 1
         call parseBIOMT(nTrans,biomtline,TFMMat,TLTMat,art)
         BIOMT = .true.
        else
         if (BIOMT.and.(StringLength(pdbline1).gt.15)) then
          biomtline(1) = pdbline1
          cmnt(i1) = biomtline(1)  
          i1 = i1 + 1
          read (1,10) biomtline(2)
          cmnt(i1) = biomtline(2)
          i1 = i1 + 1
          read (1,10) biomtline(3)
          cmnt(i1) = biomtline(3)
          i1 = i1 + 1
          linelen = StringLength(biomtline(1))
          if (linelen.gt.15.and.
     &        biomtline(1)(14:18).eq.'BIOMT') then
           nTrans = nTrans + 1
           nChain(nTrans) = nChain(nTrans-1)
           do i2 = 1,nChain(nTrans)
            BMTchain(nTrans,i2) = BMTchain(nTrans-1,i2)
           end do
           call parseBIOMT(nTrans,biomtline,TFMMat,TLTMat,art)
          end if
          BIOMT = .true.
         end if
        end if
       end if
       if (pdbline1(8:10).eq.'300') then
        call oligomercheck(pdbline1,eOlig)
        cmnt(i1) = pdbline1
        i1 = i1 + 1
       end if
       goto 100
      end if

c Again, if its one of the following types, we want to simply
c store it without manipulation
      if (t.gt.10.and.t.le.21) then
       cmnt(i1) = pdbline1
       i1 = i1 + 1
       goto 100
      end if
 
      if (t.eq.22.or.t.eq.24) then

c It's an ATOM or HETATM record - parse the line and 
c extract the data

       call parseATOM(pdbline1,tempah,tempanum,tempanam,
     &  tempmc,tempres,tempc,temptru,tempsub,tempx,tempy,
     &  tempz,temprad,tempB)

c OK, now run some tests to decide if we're going to pay
c attention to this particular atom entry

c Ignore Water
c       if (tempres.eq.'H0H'.or.
c     &     tempres.eq.'D0D'.or.
c     &     tempres.eq.'S04'.or.
c     &     tempres.eq.'HOH'.or.
c     &     tempres.eq.'DOD'.or.
c     &     tempres.eq.'SO4') then
c        goto 100
c       end if  

c Ignore Multiple entries of a given atom
c       if (atcnt.gt.1.and.
c     &    tempres.eq.resnum(atcnt).and.
c     &    tempanam.eq.atomname(atcnt).and.
c     &    tempsub.ne.subnum(atcnt)) then
c        goto 100
c       end if 
 
c Ignore Hydrogen
       if (tempanam(2:2).eq.'H') then
        goto 100
       end if
      
       if (atcnt.eq.maxat-1) then
         print *,'Too Many Atoms for pdbreader!'
         print *,'It''s probably an NMR structure'
         stop
       end if

       if (atcnt.eq.0) then
        lrnum = temptru
        frnum = temptru
        lsnum = tempsub
        nrnum = lrnum + 1
       end if
       atcnt = atcnt + 1
       atheader(atcnt) = tempah
       atomnum(atcnt)  = tempanum
       atomname(atcnt) = tempanam
       multiconf(atcnt) = tempmc
       res(atcnt) = tempres
       if (t.eq.24) then
        hetatm(atcnt) = .true.
       end if
       if (tempc.eq.' '.and.art) then
        chain(atcnt) = 'A'
       else
        chain(atcnt) = tempc
       end if
       trueResnum(atcnt) = temptru
       subnum(atcnt) = tempsub
       coord(1,atcnt) = tempx
       coord(2,atcnt) = tempy
       coord(3,atcnt) = tempz
       rad(atcnt) = temprad
       Bfac(atcnt) = tempB

c       print *,'Internal Representation'
c       print *,atheader(atcnt),atomnum(atcnt),atomname(atcnt),
c     &   multiconf(atcnt),res(atcnt),chain(atcnt),
c     &   trueResnum(atcnt),subnum(atcnt)

      end if

      if (t.eq.23) then
c it's a terminator line - do nothing!
       goto 100
      end if

      if (t.eq.25) then
c ignore it (it's a master line)
       goto 100
      end if

      if (t.eq.26) then
c ignore it (it's connectivity data)
       goto 100
      end if

      if (t.eq.27) then
c we're at the end of the pdbfile 
       goto 100
      end if

      if (atcnt.gt.0) then
       if (atcnt.gt.1) then
        if (chain(atcnt).ne.chain(atcnt-1)) then
          numchains = numchains + 1
        end if
       end if 

      if (trueResnum(atcnt).eq.lrnum.and.
     &    subnum(atcnt).eq.lsnum) then
       goto 300
      end if

      if (trueResnum(atcnt).eq.lrnum.and.
     &    subnum(atcnt).ne.lsnum) then
       numres = numres + 1
       resinc = resinc + 1
       goto 300
      end if

      if (trueResnum(atcnt).eq.nrnum.and.
     &    subnum(atcnt).eq.' ') then
       numres = numres + 1
       nrnum = nrnum + 1
       goto 300
      endif

      if (trueResnum(atcnt).gt.nrnum.and.
     &    subnum(atcnt).eq.' ') then
       resinc = resinc +
     &  (trueResnum(atcnt)-nrnum-1)
       numres = numres + 1
       nrnum = nrnum + resinc
       goto 300 
      endif

      if (trueResnum(atcnt).lt.lrnum) then
       resinc = 0
       numres = numres + 1
       frnum = trueResnum(atcnt)
       nrnum = trueResnum(atcnt) +1
       lrnum = trueResnum(atcnt)
       goto 300
      end if
      end if

300   if (atcnt.gt.0) then
      lrnum = trueResnum(atcnt)
      lsnum = subnum(atcnt)
      lanum = atomnum(atcnt)
      resnum(atcnt) = trueResnum(atcnt) +
     &  resinc - srnum
      subnum(atcnt) = ' '
      end if

      goto 100

900   close (1)
      numatm = atcnt 
      numcomlines = i1-1

      origNumchains = numchains
      goto 997

995   print *,'ERROR opening pdbfile!'
      stop


997   continue
999   end

c------------------------------------------------------------
      subroutine pdbread_no_biomt_reduced (fname,dir)

      parameter (mxBMT = 50)

c----------------------------
c      character*80 cmnt(500)
c----------------------------

      include 'parameters.h'
      include 'pdb.h'
      character*80 pdbline1
      character*80 biomtline(3)

      character*(*) fname 
      integer fnstart
      integer fnend
      character*50 tempfname
      character*120 concat
      character*(*) dir
      character*80 fname2

      integer i1,i2
      integer j1,j2,j3
      integer i,j,k,l,m
      integer t
      integer eOlig
      integer resinc
      integer srnum
      real    TFMMat(mxBMT,3,3)
      real    TLTMat(mxBMT,3)
      character BMTchain(mxBMT,20)
      character tempchain(mxBMT)
      logical BIOMTstart
      logical BIOMT
      logical art
      integer nTrans
      integer nChain(mxBMT)
      integer tempNchain
 
      integer atcnt

      character*6 tempah
      character*4 tempanam
      character*3 tempres
      character*1 tempsub
      character*1 tempc
      character*1 tempmc
      integer     tempanum
      integer     temptru
      real        tempx
      real        tempy
      real        tempz
      real        temprad
      real        tempB  

      integer frnum,nrnum,lrnum
      character*1 lsnum
      integer lanum

      integer nat
      integer origNumChains
      integer newatplace
      real    oldCoord(3)
      real    newCoord(3)
      character IncChar
      real tempMatA(3,3)
      real tempMatB(3)
      logical exist
      character*80 tempfname1
      character*80 tempfname2

      integer StringLength
      integer linelen

      resinc = 0
      srnum = 0
      do j1 = 1,mxBMT
       do j2 = 1,3
        do j3 = 1,3
         TFMMat(j1,j2,j3) = 0.0
        end do
        TLTMat(j1,j2) = 0.0
       end do
       nChain(j1) = 0
      end do
      nTrans = 0
      atcnt = 0
      art = .false.
      numchains = 1
      BIOMT = .false.

10    format (a80)

      print *,dir

c  first see if the output from Jane Richardon's "reduce"
c  exists in the directory
      fname2 = ''
      fname2(1:11) = fname(1:11)
      fname2(12:32) = '_reduced_no_hydrogens'
      call fileexist(fname2,dir,concat,exist)
c      print *,'fname = ',fname2
c      print *,'exist = ',exist
      if (.not.exist) then
       tempfname1=''
       tempfname1(1:4) = fname(4:7)
       tempfname1(5:30) = '.mmol_reduced_no_hydrogens'       
       call fileexist(tempfname1,dir,concat,exist)
c       print *,'tempfname1 = ',tempfname1
c       print *,'exist = ',exist
       if (.not.exist) then
        tempfname2 = ''
        tempfname2(1:4) = fname(4:7)
        tempfname2(5:5) = '_'
        tempfname2(6:6) = '1'
        tempfname2(7:32) = '.mmol_reduced_no_hydrogens'
        call fileexist(tempfname2,dir,concat,exist)
c        print *,'tempfname2 = ',tempfname2
c        print *,'exist = ',exist
        if (.not.exist) then
          goto 995
        end if
       end if
      end if

c      print *,'concat = ',concat

      i1 = 1
      open (unit=1,err = 995,file=concat,status='OLD')

100   read (1,10,end=900) pdbline1
c      print *,pdbline1

      call catEntry(pdbline1,t)

c      print *,'type = ',t

c OK, now if it falls into categories 1 through 9,
c IGNORE IT and add it to the cmnt array

      if ((t.ge.1.and.t.le.9).or.
     &    (t.eq.28).or.
     &    (t.eq.29).or.
     &    (t.eq.30).or.
     &    (t.eq.31).or.
     &    (t.eq.32).or.
     &    (t.eq.33).or.
     &    (t.eq.34).or.
     &    (t.eq.35).or.
     &    (t.eq.36).or.
     &    (t.eq.37).or.
     &    (t.eq.38).or.
     &    (t.eq.39).or.
     &    (t.eq.40).or.
     &    (t.eq.999)) then
       cmnt(i1) = pdbline1
       i1 = i1 + 1
       goto 100
      end if

c If t = 10, then it's a remark... 
c There are currenly two pieces of data that we want
c to retrieve from remarks:
c
c	words such as DIMER, TETRAMER, etc.
c	BIOMT data

      if (t.eq.10) then
       if (pdbline1(8:10).eq.'350') then
        if (pdbline1(12:16).eq.'APPLY') then
         call startBIOMT(pdbline1,BIOMTstart,tempNchain,
     &      tempchain,art)
         cmnt(i1) = pdbline1
         i1 = i1 + 1
        else
         BIOMTstart = .false.  
        end if
        if (BIOMTstart) then
c Now, immediately read in the next 3 lines and send this
c data to a special parser surbroutine
         nTrans = nTrans + 1
         nChain(nTrans) = tempNchain
         do i2 = 1,nChain(nTrans)
          BMTchain(nTrans,i2) = tempchain(i2)
         end do
         read (1,10) biomtline(1)
         cmnt(i1) = biomtline(1)
         i1 = i1 + 1
         read (1,10) biomtline(2)
         cmnt(i1) = biomtline(2)
         i1 = i1 + 1
         read (1,10) biomtline(3)
         cmnt(i1) = biomtline(3)
         i1 = i1 + 1
         call parseBIOMT(nTrans,biomtline,TFMMat,TLTMat,art)
         BIOMT = .true.
        else
         if (BIOMT.and.(StringLength(pdbline1).gt.15)) then
          biomtline(1) = pdbline1
          cmnt(i1) = biomtline(1)  
          i1 = i1 + 1
          read (1,10) biomtline(2)
          cmnt(i1) = biomtline(2)
          i1 = i1 + 1
          read (1,10) biomtline(3)
          cmnt(i1) = biomtline(3)
          i1 = i1 + 1
          linelen = StringLength(biomtline(1))
          if (linelen.gt.15.and.
     &        biomtline(1)(14:18).eq.'BIOMT') then
           nTrans = nTrans + 1
           nChain(nTrans) = nChain(nTrans-1)
           do i2 = 1,nChain(nTrans)
            BMTchain(nTrans,i2) = BMTchain(nTrans-1,i2)
           end do
           call parseBIOMT(nTrans,biomtline,TFMMat,TLTMat,art)
          end if
          BIOMT = .true.
         end if
        end if
       end if
       if (pdbline1(8:10).eq.'300') then
        call oligomercheck(pdbline1,eOlig)
        cmnt(i1) = pdbline1
        i1 = i1 + 1
       end if
       goto 100
      end if

c Again, if its one of the following types, we want to simply
c store it without manipulation
      if (t.gt.10.and.t.le.21) then
       cmnt(i1) = pdbline1
       i1 = i1 + 1
       goto 100
      end if
 
      if (t.eq.22.or.t.eq.24) then

c It's an ATOM or HETATM record - parse the line and 
c extract the data

       call parseATOM(pdbline1,tempah,tempanum,tempanam,
     &  tempmc,tempres,tempc,temptru,tempsub,tempx,tempy,
     &  tempz,temprad,tempB)

c OK, now run some tests to decide if we're going to pay
c attention to this particular atom entry

c Ignore Water
c       if (tempres.eq.'H0H'.or.
c     &     tempres.eq.'D0D'.or.
c     &     tempres.eq.'S04'.or.
c     &     tempres.eq.'HOH'.or.
c     &     tempres.eq.'DOD'.or.
c     &     tempres.eq.'SO4') then
c        goto 100
c       end if  

c Ignore Multiple entries of a given atom
c       if (atcnt.gt.1.and.
c     &    tempres.eq.resnum(atcnt).and.
c     &    tempanam.eq.atomname(atcnt).and.
c     &    tempsub.ne.subnum(atcnt)) then
c        goto 100
c       end if 
 
c Ignore Hydrogen
       if (tempanam(2:2).eq.'H') then
        goto 100
       end if
      
       if (atcnt.eq.maxat-1) then
         print *,'Too Many Atoms for pdbreader!'
         print *,'It''s probably an NMR structure'
         stop
       end if

       if (atcnt.eq.0) then
        lrnum = temptru
        frnum = temptru
        lsnum = tempsub
        nrnum = lrnum + 1
       end if
       atcnt = atcnt + 1
       atheader(atcnt) = tempah
       atomnum(atcnt)  = tempanum
       atomname(atcnt) = tempanam
       multiconf(atcnt) = tempmc
       res(atcnt) = tempres
       if (t.eq.24) then
        hetatm(atcnt) = .true.
       end if
       if (tempc.eq.' '.and.art) then
        chain(atcnt) = 'A'
       else
        chain(atcnt) = tempc
       end if
       trueResnum(atcnt) = temptru
       subnum(atcnt) = tempsub
       coord(1,atcnt) = tempx
       coord(2,atcnt) = tempy
       coord(3,atcnt) = tempz
       rad(atcnt) = temprad
       Bfac(atcnt) = tempB

c       print *,'Internal Representation'
c       print *,atheader(atcnt),atomnum(atcnt),atomname(atcnt),
c     &   multiconf(atcnt),res(atcnt),chain(atcnt),
c     &   trueResnum(atcnt),subnum(atcnt)

      end if

      if (t.eq.23) then
c it's a terminator line - do nothing!
       goto 100
      end if

      if (t.eq.25) then
c ignore it (it's a master line)
       goto 100
      end if

      if (t.eq.26) then
c ignore it (it's connectivity data)
       goto 100
      end if

      if (t.eq.27) then
c we're at the end of the pdbfile 
       goto 100
      end if

      if (atcnt.gt.0) then
       if (atcnt.gt.1) then
        if (chain(atcnt).ne.chain(atcnt-1)) then
          numchains = numchains + 1
        end if
       end if 

      if (trueResnum(atcnt).eq.lrnum.and.
     &    subnum(atcnt).eq.lsnum) then
       goto 300
      end if

      if (trueResnum(atcnt).eq.lrnum.and.
     &    subnum(atcnt).ne.lsnum) then
       numres = numres + 1
       resinc = resinc + 1
       goto 300
      end if

      if (trueResnum(atcnt).eq.nrnum.and.
     &    subnum(atcnt).eq.' ') then
       numres = numres + 1
       nrnum = nrnum + 1
       goto 300
      endif

      if (trueResnum(atcnt).gt.nrnum.and.
     &    subnum(atcnt).eq.' ') then
       resinc = resinc +
     &  (trueResnum(atcnt)-nrnum-1)
       numres = numres + 1
       nrnum = nrnum + resinc
       goto 300 
      endif

      if (trueResnum(atcnt).lt.lrnum) then
       resinc = 0
       numres = numres + 1
       frnum = trueResnum(atcnt)
       nrnum = trueResnum(atcnt) +1
       lrnum = trueResnum(atcnt)
       goto 300
      end if
      end if

300   if (atcnt.gt.0) then
      lrnum = trueResnum(atcnt)
      lsnum = subnum(atcnt)
      lanum = atomnum(atcnt)
      resnum(atcnt) = trueResnum(atcnt) +
     &  resinc - srnum
      subnum(atcnt) = ' '
      end if

      goto 100

900   close (1)
      numatm = atcnt 
      numcomlines = i1-1

      origNumchains = numchains
      goto 997

995   print *,'ERROR opening pdbfile!'
      stop


997   continue
999   end

c------------------------------------------------------
c------------------------------------------------------
      subroutine pdbread_w_write (fname,dir)

      parameter (mxBMT = 50)

      include 'parameters.h'
      include 'pdb.h'

      character*80 pdbline1
      character*80 biomtline(3)

      character*(*) fname 
      integer fnstart
      integer fnend
      character*50 tempfname
      character*120 concat
      character*(*) dir

      integer i1,i2
      integer j1,j2,j3
      integer i,j,k,l,m
      integer t
      integer eOlig
      integer resinc
      integer srnum
      real    TFMMat(mxBMT,3,3)
      real    TLTMat(mxBMT,3)
      character BMTchain(mxBMT,20)
      character tempchain(mxBMT)
      logical BIOMTstart
      logical BIOMT
      logical art
      integer nTrans
      integer nChain(mxBMT)
      integer tempNchain
 
      integer atcnt

      character*6 tempah
      character*4 tempanam
      character*3 tempres
      character*1 tempsub
      character*1 tempc
      character*1 tempmc
      integer     tempanum
      integer     temptru
      real        tempx
      real        tempy
      real        tempz
      real        temprad
      real        tempB  

      integer frnum,nrnum,lrnum
      character*1 lsnum
      integer lanum

      integer nat
      integer origNumChains
      integer newatplace
      real    oldCoord(3)
      real    newCoord(3)
      character IncChar
      real tempMatA(3,3)
      real tempMatB(3)
      logical exist

      integer StringLength
      integer linelen
 

      resinc = 0
      srnum = 0
      do j1 = 1,mxBMT
       do j2 = 1,3
        do j3 = 1,3
         TFMMat(j1,j2,j3) = 0.0
        end do
        TLTMat(j1,j2) = 0.0
       end do
       nChain(j1) = 0
      end do
      nTrans = 0
      atcnt = 0
      art = .false.
      numchains = 1
      BIOMT = .false.

10    format (a80)

      call fileexist(fname,dir,concat,exist)
      if (.not.exist) then
       goto 995
      end if

      i1 = 1
      open (unit=1,err = 995,file=concat,status='OLD')

100   read (1,10,end=900) pdbline1

      call catEntry(pdbline1,t)


c OK, now if it falls into categories 1 through 9,
c IGNORE IT and add it to the cmnt array

      if ((t.ge.1.and.t.le.9).or.
     &    (t.eq.28).or.
     &    (t.eq.29).or.
     &    (t.eq.30).or.
     &    (t.eq.31).or.
     &    (t.eq.32).or.
     &    (t.eq.33).or.
     &    (t.eq.34).or.
     &    (t.eq.35).or.
     &    (t.eq.36).or.
     &    (t.eq.37).or.
     &    (t.eq.38).or.
     &    (t.eq.39).or.
     &    (t.eq.999)) then
       cmnt(i1) = pdbline1
       i1 = i1 + 1
       goto 100
      end if

c If t = 10, then it's a remark... 
c There are currenly two pieces of data that we want
c to retrieve from remarks:
c
c	words such as DIMER, TETRAMER, etc.
c	BIOMT data

      if (t.eq.10) then
       if (pdbline1(8:10).eq.'350') then
        if (pdbline1(12:16).eq.'APPLY') then
         call startBIOMT(pdbline1,BIOMTstart,tempNchain,
     &      tempchain,art)
         cmnt(i1) = pdbline1
         i1 = i1 + 1
        else
         BIOMTstart = .false.  
        end if
        if (BIOMTstart) then
c Now, immediately read in the next 3 lines and send this
c data to a special parser surbroutine
         nTrans = nTrans + 1
         nChain(nTrans) = tempNchain
         do i2 = 1,nChain(nTrans)
          BMTchain(nTrans,i2) = tempchain(i2)
         end do
         read (1,10) biomtline(1)
         cmnt(i1) = biomtline(1)
         i1 = i1 + 1
         read (1,10) biomtline(2)
         cmnt(i1) = biomtline(2)
         i1 = i1 + 1
         read (1,10) biomtline(3)
         cmnt(i1) = biomtline(3)
         i1 = i1 + 1
         call parseBIOMT(nTrans,biomtline,TFMMat,TLTMat,art)
         BIOMT = .true.
        else
         if (BIOMT.and.(StringLength(pdbline1).gt.15)) then
          biomtline(1) = pdbline1
          cmnt(i1) = biomtline(1)  
          i1 = i1 + 1
          read (1,10) biomtline(2)
          cmnt(i1) = biomtline(2)
          i1 = i1 + 1
          read (1,10) biomtline(3)
          cmnt(i1) = biomtline(3)
          i1 = i1 + 1
          linelen = StringLength(biomtline(1))
          if (linelen.gt.15.and.
     &        biomtline(1)(14:18).eq.'BIOMT') then
           nTrans = nTrans + 1
           nChain(nTrans) = nChain(nTrans-1)
           do i2 = 1,nChain(nTrans)
            BMTchain(nTrans,i2) = BMTchain(nTrans-1,i2)
           end do
           call parseBIOMT(nTrans,biomtline,TFMMat,TLTMat,art)
          end if
          BIOMT = .true.
         end if
        end if
       end if
       if (pdbline1(8:10).eq.'300') then
        call oligomercheck(pdbline1,eOlig)
        cmnt(i1) = pdbline1
        i1 = i1 + 1
       end if
       goto 100
      end if

c Again, if its one of the following types, we want to simply
c store it without manipulation
      if (t.gt.10.and.t.le.21) then
       cmnt(i1) = pdbline1
       i1 = i1 + 1
       goto 100
      end if
 
      if (t.eq.22.or.t.eq.24) then

c It's an ATOM or HETATM record - parse the line and 
c extract the data

       call parseATOM(pdbline1,tempah,tempanum,tempanam,
     &  tempmc,tempres,tempc,temptru,tempsub,tempx,tempy,
     &  tempz,temprad,tempB)

c OK, now run some tests to decide if we're going to pay
c attention to this particular atom entry

c Ignore Water
c       if (tempres.eq.'H0H'.or.
c     &     tempres.eq.'D0D'.or.
c     &     tempres.eq.'S04'.or.
c     &     tempres.eq.'HOH'.or.
c     &     tempres.eq.'DOD'.or.
c     &     tempres.eq.'SO4') then
c        goto 100
c       end if  

c Ignore Multiple entries of a given atom
c       if (atcnt.gt.1.and.
c     &    tempres.eq.resnum(atcnt).and.
c     &    tempanam.eq.atomname(atcnt).and.
c     &    tempsub.ne.subnum(atcnt)) then
c        goto 100
c       end if 
 
c Ignore Hydrogen
       if (tempanam(2:2).eq.'H') then
        goto 100
       end if

       if (atcnt.eq.0) then
        lrnum = temptru
        frnum = temptru
        lsnum = tempsub
        nrnum = lrnum + 1
       end if
       atcnt = atcnt + 1
       atheader(atcnt) = tempah
       atomnum(atcnt)  = tempanum
       atomname(atcnt) = tempanam
       multiconf(atcnt) = tempmc
       res(atcnt) = tempres
       if (t.eq.24) then
         hetatm(atcnt) = .true.
       end if
       if (tempc.eq.' '.and.art) then
        chain(atcnt) = 'A'
       else
        chain(atcnt) = tempc
       end if
       trueResnum(atcnt) = temptru
       subnum(atcnt) = tempsub
       coord(1,atcnt) = tempx
       coord(2,atcnt) = tempy
       coord(3,atcnt) = tempz
       rad(atcnt) = temprad
       Bfac(atcnt) = tempB

c        print *,atheader(atcnt),atomnum(atcnt),atomname(atcnt),
c     &   multiconf(atcnt),res(atcnt),chain(atcnt),
c     &   trueResnum(atcnt),subnum(atcnt)

      end if

      if (t.eq.23) then
c it's a terminator line - do nothing!
       goto 100
      end if

      if (t.eq.25) then
c ignore it (it's a master line)
       goto 100
      end if

      if (t.eq.26) then
c ignore it (it's connectivity data)
       goto 100
      end if

      if (t.eq.27) then
c we're at the end of the pdbfile 
       goto 100
      end if

      if (atcnt.gt.0) then
       if (atcnt.gt.1) then
        if (chain(atcnt).ne.chain(atcnt-1)) then
         numchains = numchains + 1
        end if
       end if 

      if (trueResnum(atcnt).eq.lrnum.and.
     &    subnum(atcnt).eq.lsnum) then
       goto 300
      end if

      if (trueResnum(atcnt).eq.lrnum.and.
     &    subnum(atcnt).ne.lsnum) then
       numres = numres + 1
       resinc = resinc + 1
       goto 300
      end if

      if (trueResnum(atcnt).eq.nrnum.and.
     &    subnum(atcnt).eq.' ') then
       numres = numres + 1
       nrnum = nrnum + 1
       goto 300
      endif

      if (trueResnum(atcnt).gt.nrnum.and.
     &    subnum(atcnt).eq.' ') then
       resinc = resinc +
     &  (trueResnum(atcnt)-nrnum-1)
       numres = numres + 1
       nrnum = nrnum + resinc
       goto 300 
      endif

      if (trueResnum(atcnt).lt.lrnum) then
       resinc = 0
       numres = numres + 1
       frnum = trueResnum(atcnt)
       nrnum = trueResnum(atcnt) +1
       lrnum = trueResnum(atcnt)
       goto 300
      end if
      end if

300   if (atcnt.gt.0) then
      lrnum = trueResnum(atcnt)
      lsnum = subnum(atcnt)
      lanum = atomnum(atcnt)
      resnum(atcnt) = trueResnum(atcnt) +
     &  resinc - srnum
      subnum(atcnt) = ' '
      end if

      goto 100

900   close (1)
      numatm = atcnt 
      numcomlines = i1-1

      if (numatm*nTrans.gt.maxat) then
       print *,'transformations will produce too many atoms!
     &   increase maxat!'
       reject = .true.
       goto 999
      end if   

      origNumchains = numchains

      if (BIOMT) then
       print *,fname,' is a multimer'
       do i=1,nTrans
        do j=1,nChain(i)
         nat = 0
         lanum = atomnum(numatm)
         lrnum = resnum(numatm)
         do k=1,numatm
          if (chain(k).eq.BMTchain(i,j)) then
           nat = nat + 1
           newatplace = numatm + nat

           do l=1,3
            oldCoord(l) = coord(l,k)
            tempMatB(l) = TLTMat(i,l)
            do m=1,3
             tempMatA(l,m) = TFMMat(i,l,m)
            end do
           end do

           call transform(tempMatA,tempMatB,oldCoord,newCoord)

           coord(1,newatplace) = newCoord(1)
           coord(2,newatplace) = newCoord(2)
           coord(3,newatplace) = newCoord(3)

           atheader(newatplace) = atheader(k)
           atomname(newatplace) = atomname(k)
           multiconf(newatplace) = multiconf(k)
           res(newatplace) = res(k)
           subnum(newatplace) = subnum(k)

           chain(newatplace) = incChar(chain(k),numchains)

           atomnum(newatplace) = lanum + nat
           trueResnum(newatplace) = trueResnum(k)

           Bfac(newatplace) = Bfac(k)
           rad(newatplace) = rad(k)
           hetatm(newatplace) = hetatm(k)
                                         
          end if
         end do
         numchains = numchains + 1
         numatm = numatm + nat
        end do
       end do
      end if
      do i=1,numatm
       if (art.and.chain(i).eq.'A') then
        chain(i) = ' '
       end if
      end do
      if (art) then
       print *,'Artificial designation and removal of chain specifier'
      end if
      
      call findNextWord(fname,1,fnstart,fnend,exist)
      if (exist) then
       tempfname= 'mod.'//fname(fnstart:fnend)
      end if
      goto 997

995   print *,'ERROR opening pdbfile!'
      stop


997   if (BIOMT) then
       call pdbwrite (tempfname)
      end if

999    end

c------------------------------------------------------
      subroutine pdbread_clean (fname,dir)

      parameter (mxBMT = 50)

      include 'parameters.h'
      include 'pdb.h'

      character*80 pdbline1
      character*80 biomtline(3)

      character*(*) fname 
      integer fnstart
      integer fnend
      character*50 tempfname
      character*120 concat
      character*(*) dir

      integer i1,i2
      integer j1,j2,j3
      integer i,j,k,l,m
      integer t
      integer eOlig
      integer resinc
      integer srnum
      real    TFMMat(mxBMT,3,3)
      real    TLTMat(mxBMT,3)
      character BMTchain(mxBMT,20)
      character tempchain(mxBMT)
      logical BIOMTstart
      logical BIOMT
      logical art
      integer nTrans
      integer nChain(mxBMT)
      integer tempNchain
 
      integer atcnt

      character*6 tempah
      character*4 tempanam
      character*3 tempres
      character*1 tempsub
      character*1 tempc
      character*1 tempmc
      integer     tempanum
      integer     temptru
      real        tempx
      real        tempy
      real        tempz
      real        temprad
      real        tempB  

      integer frnum,nrnum,lrnum
      integer same(20)
      integer count
      logical diff
      character*1 lsnum
      integer lanum

      integer nat
      integer origNumChains
      integer newatplace
      real    oldCoord(3)
      real    newCoord(3)
      character IncChar
      real tempMatA(3,3)
      real tempMatB(3)
      logical exist

      integer StringLength
      integer linelen
 

      resinc = 0
      srnum = 0
      do j1 = 1,mxBMT
       do j2 = 1,3
        do j3 = 1,3
         TFMMat(j1,j2,j3) = 0.0
        end do
        TLTMat(j1,j2) = 0.0
       end do
       nChain(j1) = 0
      end do
      nTrans = 0
      atcnt = 0
      art = .false.
      numchains = 1
      BIOMT = .false.

10    format (a80)

      call fileexist(fname,dir,concat,exist)
      if (.not.exist) then
       goto 995
      end if

      i1 = 1
      open (unit=1,err = 995,file=concat,status='OLD')

100   read (1,10,end=900) pdbline1
      call removeTrailingSpaces(pdbline1)

      call catEntry(pdbline1,t)


c OK, now if it falls into categories 1 through 9,
c IGNORE IT and add it to the cmnt array

      if ((t.ge.1.and.t.le.9).or.
     &    (t.eq.28).or.
     &    (t.eq.29).or.
     &    (t.eq.30).or.
     &    (t.eq.31).or.
     &    (t.eq.32).or.
     &    (t.eq.33).or.
     &    (t.eq.34).or.
     &    (t.eq.35).or.
     &    (t.eq.36).or.
     &    (t.eq.37).or.
     &    (t.eq.38).or.
     &    (t.eq.39).or.
     &    (t.eq.999)) then
       cmnt(i1) = pdbline1
       i1 = i1 + 1
       goto 100
      end if

c If t = 10, then it's a remark... 
c There are currenly two pieces of data that we want
c to retrieve from remarks:
c
c	words such as DIMER, TETRAMER, etc.
c	BIOMT data

      if (t.eq.10) then
       if (pdbline1(8:10).eq.'350') then
        if (pdbline1(12:16).eq.'APPLY') then
         
c         print *,'startBIOMT found at the following line'
c         print *, pdbline1      

         call startBIOMT(pdbline1,BIOMTstart,tempNchain,
     &      tempchain,art)
         cmnt(i1) = pdbline1
         i1 = i1 + 1
        else
         BIOMTstart = .false.  
        end if
        if (BIOMTstart) then
c Now, immediately read in the next 3 lines and send this
c data to a special parser surbroutine
         nTrans = nTrans + 1
         nChain(nTrans) = tempNchain
         do i2 = 1,nChain(nTrans)
          BMTchain(nTrans,i2) = tempchain(i2)
         end do
         read (1,10) biomtline(1)
         call removeTrailingSpaces(biomtline(1))
         cmnt(i1) = biomtline(1)
         i1 = i1 + 1
         read (1,10) biomtline(2)
         call removeTrailingSpaces(biomtline(2))
         cmnt(i1) = biomtline(2)
         i1 = i1 + 1
         read (1,10) biomtline(3)
         call removeTrailingSpaces(biomtline(3))
         cmnt(i1) = biomtline(3)
         i1 = i1 + 1
         call parseBIOMT(nTrans,biomtline,TFMMat,TLTMat,art)
         call removeTrailingSpaces(pdbline1)
         BIOMT = .true.
        else
         if (BIOMT.and.(StringLength(pdbline1).gt.26)) then

c          print *,'magic point reached'
c          print *,pdbline1

          biomtline(1) = pdbline1
          cmnt(i1) = biomtline(1)  
          i1 = i1 + 1
          read (1,10) biomtline(2)
          call removeTrailingSpaces(biomtline(2))
          cmnt(i1) = biomtline(2)
          i1 = i1 + 1
          read (1,10) biomtline(3)
          call removeTrailingSpaces(biomtline(3))
          cmnt(i1) = biomtline(3)
          i1 = i1 + 1
          linelen = StringLength(biomtline(1))
          if (linelen.gt.15.and.
     &        biomtline(1)(14:18).eq.'BIOMT') then
           nTrans = nTrans + 1
           nChain(nTrans) = nChain(nTrans-1)
           do i2 = 1,nChain(nTrans)
            BMTchain(nTrans,i2) = BMTchain(nTrans-1,i2)
           end do
           call parseBIOMT(nTrans,biomtline,TFMMat,TLTMat,art)
          end if
          BIOMT = .true.
         end if
        end if
       end if
       if (pdbline1(8:10).eq.'300') then
        call oligomercheck(pdbline1,eOlig)
        cmnt(i1) = pdbline1
        i1 = i1 + 1
       end if
       goto 100
      end if

c Again, if its one of the following types, we want to simply
c store it without manipulation
      if (t.gt.10.and.t.le.21) then
       cmnt(i1) = pdbline1
       i1 = i1 + 1
       goto 100
      end if
 
      if (t.eq.22.or.t.eq.24) then

c It's an ATOM or HETATM record - parse the line and 
c extract the data

       call parseATOM(pdbline1,tempah,tempanum,tempanam,
     &  tempmc,tempres,tempc,temptru,tempsub,tempx,tempy,
     &  tempz,temprad,tempB)

c OK, now run some tests to decide if we're going to pay
c attention to this particular atom entry

c Since we're "cleaning" the pdb file here, we're going to
c discard ANY heteroatom record we come across.
       if (tempah.eq.'HETATM') then
         goto 100
       end if

c Ignore Water
       if (tempres.eq.'H0H'.or.
     &     tempres.eq.'D0D'.or.
     &     tempres.eq.'S04'.or.
     &     tempres.eq.'HOH'.or.
     &     tempres.eq.'DOD'.or.
     &     tempres.eq.'SO4') then
        goto 100
       end if  

c Ignore Multiple entries of a given atom
       if (atcnt.gt.1.and.
     &     (tempres.eq.res(atcnt)).and.
     &     (tempanam.eq.atomname(atcnt)).and.
     &     (tempsub.ne.subnum(atcnt))) then
        print *,'Found multiple entry for residue ',tempres,temptru
        goto 100
       end if 
 
c Ignore Hydrogen
       if (tempanam(2:2).eq.'H') then
        goto 100
       end if

       if (atcnt.eq.0) then
        lrnum = temptru
        frnum = temptru
        lsnum = tempsub
        nrnum = lrnum + 1
       end if
       atcnt = atcnt + 1
       atheader(atcnt) = tempah
       atomnum(atcnt)  = tempanum
       atomname(atcnt) = tempanam
       multiconf(atcnt) = tempmc
       res(atcnt) = tempres
       if (t.eq.24) then
         hetatm(atcnt) = .true.
       end if
       if (tempc.eq.' '.and.art) then
        chain(atcnt) = 'A'
       else
        chain(atcnt) = tempc
       end if
       trueResnum(atcnt) = temptru
       subnum(atcnt) = tempsub
       coord(1,atcnt) = tempx
       coord(2,atcnt) = tempy
       coord(3,atcnt) = tempz
       rad(atcnt) = temprad
       Bfac(atcnt) = tempB

c        print *,atheader(atcnt),atomnum(atcnt),atomname(atcnt),
c     &   multiconf(atcnt),res(atcnt),chain(atcnt),
c     &   trueResnum(atcnt),subnum(atcnt)

      end if

      if (t.eq.23) then
c it's a terminator line - do nothing!
       goto 100
      end if

      if (t.eq.25) then
c ignore it (it's a master line)
       goto 100
      end if

      if (t.eq.26) then
c ignore it (it's connectivity data)
       goto 100
      end if

      if (t.eq.27) then
c we're at the end of the pdbfile 
       goto 100
      end if

      if (atcnt.gt.0) then
       if (atcnt.gt.1) then
        if (chain(atcnt).ne.chain(atcnt-1)) then
         numchains = numchains + 1
        end if
       end if 

      if (trueResnum(atcnt).eq.lrnum.and.
     &    subnum(atcnt).eq.lsnum) then
       goto 300
      end if

      if (trueResnum(atcnt).eq.lrnum.and.
     &    subnum(atcnt).ne.lsnum) then
       numres = numres + 1
       resinc = resinc + 1
       goto 300
      end if

      if (trueResnum(atcnt).eq.nrnum.and.
     &    subnum(atcnt).eq.' ') then
       numres = numres + 1
       nrnum = nrnum + 1
       goto 300
      endif

      if (trueResnum(atcnt).gt.nrnum.and.
     &    subnum(atcnt).eq.' ') then
       resinc = resinc +
     &  (trueResnum(atcnt)-nrnum-1)
       numres = numres + 1
       nrnum = nrnum + resinc
       goto 300 
      endif

      if (trueResnum(atcnt).lt.lrnum) then
       resinc = 0
       numres = numres + 1
       frnum = trueResnum(atcnt)
       nrnum = trueResnum(atcnt) +1
       lrnum = trueResnum(atcnt)
       goto 300
      end if
      end if

300   if (atcnt.gt.0) then
      lrnum = trueResnum(atcnt)
      lsnum = subnum(atcnt)
      lanum = atomnum(atcnt)
      resnum(atcnt) = trueResnum(atcnt) +
     &  resinc - srnum
      subnum(atcnt) = ' '
      end if

      goto 100

900   close (1)
      numatm = atcnt 
      numcomlines = i1-1

      if (numatm*nTrans.gt.maxat) then
       print *,'transformations will produce too many atoms!
     &   increase maxat!'
       reject = .true.
       goto 999
      end if   

      origNumchains = numchains

      if (BIOMT) then

c First, check to see how many 'real' transformations exist
c (some pdb files list the transformations twice (for some
c silly reason)

      if (art) then
       count = 0
       if (nTrans.gt.1) then 
         do itemp=1,nTrans-1
          do jtemp= nTrans, itemp+1 ,-1
           if (BMTchain(jtemp,1).eq.BMTchain(itemp,1)) then
             diff = .false.
              do itemp2 = 1,3
               if (TLTMat(itemp,itemp2).ne.TLTMat(jtemp,itemp2)) then
                 diff = .true.
               end if
               do jtemp2 = 1,3
                 if (TFMMat(itemp,itemp2,jtemp2).ne.
     &               TFMMat(jtemp,itemp2,jtemp2)) then
                  diff = .true.
                 end if
               end do
              end do
             if (.not.diff) then
               count = count + 1
               same(count) = jtemp;
             end if 
           end if 
          end do
         end do
        end if
       end if

c       print *,( same(count),count=1,20)
            
       print *,fname,' is a multimer'
       do i=1,nTrans
        diff = .true.
        do i2 = 1,count
         if (i.eq.same(i2)) then
           diff = .false.
         end if
        end do
        if (diff) then
        do j=1,nChain(i)
         nat = 0
         lanum = atomnum(numatm)
         lrnum = resnum(numatm)
         do k=1,numatm
          if (chain(k).eq.BMTchain(i,j)) then
           nat = nat + 1
           newatplace = numatm + nat

           do l=1,3
            oldCoord(l) = coord(l,k)
            tempMatB(l) = TLTMat(i,l)
            do m=1,3
             tempMatA(l,m) = TFMMat(i,l,m)
            end do
           end do

           call transform(tempMatA,tempMatB,oldCoord,newCoord)

           coord(1,newatplace) = newCoord(1)
           coord(2,newatplace) = newCoord(2)
           coord(3,newatplace) = newCoord(3)

           atheader(newatplace) = atheader(k)
           atomname(newatplace) = atomname(k)
           multiconf(newatplace) = multiconf(k)
           res(newatplace) = res(k)
           subnum(newatplace) = subnum(k)

           chain(newatplace) = incChar(chain(k),numchains+1)

           atomnum(newatplace) = lanum + nat
           trueResnum(newatplace) = trueResnum(k)

           Bfac(newatplace) = Bfac(k)
           rad(newatplace) = rad(k)
           hetatm(newatplace) = hetatm(k)
                                         
          end if
         end do
         numchains = numchains + 1
         numatm = numatm + nat
        end do
        end if
       end do
      end if
      do i=1,numatm
       if (art.and.chain(i).eq.'A') then
        chain(i) = ' '
       end if
      end do
      if (art) then
       print *,'Artificial designation and removal of chain specifier'
      end if
      
      call findNextWord(fname,1,fnstart,fnend,exist)
      if (exist) then
       tempfname= 'mod.'//fname(fnstart:fnend)
      end if
      goto 997

995   print *,'ERROR opening pdbfile!'
      stop


997   if (BIOMT) then
c       call pdbwrite (tempfname)
      end if

999    end

c------------------------------------------------------
      subroutine pdbread_clean_nobiomt_with_h (fname,dir)

      parameter (mxBMT = 50)

      include 'parameters.h'
      include 'pdb.h'

      character*80 pdbline1
      character*80 biomtline(3)

      character*(*) fname 
      integer fnstart
      integer fnend
      character*50 tempfname
      character*120 concat
      character*(*) dir

      integer i1,i2
      integer j1,j2,j3
      integer i,j,k,l,m
      integer t
      integer eOlig
      integer resinc
      integer srnum
      real    TFMMat(mxBMT,3,3)
      real    TLTMat(mxBMT,3)
      character BMTchain(mxBMT,20)
      character tempchain(mxBMT)
      logical BIOMTstart
      logical BIOMT
      logical art
      integer nTrans
      integer nChain(mxBMT)
      integer tempNchain
 
      integer atcnt

      character*6 tempah
      character*4 tempanam
      character*3 tempres
      character*1 tempsub
      character*1 tempc
      character*1 tempmc
      integer     tempanum
      integer     temptru
      real        tempx
      real        tempy
      real        tempz
      real        temprad
      real        tempB  

      integer frnum,nrnum,lrnum
      integer same(20)
      integer count
      logical diff
      character*1 lsnum
      integer lanum

      integer nat
      integer origNumChains
      integer newatplace
      real    oldCoord(3)
      real    newCoord(3)
      character IncChar
      real tempMatA(3,3)
      real tempMatB(3)
      logical exist

      integer StringLength
      integer linelen
 

      resinc = 0
      srnum = 0
      do j1 = 1,mxBMT
       do j2 = 1,3
        do j3 = 1,3
         TFMMat(j1,j2,j3) = 0.0
        end do
        TLTMat(j1,j2) = 0.0
       end do
       nChain(j1) = 0
      end do
      nTrans = 0
      atcnt = 0
      art = .false.
      numchains = 1
      BIOMT = .false.

10    format (a80)

      call fileexist(fname,dir,concat,exist)
      if (.not.exist) then
       goto 995
      end if

      i1 = 1
      open (unit=1,err = 995,file=concat,status='OLD')

100   read (1,10,end=900) pdbline1
      call removeTrailingSpaces(pdbline1)

      call catEntry(pdbline1,t)


c OK, now if it falls into categories 1 through 9,
c IGNORE IT and add it to the cmnt array

      if ((t.ge.1.and.t.le.9).or.
     &    (t.eq.28).or.
     &    (t.eq.29).or.
     &    (t.eq.30).or.
     &    (t.eq.31).or.
     &    (t.eq.32).or.
     &    (t.eq.33).or.
     &    (t.eq.34).or.
     &    (t.eq.35).or.
     &    (t.eq.36).or.
     &    (t.eq.37).or.
     &    (t.eq.38).or.
     &    (t.eq.39).or.
     &    (t.eq.999)) then
       cmnt(i1) = pdbline1
       i1 = i1 + 1
       goto 100
      end if

c If t = 10, then it's a remark... 
c There are currenly two pieces of data that we want
c to retrieve from remarks:
c
c	words such as DIMER, TETRAMER, etc.
c	BIOMT data

      if (t.eq.10) then
       if (pdbline1(8:10).eq.'350') then
        if (pdbline1(12:16).eq.'APPLY') then
         
c         print *,'startBIOMT found at the following line'
c         print *, pdbline1      

         call startBIOMT(pdbline1,BIOMTstart,tempNchain,
     &      tempchain,art)
         cmnt(i1) = pdbline1
         i1 = i1 + 1
        else
         BIOMTstart = .false.  
        end if
        if (BIOMTstart) then
c Now, immediately read in the next 3 lines and send this
c data to a special parser surbroutine
         nTrans = nTrans + 1
         nChain(nTrans) = tempNchain
         do i2 = 1,nChain(nTrans)
          BMTchain(nTrans,i2) = tempchain(i2)
         end do
         read (1,10) biomtline(1)
         call removeTrailingSpaces(biomtline(1))
         cmnt(i1) = biomtline(1)
         i1 = i1 + 1
         read (1,10) biomtline(2)
         call removeTrailingSpaces(biomtline(2))
         cmnt(i1) = biomtline(2)
         i1 = i1 + 1
         read (1,10) biomtline(3)
         call removeTrailingSpaces(biomtline(3))
         cmnt(i1) = biomtline(3)
         i1 = i1 + 1
         call parseBIOMT(nTrans,biomtline,TFMMat,TLTMat,art)
         call removeTrailingSpaces(pdbline1)
         BIOMT = .true.
        else
         if (BIOMT.and.(StringLength(pdbline1).gt.26)) then

c          print *,'magic point reached'
c          print *,pdbline1

          biomtline(1) = pdbline1
          cmnt(i1) = biomtline(1)  
          i1 = i1 + 1
          read (1,10) biomtline(2)
          call removeTrailingSpaces(biomtline(2))
          cmnt(i1) = biomtline(2)
          i1 = i1 + 1
          read (1,10) biomtline(3)
          call removeTrailingSpaces(biomtline(3))
          cmnt(i1) = biomtline(3)
          i1 = i1 + 1
          linelen = StringLength(biomtline(1))
          if (linelen.gt.15.and.
     &        biomtline(1)(14:18).eq.'BIOMT') then
           nTrans = nTrans + 1
           nChain(nTrans) = nChain(nTrans-1)
           do i2 = 1,nChain(nTrans)
            BMTchain(nTrans,i2) = BMTchain(nTrans-1,i2)
           end do
           call parseBIOMT(nTrans,biomtline,TFMMat,TLTMat,art)
          end if
          BIOMT = .true.
         end if
        end if
       end if
       if (pdbline1(8:10).eq.'300') then
        call oligomercheck(pdbline1,eOlig)
        cmnt(i1) = pdbline1
        i1 = i1 + 1
       end if
       goto 100
      end if

c Again, if its one of the following types, we want to simply
c store it without manipulation
      if (t.gt.10.and.t.le.21) then
       cmnt(i1) = pdbline1
       i1 = i1 + 1
       goto 100
      end if
 
      if (t.eq.22.or.t.eq.24) then

c It's an ATOM or HETATM record - parse the line and 
c extract the data

       call parseATOM(pdbline1,tempah,tempanum,tempanam,
     &  tempmc,tempres,tempc,temptru,tempsub,tempx,tempy,
     &  tempz,temprad,tempB)

c OK, now run some tests to decide if we're going to pay
c attention to this particular atom entry

c Since we're "cleaning" the pdb file here, we're going to
c discard ANY heteroatom record we come across.
       if (tempah.eq.'HETATM') then
         goto 100
       end if

c Ignore Water
       if (tempres.eq.'H0H'.or.
     &     tempres.eq.'D0D'.or.
     &     tempres.eq.'S04'.or.
     &     tempres.eq.'HOH'.or.
     &     tempres.eq.'DOD'.or.
     &     tempres.eq.'SO4') then
        goto 100
       end if  

c Ignore Multiple entries of a given atom
       if (atcnt.gt.1.and.
     &     (tempres.eq.res(atcnt)).and.
     &     (tempanam.eq.atomname(atcnt)).and.
     &     (tempsub.ne.subnum(atcnt))) then
        print *,'Found multiple entry for residue ',tempres,temptru
        goto 100
       end if 
 
c Ignore Hydrogen
c       if (tempanam(2:2).eq.'H') then
c        goto 100
c       end if

       if (atcnt.eq.0) then
        lrnum = temptru
        frnum = temptru
        lsnum = tempsub
        nrnum = lrnum + 1
       end if
       atcnt = atcnt + 1
       atheader(atcnt) = tempah
       atomnum(atcnt)  = tempanum
       atomname(atcnt) = tempanam
       multiconf(atcnt) = tempmc
       res(atcnt) = tempres
       if (t.eq.24) then
         hetatm(atcnt) = .true.
       end if
       if (tempc.eq.' '.and.art) then
        chain(atcnt) = 'A'
       else
        chain(atcnt) = tempc
       end if
       trueResnum(atcnt) = temptru
       subnum(atcnt) = tempsub
       coord(1,atcnt) = tempx
       coord(2,atcnt) = tempy
       coord(3,atcnt) = tempz
       rad(atcnt) = temprad
       Bfac(atcnt) = tempB

c        print *,atheader(atcnt),atomnum(atcnt),atomname(atcnt),
c     &   multiconf(atcnt),res(atcnt),chain(atcnt),
c     &   trueResnum(atcnt),subnum(atcnt)

      end if

      if (t.eq.23) then
c it's a terminator line - do nothing!
       goto 100
      end if

      if (t.eq.25) then
c ignore it (it's a master line)
       goto 100
      end if

      if (t.eq.26) then
c ignore it (it's connectivity data)
       goto 100
      end if

      if (t.eq.27) then
c we're at the end of the pdbfile 
       goto 100
      end if

      if (atcnt.gt.0) then
       if (atcnt.gt.1) then
        if (chain(atcnt).ne.chain(atcnt-1)) then
         numchains = numchains + 1
        end if
       end if 

      if (trueResnum(atcnt).eq.lrnum.and.
     &    subnum(atcnt).eq.lsnum) then
       goto 300
      end if

      if (trueResnum(atcnt).eq.lrnum.and.
     &    subnum(atcnt).ne.lsnum) then
       numres = numres + 1
       resinc = resinc + 1
       goto 300
      end if

      if (trueResnum(atcnt).eq.nrnum.and.
     &    subnum(atcnt).eq.' ') then
       numres = numres + 1
       nrnum = nrnum + 1
       goto 300
      endif

      if (trueResnum(atcnt).gt.nrnum.and.
     &    subnum(atcnt).eq.' ') then
       resinc = resinc +
     &  (trueResnum(atcnt)-nrnum-1)
       numres = numres + 1
       nrnum = nrnum + resinc
       goto 300 
      endif

      if (trueResnum(atcnt).lt.lrnum) then
       resinc = 0
       numres = numres + 1
       frnum = trueResnum(atcnt)
       nrnum = trueResnum(atcnt) +1
       lrnum = trueResnum(atcnt)
       goto 300
      end if
      end if

300   if (atcnt.gt.0) then
      lrnum = trueResnum(atcnt)
      lsnum = subnum(atcnt)
      lanum = atomnum(atcnt)
      resnum(atcnt) = trueResnum(atcnt) +
     &  resinc - srnum
      subnum(atcnt) = ' '
      end if

      goto 100

900   close (1)
      numatm = atcnt 
      numcomlines = i1-1

      if (numatm*nTrans.gt.maxat) then
       print *,'transformations will produce too many atoms!
     &   increase maxat!'
       reject = .true.
       goto 999
      end if   

      origNumchains = numchains

c  we're going the skip the generation of the biomolecule
      goto 999

      if (BIOMT) then

c First, check to see how many 'real' transformations exist
c (some pdb files list the transformations twice (for some
c silly reason)

      if (art) then
       count = 0
       if (nTrans.gt.1) then 
         do itemp=1,nTrans-1
          do jtemp= nTrans, itemp+1 ,-1
           if (BMTchain(jtemp,1).eq.BMTchain(itemp,1)) then
             diff = .false.
              do itemp2 = 1,3
               if (TLTMat(itemp,itemp2).ne.TLTMat(jtemp,itemp2)) then
                 diff = .true.
               end if
               do jtemp2 = 1,3
                 if (TFMMat(itemp,itemp2,jtemp2).ne.
     &               TFMMat(jtemp,itemp2,jtemp2)) then
                  diff = .true.
                 end if
               end do
              end do
             if (.not.diff) then
               count = count + 1
               same(count) = jtemp;
             end if 
           end if 
          end do
         end do
        end if
       end if

c       print *,( same(count),count=1,20)
            
       print *,fname,' is a multimer'
       do i=1,nTrans
        diff = .true.
        do i2 = 1,count
         if (i.eq.same(i2)) then
           diff = .false.
         end if
        end do
        if (diff) then
        do j=1,nChain(i)
         nat = 0
         lanum = atomnum(numatm)
         lrnum = resnum(numatm)
         do k=1,numatm
          if (chain(k).eq.BMTchain(i,j)) then
           nat = nat + 1
           newatplace = numatm + nat

           do l=1,3
            oldCoord(l) = coord(l,k)
            tempMatB(l) = TLTMat(i,l)
            do m=1,3
             tempMatA(l,m) = TFMMat(i,l,m)
            end do
           end do

           call transform(tempMatA,tempMatB,oldCoord,newCoord)

           coord(1,newatplace) = newCoord(1)
           coord(2,newatplace) = newCoord(2)
           coord(3,newatplace) = newCoord(3)

           atheader(newatplace) = atheader(k)
           atomname(newatplace) = atomname(k)
           multiconf(newatplace) = multiconf(k)
           res(newatplace) = res(k)
           subnum(newatplace) = subnum(k)

           chain(newatplace) = incChar(chain(k),numchains+1)

           atomnum(newatplace) = lanum + nat
           trueResnum(newatplace) = trueResnum(k)

           Bfac(newatplace) = Bfac(k)
           rad(newatplace) = rad(k)
           hetatm(newatplace) = hetatm(k)
                                         
          end if
         end do
         numchains = numchains + 1
         numatm = numatm + nat
        end do
        end if
       end do
      end if
      do i=1,numatm
       if (art.and.chain(i).eq.'A') then
        chain(i) = ' '
       end if
      end do
      if (art) then
       print *,'Artificial designation and removal of chain specifier'
      end if
      
      call findNextWord(fname,1,fnstart,fnend,exist)
      if (exist) then
       tempfname= 'mod.'//fname(fnstart:fnend)
      end if
      goto 997

995   print *,'ERROR opening pdbfile!'
      stop


997   if (BIOMT) then
c       call pdbwrite (tempfname)
      end if

999    end
c----------------------------------------------------------
      subroutine catEntry(line,type)

c these are the different values for pdb line 'types'
c--------------
c 1 - HEADER		When a line is read in, we want to
c 2 - TITLE		be able to assign a 'type' to that
c 3 - COMPND		that line, which will be sent back
c 4 - SOURCE		to the calling subroutine.
c 5 - KEYWDS	
c 6 - EXPDTA		The type will be an integer value.
c 7 - AUTHOR
c 8 - REVDAT
c 9 - JRNL
c 10- REMARK
c 11- DBREF
c 12- SEQRES
c 13- FORMUL
c 14- HELIX
c 15- SHEET
c 16- CISPEP
c 17- SITE
c 18- CRYST
c 19- ORIGX
c 20- SCALE
c 21- MTRIX
c 22- ATOM
c 23- TER
c 24- HETATM
c 25- MASTER
c 26- CONECT
c 27- END
c 28- FTNOTE
c 29- SEQADV
c 30- TURN
c 31- HET
c 32- HETNAM 
c 33- SSBOND 
c 34- LINK
c 35- MODRES
c 36- SLTBRG
c 37- SPRSDE
c 38- MODEL
c 39- ANISOU 
c 40- USER - this is specific to output from Jane Richardson's reduce program
c--------------

      character*(*) line
      integer type
      integer start
      integer finish
      logical exist
 
c Find the start and end of the first word

      call findNextWord(line,1,start,finish,exist)

      if (start.ne.1.or.(.not.exist)) then
       print *,'error reading pdbfile'
       stop
      end if

        if (line(start:start+5).eq.'HEADER') then
         type = 1
         return
        end if
        if (line(start:start+5).eq.'COMPND') then
         type = 3
         return
        end if
        if (line(start:start+5).eq.'SOURCE') then
         type = 4
         return
        end if 
        if (line(start:start+5).eq.'KEYWDS') then
         type = 5
         return
        end if 
        if (line(start:start+5).eq.'EXPDTA') then
         type = 6
         return
        end if 
        if (line(start:start+5).eq.'AUTHOR') then
         type = 7
         return
        end if 
        if (line(start:start+5).eq.'REVDAT') then
         type = 8
         return
        end if 
        if (line(start:start+5).eq.'REMARK') then
         type = 10
         return
        end if 
        if (line(start:start+5).eq.'SEQRES') then
         type = 12
         return
        end if 
        if (line(start:start+5).eq.'FORMUL') then
         type = 13
         return
        end if
        if (line(start:start+5).eq.'CISPEP') then
         type = 16
         return
        end if 
        if (line(start:start+5).eq.'HETATM') then
         type = 24
         return
        end if 
        if (line(start:start+5).eq.'MASTER') then
         type = 25
         return
        end if 
        if (line(start:start+5).eq.'CONECT') then
         type = 26
         return
        end if
        if (line(start:start+4).eq.'CRYST') then
         type = 18
         return
        end if
        if (line(start:start+4).eq.'ORIGX') then
         type = 19
         return
        end if
    
        if (line(start:start+4).eq.'TITLE') then
         type = 2
         return
        end if 
        if (line(start:start+4).eq.'DBREF') then
         type = 11
         return
        end if 
        if (line(start:start+4).eq.'HELIX') then
         type = 14
         return
        end if 
        if (line(start:start+4).eq.'SHEET') then
         type = 15
         return
        end if 
        if (line(start:start+4).eq.'SCALE') then
         type = 20
         return
        end if 
        if (line(start:start+4).eq.'MTRIX') then
         type = 21
         return
        end if 

        if (line(start:start+3).eq.'JRNL'.or.
     &      line(start:start+3).eq.'JNRL') then
         type = 9
         return
        end if 
        if (line(start:start+3).eq.'SITE') then
         type = 17
         return
        end if 
        if (line(start:start+3).eq.'ATOM') then
         type = 22
         return
        end if 

        if (line(start:start+2).eq.'TER') then
         type = 23
         return
        end if 
        if (line(start:start+2).eq.'END') then
         type = 27
         return
        end if 
        if (line(start:start+5).eq.'FTNOTE') then
         type = 28
         return 
        end if
        if (line(start:start+5).eq.'SEQADV') then
         type = 29
         return
        end if  
        if (line(start:start+3).eq.'TURN') then
         type = 30
         return
        end if
        if (line(start:start+2).eq.'HET') then
         type = 31
         return
        end if 
        if (line(start:start+5).eq.'HETNAM') then
         type = 32
         return
        end if
        if (line(start:start+5).eq.'SSBOND') then
         type = 33
         return
        end if 
        if (line(start:start+3).eq.'LINK') then
         type = 34
         return
        end if 
        if (line(start:start+5).eq.'MODRES') then
         type = 35
         return
        end if 
        if (line(start:start+5).eq.'SLTBRG') then
         type = 36
         return
        end if 
        if (line(start:start+5).eq.'SPRSDE') then
         type = 37
         return
        end if 
        if (line(start:start+4).eq.'MODEL') then
         type = 38
         return
        end if
        if (line(start:start+5).eq.'ANISOU') then
         type = 39
         return
        end if  
        if (line(start:start+3).eq.'USER') then
         type = 40 
         return
        end if  

c      print *,'unable to determine the entry type of ',
c     &  line(start:finish)
       type = 999
c      stop
      end

c-----------------------------------------------------------------
      subroutine oligomercheck (line,type)

c This subroutine parses REMARK 300 entries to find the expected
c oligomeric state of the biomolecule

      character*(*) line
      character*10 string
      integer type
      integer dummy
      integer length
      integer StringLength

      logical found

      found = .false.
      length = StringLength(line)
 
c      print *,'in oligomercheck'
c      print *,'length = ',length

      if (length.lt.15) then
c       print *,'Too short!'
       return
      end if

      string = 'DIMER'
      call searchString(line,string,found,dummy)
       if (found) then
        print *,line
        type = 2
        return
       end if

      string = 'TRIMER'
      call searchString(line,string,found,dummy)
      if (found) then
       print *,line
       type = 3
       return
      end if

      string = 'TETRAMER'
      call searchString(line,string,found,dummy)
      if (found) then
       print *,line
       type = 4
       return
      end if

      string = 'PENTAMER'
      call searchString(line,string,found,dummy)
      if (found) then
       print *,line
       type = 5
       return
      end if

      string = 'HEXAMER'
      call searchString(line,string,found,dummy)
      if (found) then
       print *,line
       type = 6
       return
      end if

      string = 'HEPTAMER'
      call searchString(line,string,found,dummy)
      if (found) then
       print *,line
       type = 7
       return
      end if 

      string = 'OCTAMER'
      call searchString(line,string,found,dummy)
      if (found) then
       print *,line
       type = 8
       return
      end if

      string = 'MER'
      call searchString(line,string,found,dummy)
      if (found) then
       print *,'ERROR 5 - I can''t correctly parse:'
       print *,line
      end if

      end

c--------------------------------------------------------------
      subroutine startBIOMT(line,BIOMTstart,nChain,tempchain,art)

      parameter (mxBMT = 50)

      character*(*) line

      integer i1,i2
      integer linelen
      integer StringLength
      integer x
      integer start
      integer finish
      integer colonplace
      integer nChain

      character tempchain(mxBMT)

      logical colonfound
      logical chainfound
      logical isitAlpha
      logical BIOMTstart
      logical art
      logical exist

      linelen = StringLength(line)
      BIOMTstart = .false.
      if (linelen.le.26) then
       return
      end if
 
      nChain = 0
c starting at placeholder 11 in the string - search for the next word
      x = 11      
      call findNextWord(line,x,start,finish,exist)
      if (exist) then
       if (line(start:finish).eq.'APPLY') then
        call searchString(line,':',colonfound,colonplace)     
        if (colonplace.ne.999) then
         i2 = 0
         do i1 = colonplace+1,linelen
          chainfound = isitAlpha(line(i1:i1))
          if (chainfound) then
            i2 = i2 + 1
            tempchain(i2) = line(i1:i1)
            chainfound = .false.
            BIOMTstart = .true.
          end if
         end do
         nChain = i2
        end if
       end if
      end if

      if (nChain.eq.4) then
       if (tempchain(1).eq.'N') then
        if (tempchain(2).eq.'U') then
         if (tempchain(3).eq.'L') then
          if (tempchain(4).eq.'L') then
            nChain = 1
            art = .true.
            do i1 = 1,4
             tempchain(i1) = 'A'
            end do
          end if
         end if
        end if
       end if
      end if

      if (nChain.eq.0) then
c       print *, 'no formal chain designation in BIOMT line'
       BIOMTstart = .true.
       art = .true.
       tempchain(1) = 'A'
       nChain = 1
      end if

c      print *,nChain,(tempchain(i),i=1,nChain)
      end

c-----------------------------------------------------
      subroutine parseBIOMT (nTrans,biomtline,TFMMat,
     &   TLTMat,art)

      parameter (mxBMT = 50)
      character*80 biomtline(3)
      integer nTrans
      integer StringLength
      integer len
      integer x
      integer begin,finish
      real temp
      real    TFMMat(mxBMT,3,3)
      real    TLTMat(mxBMT,3)    
      logical exist
      logical art
      logical success

      success = .false.
      do i = 1,3
       len = StringLength(biomtline(i))
       if (len.lt.26) then
        print *,'Error parsing BIOMT line #',i,':'
        print *, biomtline(i)
        print *
        return
       end if
       x = 24 
       do j = 1,3
        call findNextWord(biomtline(i),x,begin,finish,exist)
        if (exist) then
         call strtoreal(biomtline(i)(begin:finish),temp)
         TFMMat(nTrans,i,j) = temp         
         x = finish + 1
        else
         print *,'ERRoR ReaDIng BIOmt rEcoRd!'
        end if
       end do
       call findNextWord(biomtline(i),x,begin,finish,exist)
       if (exist) then
        call strtoreal(biomtline(i)(begin:finish),temp)
        TLTMat(nTrans,i) = temp
       end if
      end do
      success = .true.
      end

c--------------------------------------------------------
      subroutine parseATOM(templine,tempah,tempanum,tempanam,
     &  tempmc,tempres,tempc,temptru,tempsub,tempx,tempy,
     &  tempz,temprad,tempB) 

      character*(*) templine
      character*8 tempstring
      character*6 tempah
      character*4 tempanam
      character*3 tempres
      character*1 tempsub
      character*1 tempc
      character*1 tempmc
      integer     tempanum
      integer     temptru
      real        tempx
      real        tempy
      real        tempz
      real        temprad
      real        tempB

      integer StringLength
      integer linelen

      linelen = StringLength(templine)
      if (linelen.lt.50) then
       print *,'Error parsing ATOM record'
       print *,'linelen = ',linelen
       stop
      end if

      tempah = templine(1:6)
      tempstring = templine(7:11) 
      call strtoint(tempstring,tempanum)
      tempanam = templine(13:16)
      tempmc   = templine(17:17)
      tempres  = templine(18:20)
      tempc    = templine(22:22)
      tempstring = templine(23:26)
      call strtoint(tempstring,temptru)
      tempsub  = templine(27:27)
      tempstring = templine(31:38)
      call strtoreal(tempstring,tempx)
      tempstring = templine(39:46)
      call strtoreal(tempstring,tempy)
      tempstring = templine(47:54)
      call strtoreal(tempstring,tempz)

      if (linelen.ge.60) then
       tempstring = templine(57:60)
       call strtoreal(tempstring,temprad)
      else
       temprad = 1.00
      end if
 
      if (linelen.ge.65) then
       tempstring = templine(61:65)
       call strtoreal(tempstring,tempB)
      else 
       tempB = 0.00
      end if
  
      end
c-----------------------------------------------------------
      subroutine pdbwrite (filename)

      include 'parameters.h'
      include 'pdb.h'
      character*(*) filename
      integer i,j

9     format (a80)
10    format (a6,i5,1x,a4,1x,a3,1x,a1,i4,a1,3x,f8.3,f8.3,
     +    f8.3,2x,f4.2,f6.2)
11    format (a6,1x,i4,1x,i4,1x,i4,1x,i4,1x,i4)

      open (unit=1,file=filename)

      do i = 1,numcomlines
       write (1,9) cmnt(i)
      end do
      do j = 1,numatm
      write (1,10) atheader(j),atomnum(j),atomname(j),
     +  res(j),chain(j),trueResnum(j),subnum(j),coord(1,j),
     +  coord(2,j),coord(3,j),rad(j),Bfac(j)
      end do

      close (1)

      end 
c-----------------------------------------------------------
      subroutine pdbwrite_w_ter (filename)

      include 'parameters.h'
      include 'pdb.h'
      character*(*) filename
      integer i,j

9     format (a80)
10    format (a6,i5,1x,a4,1x,a3,1x,a1,i4,a1,3x,f8.3,f8.3,
     +    f8.3,2x,f4.2,f6.2)
11    format (a6,1x,i4,1x,i4,1x,i4,1x,i4,1x,i4)

      open (unit=1,file=filename)

      do i = 1,numcomlines
       write (1,9) cmnt(i)
      end do
      do j = 1,numatm
       if (j.ne.1) then
        if (chain(j).ne.chain(j-1).or.atomname(j-1).eq.' OXT') then
          write(1,10) 'TER   '
        end if
         write (1,10) atheader(j),atomnum(j),atomname(j),
     +     res(j),chain(j),trueResnum(j),subnum(j),coord(1,j),
     +     coord(2,j),coord(3,j),rad(j),Bfac(j)
       else
       write (1,10) atheader(j),atomnum(j),atomname(j),
     +  res(j),chain(j),trueResnum(j),subnum(j),coord(1,j),
     +  coord(2,j),coord(3,j),rad(j),Bfac(j)
       end if
      end do
      write(1,10) 'TER   '

      close (1)

      end 
c-------------------------------------------------------------
      subroutine pdbwrite_fragments (filename,fragments,
     &    numfragments,rvchain)

      include 'parameters.h'
      include 'pdb.h'
      character*(*) filename
      integer i,j
      integer numfragments
      character*1 rvchain
      integer fragments(20,2)

9     format (a80)
10    format (a6,i5,1x,a4,1x,a3,1x,a1,i4,a1,3x,f8.3,f8.3,
     +    f8.3,2x,f4.2,f6.2)
11    format (a6,1x,i4,1x,i4,1x,i4,1x,i4,1x,i4)

      open (unit=1,file=filename)

c      do i = 1,numcomlines
c       write (1,9) cmnt(i)
c      end do

      do i = 1,numfragments
        do j = 1,numatm
          if (resnum(j).ge.fragments(i,1).and.
     &        resnum(j).le.fragments(i,2).and.
     &        chain(j).eq.rvchain) then
            write (1,10) atheader(j),atomnum(j),atomname(j),
     &         res(j),chain(j),trueResnum(j),subnum(j),coord(1,j),
     &         coord(2,j),coord(3,j),rad(j),Bfac(j)
          end if
        end do
      end do

      close (1)

      end 
c-------------------------------------------------------------
      subroutine transform (matrix1,matrix2,oldCoord,newCoord)

c Expects to be few a 3x3 transformation matrix and a 1x3 translation
c matrix
c
c Author: CM Summa
c Last Modified; 6/25/98

      real matrix1(3,3)
      real matrix2(3)
      real oldCoord(3)
      real newCoord(3)
      real tx
      real ty
      real tz

      tx = (matrix1(1,1)*oldCoord(1)) + (matrix1(1,2)*oldCoord(2)) +
     &     (matrix1(1,3)*oldCoord(3))
      ty = (matrix1(2,1)*oldCoord(1)) + (matrix1(2,2)*oldCoord(2)) +
     &     (matrix1(2,3)*oldCoord(3))
      tz = (matrix1(3,1)*oldCoord(1)) + (matrix1(3,2)*oldCoord(2)) +
     &     (matrix1(3,3)*oldcoord(3))
     
       newCoord(1) = tx + matrix2(1)
       newCoord(2) = ty + matrix2(2)
       newCoord(3) = tz + matrix2(3)

      end

c---------------------------------------------------------------------
      subroutine pdbcheck (contactdir,pdbname,reason,firsttime)

      include 'parameters.h'
      include 'pdb.h'
      
      integer  AtomsInRes(maxNumRes)
      integer  IResArray(maxNumRes)
      integer  rescounter
      integer  atomcounter
      integer  example(maxNumRes)
      character*3   ResArray(maxNumRes)
      character*20  reason
      character*(*)  contactdir
      character*(*)  pdbname
      character*120  concat
      logical  protein(maxNumRes)
      logical  badAtom
      logical  firsttime
      logical  exist

12    format (i2)
13    format (a3)
14    format (a4)
15    format (a20,1x,a3,a1,a1,i4,1x,a33)

      if (.not.readStdAminos) then
       call fileexist('prot.std',contactdir,concat,exist)
       if (.not.exist) then
        print *,'Can''t find prot.std in appropriate place'
        print *,'Stopping'
        stop
       end if
 
       open (unit=1,file=concat,status='OLD')
 
       do j=1,20
        read (1,12) atomsInAA(j)
        read (1,13) AAname(j)
        do i=1,atomsInAA(j)
         read (1,14) AtName(j,i)
        end do
       end do
       close (1)
      end if

      badAtom = .false.

c      print *,'numres = ',numres
      count1 = 1
      count2 = 1
      count3 = 1
      do i=1,maxNumRes
       AtomsInRes(i) = 0
       protein(i) = .false.
      end do
  

      rescounter = 1
      atomcounter = 1
      do i=1,numatm
       if (.not.unknown(i)) then
       if (atheader(i)(1:4).eq.'ATOM') then
        do j=1,20
         if(res(i).eq.AAname(j)) then
          protein(rescounter) = .true.
         end if
        end do
       end if
  

        if (protein(rescounter)) then
       
c code to force reject if there are gaps in the chain....
c doesn't work when heteroatoms are read in (not currently used!)
         goto 125
         if (rescounter.gt.1) then
          if (protein(rescounter-1)) then
           if (trueResnum(i).ne.trueResnum(i-1).and.
     &       trueResnum(i).ne.trueResnum(i-1)+1.and.
     &       chain(i).eq.chain(i-1)) then
            print *,'GAP between atoms',atomnum(i-1),
     &        atomnum(i),'  ',pdbname(1:11)
            reject = .true.
            reason(1:5) = 'GAP  '
           end if
          end if
         end if
125      continue

c        print *,res(i),rescounter,atomcounter,pdbname
        ResArray(rescounter) = res(i)
c        AtomArray(rescounter,atomcounter) = atomname(i)
        AtomsInRes(rescounter) = AtomsInRes(rescounter) + 1
        if (resnum(i+1).ne.resnum(i)) then
c          print *,'AtomsInRes = ',AtomsInRes(rescounter)
          example(rescounter) = i
          rescounter = rescounter + 1
          atomcounter = 1
        else
          atomcounter = atomcounter + 1
        end if
       end if
       end if
      end do

150   do i=1,20
       do j=1,numres
        if(AAname(i).eq.ResArray(j)) then
         IResArray(j) = i
        end if
       end do
      end do

200   do i=1,numres
c       print *,i,protein(i),IResArray(i),
c     &  AtomsInRes(i),trueResnum(i)
       if (protein(i).and.IResArray(i).ne.0) then
        if (AtomsInRes(i).ne.atomsInAA(IResArray(i))) then
         if ((AtomsInRes(i).ne.atomsInAA(IResArray(i))+1).and.
     &      (i.ne.numres)) then
          if (AtomsInRes(i).lt.atomsInAA(IResArray(i))) then
c            if (.not.reject) then
c             temptext = pdbname(5:7)
c             call toUpperCase(temptext,3)
c             print 15,'Color Molecule Atoms',temptext,':',
c     &       chain(example(i)),trueResnum(example(i)),
c     &       'Specified Specification 255,127,0'
c            end if
c            reject = .true.
c            reason(6:10) = 'MISAT'
          end if
          if (AtomsInRes(i).gt.atomsInAA(IResArray(i))) then
           if (.not.firsttime) then
            print *,'too many atoms in res ',i,'  ',pdbname(1:11)
           end if
           badAtom = .true.
          end if
         end if
        end if
       end if
      end do

c      print *,'calling oligomercheck'
c      call oligomercheck(cmnt,numcomlines,numchains,reject,reason)
c      print *,'finished oligomercheck'

999   end
c--------------------------------------------------------------
      subroutine findRelevantChain (relchainstr)

      include 'parameters.h'
      include 'pdb.h'
      character*1 relchainstr
      logical startFound
      logical endFound

      startFound = .false.
      endFound = .false.

c Case where relchainstr is 'defined'
      if (relchainstr.ne.' ') then
 
        i = 1
         do while ((.not.startFound).and.i.lt.numatm)
          if (chain(i).eq.relchainstr.and.
     &        atheader(i)(1:4).eq.'ATOM') then
           startFound = .true.
           chainStart = i
          end if
          i = i + 1
         end do 

        do while (.not.endFound.and.i.le.numatm)
         if (chain(i).ne.relchainstr.or.
     &     atheader(i)(1:4).ne.'ATOM') then
           endFound = .true.
           chainEnd = i - 1
         end if
         if (i.eq.numatm.and.(.not.endFound)) then
          endFound = .true.
          chainEnd = i
         end if
c         print *,i,relchainstr,chain(i),endFound
         i = i + 1
        end do

        if (startFound.and.endFound) then
         return
        else
         print *,'Error - chain start or end not found!'
         stop
        end if
      end if

c Case where relchainstr is 'undefined'
      if (relchainstr.eq.' ') then
       print *,'Null chainstring...'
        i = 1
        do while (.not.startFound)
          if (atheader(i)(1:4).eq.'ATOM') then
           startFound = .true.
           chainStart = i
          end if
          i = i + 1
         end do

        do while (.not.endFound.and.i.le.numatm)
         if ((atheader(i)(1:4).ne.'ATOM').or.
     &       (atheader(i)(1:4).eq.'ATOM'.and.
     &        chain(i).ne.relchainstr)) then
           endFound = .true.
           chainEnd = i - 1
         end if
         if (i.eq.numatm.and.(.not.endFound)) then
          endFound = .true.
          chainEnd = i
         end if 
         i = i + 1
        end do
        if (startFound.and.endFound) then
         return
        else
         print *,'Error - chain start or end not found!'
         stop
        end if
      end if

      end
    
c--------------------------------------------------------- 
      subroutine pdbcollapse(category)


      include 'parameters.h'
      include 'pdb.h'

      logical collapse(maxat)
      integer category(maxat,3)
      integer i,j,k

      do i=1,maxat
       collapse(i) = .false.
      end do

      do i=1,numatm
c collapse out all unknown atoms
c       if (unknown(i)) then
c        collapse(i) = .true.
c       end if

        if (res(i+1).eq.res(i).and.
     &   atomname(i+1).eq.atomname(i).and.
     &   trueResnum(i+1).eq.trueResnum(i)) then
          if (rad(i).ge.rad(i+1)) then
            collapse(i) = .true.
          else
            collapse(i+1) = .true.
          end if
        end if

      end do

c      do i=1,numatm
c      print *,i,unknown(i),collapse(i)
c      end do

      i = 1
      do while (i.le.numatm)
       if (collapse(i)) then
        atheader(i) = atheader(i+1)
        atomname(i) = atomname(i+1)
        res(i) = res(i+1)
        subnum(i) = subnum(i+1)
        chain(i) = chain(i+1)
        multiconf(i) = multiconf(i+1)
        atomnum(i) = atomnum(i+1)
        resnum(i) = resnum(i+1)
        trueResnum(i) = trueResnum(i+1)
        do j=1,3
         coord(j,i) = coord(j,i+1)
         category(i,j) = category(i+1,j)
        end do
        rad(i) = rad(i+1)
        Bfac(i) = Bfac(i+1)
        unknown(i) = unknown(i+1)
        collapse(i) = collapse(i+1)
        do k =i+1,numatm
         atheader(k) = atheader(k+1)
         atomname(k) = atomname(k+1)
         res(k) = res(k+1)
         subnum(k) = subnum(k+1)
         chain(k) = chain(k+1)
         multiconf(k) = multiconf(k+1)
         atomnum(k) = atomnum(k+1)
         resnum(k) = resnum(k+1)
         trueResnum(k) = trueResnum(k+1)
         do j=1,3
          coord(j,k) = coord(j,k+1)
          category(k,j) = category(k+1,j)
         end do
         rad(k) = rad(k+1)
         Bfac(k) = Bfac(k+1)
         unknown(k) = unknown(k+1)
         collapse(k) = collapse(k+1)
        end do
        numatm = numatm - 1
       else
        i = i + 1
       end if
      end do 

      end

c--------------------------------------------------------

      function Distance(x1,y1,z1,x2,y2,z2)

      real x1,y1,z1,x2,y2,z2,Distance

      Distance = sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)

      end

c--------------------------------------------------------           

      subroutine initialize_pdb_data()

      include 'parameters.h'
      include 'pdb.h'

      integer i,j

      do i=1,maxnumcomments
       cmnt(i) = ''
      end do

      do i=1,maxat
       atheader(i) = ''
       atomname(i) = ''
       res(i)      = ''
       subnum(i)   = ''
       chain(i)    = ''
       multiconf(i)= ''
       atomnum(i) = 0
       resnum(i) = 0
       rad(i)  = 0.0
       Bfac(i) = 0.0
       unknown(i) = .false.
       hetatm(i) = .false.
       do j=1,5
        connectivity(j,i) = 0
       end do
       do j=1,3
        coord(j,i) = 0.0
       end do
      end do

      numatm = 0
      chainStart = 0
      chainEnd = 0
      numres = 0
      numchains = 0
      relevantChain = 0
      numcomlines = 0
      reject = .true.

      end

c----------------------------------------------------------
      function testMatrixEquality(x,y)

      logical testMatrixEquality

      real x(3,3),y(3,3)
      integer i,j,k

      testMatrixEquality = .true.
      do i=1,3
       do j=1,3
        if (x(i,j).ne.y(i,j)) then
          testMatrixEquality = .false.
          return
        end if
       end do
      end do
      end
       
c----------------------------------------------------------
      function testVectorEquality(x,y)

      logical testVectorEquality

      real x(3), y(3)
      integer i,j,k

      testVectorEquality = .true.
      do i=1,3
        if (x(i).ne.y(i)) then
          testVectorEquality = .false.
          return
        end if
      end do
       
      end
