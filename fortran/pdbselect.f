c-------------------------------------------------------------------------------
      subroutine Read_PDBselect(dbconcat,pdbnames,pdbChain,numpdbnames,
     &    nHlx,nBta,thrsh)

c      Author: Chris Summa
c      Last Modified: 5/28/98

c     Calls subroutine toLowerCase
c     Calls subroutine toUpperCase

      include 'parameters.h'
      include 'pdbselect.h'

       character*25 tempFtxt
       character*120 dbconcat
       character*5  tempID
       character*1  tempMethd
       integer      tempThrsh
       integer      tempNAA
       real         tempResolution
       real         tempRfac
       integer      tempNSid
       integer      tempNBck
       integer      tempNNaa
       integer      tempNHlx
       integer      tempNBta

      character*13 junk1
      character*5  junk2
      integer      j
      integer      goodRes,goodRfac,Great
      logical      keep
      

10    format (a13)
11    format (a5)
20    format (i5,3x,a5,2x,i4,2x,f4.2,2x,f4.2,5x,
     &  a1,3x,i4,3x,i4,5x,i2,3x,i4,3x,i4,a25)

      goodRfac = 0
      goodRes = 0
      Great = 0

      open (unit = 2,file=dbconcat,status='OLD')

c Start reading and eating lines until you 
c get to the actual data

100   read (2,10,end=999) junk1
      if (junk1.ne.'thrsh      ID') then
       goto 100
      end if

c Check to see that the next line actually
c contains data

      j= 1
150   read (2,11,end=999) junk2
      if (junk2.ne.'   25'.and.junk2.ne.'   35'.and.
     &  junk2.ne.'   45'.and.junk2.ne.'   30'.and.
     &  junk2.ne.'   95') then
       goto 150
      end if
 
      backspace(2) 
 

200   read (2,20,end = 999) tempThrsh,tempID,
     &  tempNAA,tempResolution,tempRfac,
     &  tempMethd,tempNSid,tempNBck,
     &  tempNNaa,tempNHlx,tempNBta,tempFtxt


     
      call pdbSelectFilter(tempNAA,tempResolution,tempRfac,
     &  tempMethd,tempNSid,tempNBck,tempNNaa,keep)

      if (keep) then
       thrsh(j) = tempThrsh
       ID(j) = tempID
       naa(j) = tempNAA
       resolution(j) = tempResolution
       Rfac(j) = tempRfac
       Methd(j) = tempMethd
       nSid(j) = tempNSid
       nBck(j) = tempNBck
       nNAA(j) = tempNNaa
       nHlx(j) = tempNHlx
       nBta(j) = tempNBta
       ftxt(j) = tempFtxt
       call toLowerCase(ID(j),5)
       j = j + 1
      end if


      goto 150

999    continue
       numpdbnames = j-1

       close (2)
      
      do i=1,numpdbnames
       if (ID(i)(5:5).eq.'_') then
        pdbChain(i) = ' '
       else
        pdbChain(i) = ID(i)(5:5)
        call toUpperCase(pdbChain(i),1)
       end if
       pdbnames(i)(1:3) = 'pdb'
       do j=4,7
        pdbnames(i)(j:j) = ID(i)(j-3:j-3)
       end do
       pdbnames(i)(8:11) = '.ent'
      end do 

c      do i=1,numpdbnames
c       print *,pdbnames(i),pdbChain(i)
c      end do

c      print *,'goodRes = ',goodRes
c      print *,'goodRfac = ',goodRfac 
c      print *,'Great = ',Great


      end

      subroutine colpdblist(rejectflag,pdbnames,pdbChain,numpdbnames,
     &   nHlx,nBta,thrsh)

      include 'parameters.h'
 
      character*11 pdbnames(maxNumPDBs)
      character*1  pdbChain(maxNumPDBs)
      integer      thrsh(maxNumPDBs)
      integer      nHlx(maxNumPDBs)
      integer      nBta(maxNumPDBs)
      integer      numpdbnames
      integer i,k
      logical rejectflag(maxNumPDBs)

      i = 1
      do while (i.le.numpdbnames)
       if (rejectflag(i)) then
        print *,'Rejecting ',pdbnames(i)
        pdbnames(i) = pdbnames(i+1)
        pdbChain(i) = pdbChain(i+1)
        nHlx(i) = nHlx(i+1)
        nBta(i) = nBta(i+1)
        thrsh(i) = thrsh(i+1)
        rejectflag(i) = rejectflag(i+1)
        do k=i+1,numpdbnames
         pdbnames(k) = pdbnames(k+1)
         pdbChain(k) = pdbChain(k+1)
         nHlx(k) = nHlx(k+1)
         nBta(k) = nBta(k+1)
         thrsh(k) = thrsh(k+1)
         rejectflag(k) = rejectflag(k+1)
        end do
        numpdbnames = numpdbnames -1
       else
        i = i + 1
       end if
      end do

      end
c---------------------------------------------------------------
      subroutine pdbSelectFilter (tempNAA,tempResolution,tempRfac,
     &  tempMethd,tempNSid,tempNBck,tempNNAA,keep)

       character*1  tempMethd
       integer      tempNAA
       real         tempResolution
       real         tempRfac
       integer      tempNSid
       integer      tempNBck
       integer      tempNNAA
       logical      keep


       keep = .true.

       if (tempMethd.eq.'N') then
        keep = .false.
        return
       end if

       if (tempNAA.lt.50) then
        keep = .false.
        return
       end if

       if (tempResolution.gt.2.60) then
        keep = .false.
        return
       end if

       if (tempRfac.gt.0.22) then
        keep = .false.
        return
       end if

       if (tempNSid.ne.tempNBck) then
        keep = .false.
        return
       end if
 
       if (tempNNAA.ne.0) then
        keep = .false.
        return
       end if

       end
c------------------------------------------------
      subroutine Read_pdblist(dbconcat)

c      Author: Chris Summa

      include 'parameters.h'
      include 'encad.h'

      character*25 tempFtxt
      character*120 dbconcat
      integer      j
      

10    format (a20)

      open (unit = 2,file=dbconcat,status='OLD')

c Start reading and eating lines until you 
c get to the actual data

      j = 1
100   read (2,10,end=800) encadpdbnames(j)
      j = j + 1
      goto 100
800   close(2)
      numencadpdbnames = j-1
      end
