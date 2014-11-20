
      subroutine relevantDistance(i1,i2,numskip,relevant)
c     & numskip,residuename1,residuename2,atomname1,atomname2,
c     & chain1,chain2,resnum,relevant)

      include 'parameters.h'
      include 'pdb.h'
      integer i1,i2,numskip
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

      if (chain(i1).ne.chain(i2)) then
       sameChain = .false.
      end if

	if (unknown(i1)) then
	 relevant = .false.
        return
       end if

       if (unknown(i2)) then
        relevant = .false.
        return
       end if

      if (numskip.eq.0) then
       relevant = .true.
       return
      end if

      if (.not.sameChain) then
       relevant = .true.
       return
      end if

c--Case where resnum(i2) is less than resnum(i1)

      if (resnum(i2).lt.resnum(i1)) then
       if (resnum(i2).ge.resnum(i1)-numskip) then
         relevant = .false.
         return
       end if      

       if (resnum(i1)-resnum(i2).eq.1.and.atomname(i1).eq.
     +  ' N  '.and.atomname(i2).eq.' C  ') then
          relevant = .false.
          return
       end if
      end if

c--Case where resnum(i2) is greater than resnum(i1)

      if (resnum(i2).ge.resnum(i1)) then
       if (resnum(i2).le.resnum(i1)+numskip) then
          relevant = .false.
          return
       end if
       if (resnum(i2)-resnum(i1).eq.1.and.atomname(i1).eq.
     +   ' C  '.and.atomname(i2).eq.' N  ') then
         relevant = .false.
         return
       end if
      end if

      end

c--------------------------------------------------------------

      subroutine relevantDistance2(i1,i2,numskip,relevant)
c     & numskip,residuename1,residuename2,atomname1,atomname2,
c     & chain1,chain2,resnum,relevant)

      include 'parameters.h'
      include 'pdb.h'
      integer i1,i2,numskip
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

      if (chain(i1).ne.chain(i2)) then
       sameChain = .false.
      end if

      if (unknown(i1)) then
        relevant = .false.
        return
      end if

      if (unknown(i2)) then
       relevant = .false.
       return
      end if

      if (numskip.eq.0) then
       relevant = .true.
       return
      end if

      if (.not.sameChain) then
       relevant = .true.
       return
      end if

c--Case where resnum(i2) is less than resnum(i1)

      if (resnum(i2).lt.resnum(i1)) then
       if (resnum(i2).gt.resnum(i1)-numskip) then
         relevant = .false.
         return
       end if      

       if (resnum(i1)-resnum(i2).eq.1.and.atomname(i1).eq.
     &  ' N  '.and.atomname(i2).eq.' C  '
     &  .and.numskip.eq.1) then
          relevant = .false.
          return
       end if
      end if

c--Case where resnum(i2) is greater than resnum(i1)

      if (resnum(i2).ge.resnum(i1)) then
       if (resnum(i2).lt.resnum(i1)+numskip) then
          relevant = .false.
          return
       end if
       if (resnum(i2)-resnum(i1).eq.1.and.atomname(i1).eq.
     &   ' C  '.and.atomname(i2).eq.' N  '
     &   .and.numskip.eq.1) then
         relevant = .false.
         return
       end if
      end if

      end
