      subroutine read_encad_namelist()
      implicit none
      
      include 'encad.h'

      character*80 filename
      integer counter

10    format(a1,a3,3x,a2)

      filename = './Full_list_encad.dat'
      open(unit =1, file=filename,err=180)
      encadnumentries = 1
50    read(1,10,end=200) encadresname(encadnumentries),
     & encadatomname(encadnumentries),encadatomtype(encadnumentries)
c      print *,'read a line'
      encadnumentries = encadnumentries + 1
      goto 50

180   print *, 'Error opening file named ',filename
200   close (1)
      encadnumentries = encadnumentries -1
      end

c----------------------------------------------------------------
      function convertFromSingleLetterCode( slc )
      implicit none

      character*1 slc
      character*3 convertFromSingleLetterCode

      
      if (slc.eq.'.') then
        convertFromSingleLetterCode = '***'
	return
      end if
      if (slc.eq.'A') then
        convertFromSingleLetterCode = 'ALA'
	return
      end if
      if (slc.eq.'C') then
        convertFromSingleLetterCode = 'CYS'
	return
      end if
      if (slc.eq.'D') then
        convertFromSingleLetterCode = 'ASP'
	return
      end if
      if (slc.eq.'E') then
        convertFromSingleLetterCode = 'GLU'
	return
      end if
      if (slc.eq.'F') then
        convertFromSingleLetterCode = 'PHE'
	return
      end if
      if (slc.eq.'G') then
        convertFromSingleLetterCode = 'GLY'
	return
      end if
      if (slc.eq.'H') then
        convertFromSingleLetterCode = 'HIS'
	return
      end if
      if (slc.eq.'I') then
        convertFromSingleLetterCode = 'ILE'
	return
      end if
      if (slc.eq.'K') then
        convertFromSingleLetterCode = 'LYS'
	return
      end if
      if (slc.eq.'L') then
        convertFromSingleLetterCode = 'LEU'
	return
      end if
      if (slc.eq.'M') then
        convertFromSingleLetterCode = 'MET'
	return
      end if
      if (slc.eq.'N') then
        convertFromSingleLetterCode = 'ASN'
	return
      end if
      if (slc.eq.'P') then
        convertFromSingleLetterCode = 'PRO'
	return
      end if
      if (slc.eq.'Q') then
        convertFromSingleLetterCode = 'GLN'
	return
      end if
      if (slc.eq.'R') then
        convertFromSingleLetterCode = 'ARG'
	return
      end if
      if (slc.eq.'S') then
        convertFromSingleLetterCode = 'SER'
	return
      end if
      if (slc.eq.'T') then
        convertFromSingleLetterCode = 'THR'
	return
      end if
      if (slc.eq.'V') then
        convertFromSingleLetterCode = 'VAL'
	return
      end if
      if (slc.eq.'W') then
        convertFromSingleLetterCode = 'TRP'
	return
      end if
      if (slc.eq.'Y') then
        convertFromSingleLetterCode = 'TYR'
	return
      end if
      end

c----------------------------------------------------------------
      function convertToSingleLetterCode( name )
      implicit none

      character*3 name
      character*1 convertToSingleLetterCode

      
      if (name.eq.'***') then
        convertToSingleLetterCode = '.'
	return
      end if
      if (name.eq.'ALA') then
        convertToSingleLetterCode = 'A'
	return
      end if
      if (name.eq.'CYS') then
        convertToSingleLetterCode = 'C'
	return
      end if
      if (name.eq.'ASP') then
        convertToSingleLetterCode = 'D'
	return
      end if
      if (name.eq.'GLU') then
        convertToSingleLetterCode = 'E'
	return
      end if
      if (name.eq.'PHE') then
        convertToSingleLetterCode = 'F'
	return
      end if
      if (name.eq.'GLY') then
        convertToSingleLetterCode = 'G'
	return
      end if
      if (name.eq.'HIS') then
        convertToSingleLetterCode = 'H'
	return
      end if
      if (name.eq.'ILE') then
        convertToSingleLetterCode = 'I'
	return
      end if
      if (name.eq.'LYS') then
        convertToSingleLetterCode = 'K'
	return
      end if
      if (name.eq.'LEU') then
        convertToSingleLetterCode = 'L'
	return
      end if
      if (name.eq.'MET') then
        convertToSingleLetterCode = 'M'
	return
      end if
      if (name.eq.'ASN') then
        convertToSingleLetterCode = 'N'
	return
      end if
      if (name.eq.'PRO') then
        convertToSingleLetterCode = 'P'
	return
      end if
      if (name.eq.'GLN') then
        convertToSingleLetterCode = 'Q'
	return
      end if
      if (name.eq.'ARG') then
        convertToSingleLetterCode = 'R'
	return
      end if
      if (name.eq.'SER') then
        convertToSingleLetterCode = 'S'
	return
      end if
      if (name.eq.'THR') then
        convertToSingleLetterCode = 'T'
	return
      end if
      if (name.eq.'VAL') then
        convertToSingleLetterCode = 'V'
	return
      end if
      if (name.eq.'TRP') then
        convertToSingleLetterCode = 'W'
	return
      end if
      if (name.eq.'TYR') then
        convertToSingleLetterCode = 'Y'
	return
      end if
      end
c----------------------------------------------------
c234567
      subroutine generateCategoryMap()
      implicit none
      include 'encad.h'
      
      integer i,j,k
      integer currentnumentries
      logical found
     
      currentnumentries = 0
      do i=1,encadnumentries
        if (currentnumentries.eq.0) then
	  currentnumentries = currentnumentries + 1
	  encadcategorymapchar(currentnumentries) =
     &      encadatomtype(currentnumentries)
        else
	  found = .false.
	  do j=1,currentnumentries
	    if (encadcategorymapchar(j).eq.encadatomtype(i)) then
	      found = .true.
	    end if
	  end do
	  if (.not.found) then
	    currentnumentries = currentnumentries + 1
	    encadcategorymapchar(currentnumentries) =
     &        encadatomtype(i)
          end if
	end if
       end do

       encadcategorymapsize = currentnumentries

c       do i=1,encadcategorymapsize
c         print *, i, encadcategorymapchar(i)
c       end do

       do i=1,encadnumentries
        do j=1,encadcategorymapsize
	 if (encadatomtype(i).eq.encadcategorymapchar(j)) then
	   encadatomtypeint(i) = j
	 end if
	end do
       end do

       end
