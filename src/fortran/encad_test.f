c234567
      program testit
      implicit none
      include 'encad.h'
      character*3 convertFromSingleLetterCode
      integer i
      logical compareIgnoringWhitespace
      logical flag


      call read_encad_namelist()
      call generateCategoryMap()
      do i=1,encadnumentries
        print *, convertFromSingleLetterCode(encadresname(i)),
     &  encadatomname(i), encadatomtype(i), encadatomtypeint(i)
      end do

      flag = compareIgnoringWhitespace(" H ", "H  ")
      print *,flag      
      flag = compareIgnoringWhitespace(" H2", "H2")
      print *,flag      
      flag = compareIgnoringWhitespace("3H ", "  3H  ")
      print *,flag      
      flag = compareIgnoringWhitespace("1H ", "H1  ")
      print *,flag      

      end
