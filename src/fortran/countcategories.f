      subroutine countcategories (category,numatomsInCategory)

      include 'parameters.h'
      include 'pdb.h'

      real numatomsInCategory(maxnumatomtypes,3)
      integer  category(maxat,3)
      integer  i

      do i=1,maxnumatomtypes
       do j=1,3
        numatomsInCategory(i,j) = 0.0
       end do
      end do
   
      print *,'chainStart = ',atomnum(chainStart)
      print *,'chainEnd = ',atomnum(chainEnd)
      
      do i=chainStart,chainEnd
       if (category(i,1).ne.0.and.category(i,2).ne.0.and.
     &     category(i,3).ne.0) then
        numatomsInCategory(category(i,1),1) = 
     &   numatomsInCategory(category(i,1),1) + 1.0
        numatomsInCategory(category(i,2),2) = 
     &   numatomsInCategory(category(i,2),2) + 1.0
        numatomsInCategory(category(i,3),3) =
     &   numatomsInCategory(category(i,3),3) + 1.0
       end if
      end do

      end 
