
      subroutine BubbleSort(energy,best,bestmap,worst,worstmap,
     & numatomtypes,comboTypes,sig)

      include 'parameters.h'

      real energy(maxnumatomtypes,maxNumCoord,570)
      real best(maxnumatomtypes,20)
      real worst(maxnumatomtypes,20)
      logical sig(maxnumatomtypes,maxNumCoord,570)
      integer numatomtypes
      integer comboTypes(maxNumCoord)
      integer bestmap(20,20,2)
      integer worstmap(20,20,2)

c-----------------------------------------------------------
c Initialize the best and worst arrays with dummy values
c which will all be scrapped after the sort has run. Since
c a good energy (in our case) is defined by a negative value,
c and a bad energy is defined as one having a positive value,
c the bubble sort is 'backwards' from what we are normally
c used to seeing.....

      do i=1,numatomtypes
       do j=1,20
        best(i,j) = 5.0
        worst(i,j) = -5.0
        do k=1,2
         bestmap(i,j,k) = 0
         worstmap(i,j,k) = 0
        end do
       end do
      end do 

c-------------------------------------------------------------
c Actual sort begins here

      do c1 = 1,numatomtypes
       do c2 = 1,maxNumCoord
        do c3 = 1,comboTypes(c2)

         if ((energy(c1,c2,c3).le.best(c1,20)).and.
     &       (sig(c1,c2,c3)).and.
     &       (energy(c1,c2,c3).ne.33.0)) then
          best(c1,20) = energy(c1,c2,c3)
          bestmap(c1,20,1) = c2
          bestmap(c1,20,2) = c3
          do c4 = 20,2,-1
           if (best(c1,c4).le.best(c1,c4-1)) then
            call rswap(best(c1,c4),best(c1,c4-1))
            call iswap(bestmap(c1,c4,1),bestmap(c1,c4-1,1))
            call iswap(bestmap(c1,c4,2),bestmap(c1,c4-1,2))
           end if
          end do 
         end if

         if ((energy(c1,c2,c3).ge.worst(c1,20)).and.
     &       (sig(c1,c2,c3)).and.
     &       (energy(c1,c2,c3).ne.33.0)) then
          worst(c1,20) = energy(c1,c2,c3)
          worstmap(c1,20,1) = c2
          worstmap(c1,20,2) = c3
          do c4 = 20,2,-1
           if (worst(c1,c4).ge.worst(c1,c4-1)) then
            call rswap(worst(c1,c4),worst(c1,c4-1))
            call iswap(worstmap(c1,c4,1),worstmap(c1,c4-1,1))
            call iswap(worstmap(c1,c4,2),worstmap(c1,c4-1,2))
           end if
          end do
         end if

        end do
       end do
      end do

      end
c----------------------------------------------------------------
      subroutine bsort(realArray,sortedMap,num)

      real realArray(*)
      real tempArray(100000)
      integer sortedMap(*)
      integer num

      do i=1,num
       tempArray(i) = realArray(i)
       sortedMap(i) = i
      end do

      do i=1,num-1
       do j=i+1,num
        if (tempArray(i).lt.tempArray(j)) then
         call rswap(tempArray(i),tempArray(j))
         call iswap(sortedMap(i),sortedMap(j))
        end if
       end do
      end do

      end

c---------------------------------------------------------------
      subroutine iswap (a,b)
      integer a,b,c

      c = a
      a = b
      b = c

      end
c----------------------------------------------------------------
      subroutine rswap (a,b)
      real a,b,c

      c = a
      a = b
      b = c

      end
