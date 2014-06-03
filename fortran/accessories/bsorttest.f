      program bsorttest

      real array(10)
      integer map(10)
      integer num

      num = 10

      array(1) = 5
      array(2) = 11
      array(3) = 100
      array(4) = 1
      array(5) = 16
      array(6) = 4
      array(7) = 33
      array(8) = 15
      array(9) = 200
      array(10) = 3
  
      call bsort(array,map,num)
      end

c----------------------------------------------------------------
      subroutine bsort(realArray,sortedMap,num)

      real realArray(*)
      real tempArray(100000)
      real target
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

      do i=1,num
       print *,tempArray(i),sortedMap(i)
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
