      subroutine random_subset(numatm,numToSample,iseed,subset)

c calls iswap,randm,ran


      include 'parameters.h'

      real randm
      real numatmReal
      integer numatm
      integer iseed
      integer j
      integer i
      integer numToSample
      integer subset(maxat)
      integer tempint
      logical found

      numatmReal = float(numatm)

      do i=1,numToSample
       found = .true.
       do while (found)
        found = .false.
        tempint = anint(randm(iseed,numatmReal))
        do j=1,i-1
         if (tempint.eq.subset(j)) then
          found = .true.
         end if
        end do
       end do
       subset(i) = tempint
      end do

      do i=1,numToSample-1
       do j=i,numToSample
        if (subset(i).gt.subset(j)) then
        call iswap(subset(i),subset(j))
        end if
       end do
      end do

c      print *,'ordered subset'
c      do i=1,numToSample
c       print *,subset(i)
c      end do

      end

c
c *********************************
c
      function randm(iseed,max)

      integer iseed
       real max,randm

       randm = (max*ran1(iseed))+0.5

       return
       end
c
c
c *********************************
c
       function ran1(iseed)
       real ran1
         iseed = iseed*69069 + 1
         ran1 = float(ishft(iseed,-8))*0.5**24
       return
       end
