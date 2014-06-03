      subroutine correlate(x,y,n,r)

c x and y are 1D corresponding arrays of type real (maximum of 'bignum' units)
c n is an integer value designating the actual length of the x and y arrays
c r is the returned correlation coefficient

      parameter (bignum=500)
      integer i,n
      real x(*)
      real y(*)
      real r
      real xdiff(bignum)
      real ydiff(bignum)
      real*8 xbar,ybar,numerator,xdenom,ydenom,denominator

c     initialize varibles

      do i=1,n
       xdiff(i) = 0.0
       ydiff(i) = 0.0
      end do
      r = 0.0
      xbar = 0.0
      ybar = 0.0
      xdenom = 0.0
      ydenom = 0.0
      numerator = 0.0
      denominator = 0.0

      do i=1,n
       xbar = xbar + x(i)
       ybar = ybar + y(i)
      end do
      xbar = xbar / n
      ybar = ybar / n

      do i=1,n
       xdiff(i) = x(i) - xbar
       ydiff(i) = y(i) - ybar
      end do

      do i=1,n
       numerator = numerator + (xdiff(i) * ydiff(i))
       xdenom = xdenom + (xdiff(i)**2)
       ydenom = ydenom + (ydiff(i)**2)
      end do
      denominator = sqrt(xdenom * ydenom)
      r = numerator/denominator
      end
