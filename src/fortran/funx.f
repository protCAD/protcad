      function func1(middle,range,value)

      real func1
      real middle
      real range
      real value
      real tvalue
      real transform
      real slope

      transform = middle - range
      tvalue = value - transform
      slope = -1/(range*2)

      func1 = slope*tvalue + 1

      end

c--------------------------------
      function Fact(foobar)

      double precision Fact
      integer i,foobar

      Fact = 1
      if (foobar.ge.2) then
       do i=2,foobar
        Fact = Fact*i
       end do
      end if
      end
