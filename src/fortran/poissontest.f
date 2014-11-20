      program poissontest

      real Poisson
      real expected
      real true
      real sig
      
      true = 1.0
      do expected = 0.0,5.0,0.1
       sig = 1-Poisson(true,expected)
       print *,expected,sig
      end do 
      end

c---------------------------------------
      function Poisson(x1,mu)

      real x1
      real mu
      real Poisson
      double precision Fact

      Poisson = (exp(-mu)*(mu**x1))/Fact(int(x1))
      end
c--------------------------------------
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

