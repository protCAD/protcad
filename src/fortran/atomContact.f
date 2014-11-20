      subroutine atomContact(dist,type1,type2,bool)

      include 'parameters.h'
      include 'pdb.h'
      integer i1,i2
      integer type1,type2
      logical bool
      real distance
      real    radius(4)

c------------------------------------------------------
c    type 1 = N
c    type 2 = C
c    type 3 = O
c    type 4 = S
c------------------------------------------------------
c   set up matrix

      radius(1) = 3.5
      radius(2) = 3.7
      radius(3) = 3.2
      radius(4) = 4.0

      bool = .false.
      distance = ((radius(type1) + radius(type2))/2) + 0.3
      if (dist.lt.distance) then
        bool = .true.
        return
      end if

      end
