      subroutine Make_Bin_Decision(tempdistance,maxdistance,
     &    mindistance,distanceStepsize,whereToBin)

      real tempdistance,maxdistance,mindistance,distanceStepsize,
     +     templowerlimit,tempupperlimit
      integer whereToBin

      templowerlimit = mindistance

      if (tempdistance.le.mindistance) then
       goto 500
      end if

      whereToBin = 1
      do tempupperlimit = mindistance+distanceStepsize,maxdistance,
     +   distanceStepsize
       if (tempdistance.gt.templowerlimit.and.tempdistance.le.
     +     tempupperlimit) then
         return
       end if
       whereToBin = whereToBin + 1
       templowerlimit = tempupperlimit
      end do

c------------------------
c If subroutine reaches this point, interaction falls
c outside the bounds of the specified max and min distances.
c A value of 999 acts a flag to tell the calling routine that
c this is a special case.

500   whereToBin = 999

      end
