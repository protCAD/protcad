      subroutine gettheinteger(istr,str)

      character*80 str
        character ch
        integer numvar,istr,ich,power
        integer begin(10),end(10)

        call parseline(str,begin,end,numvar)
        istr = 0
        if(numvar.gt.2)then
           write(6,*)'Error parsing  argument in gettheinteger()'
      else
           do i = end(1),begin(1),-1
                ch = str(i:i)
                ich = ichar(ch) - 48
                power = 10**(end(1)-i)
                istr = istr + power*ich
         enddo
      endif

        return
        end
c
c *************************************************************
c
        subroutine parseline(longstring,begin,end,numvar)

        character*80 longstring
        character*1 char
        integer mxwordsinstring
        parameter (mxwordsinstring=10)
        integer begin(mxwordsinstring),end(mxwordsinstring)
        integer begincount,endcount,numvar

        begincount = 0
        endcount = 0
        do i = 1,80
           char = longstring(i:i)
           if(begincount .eq. endcount)then
c     -     I am looking for the beginning of an entry.
                if(char.ne.' ')then
                   begincount = begincount + 1
                   begin(begincount) = i
            endif
         elseif (begincount .eq. (endcount + 1))then
c     -     I am looking for the end of an entry.
                if(char.eq.' ')then
                   endcount = endcount + 1
                   end(endcount) = i - 1
            else
                   if(i.eq.80)then
                        endcount = endcount + 1
                        end(endcount) = i
               endif
            endif
         else
                write(6,*)'Error in subroutine parseline().'
                goto 20
         endif
      enddo

        numvar = begincount

        return
 20   end
c
c ****************************************************************
c
      subroutine getthereal(rstr,str)

      integer mxcmdlinevar
        parameter (mxcmdlinevar=10)
        integer begin(mxcmdlinevar),end(mxcmdlinevar)
        integer ich,power,decimalptindx,tempint,numvar
        character*80 str
        character ch,chch
        real parity
        real rstr,leftofdecimal,rightofdecimal

        call parseline(str,begin,end,numvar)

c     Begin by checking if the argument is negative.
      i = begin(1)
      chch = str(i:i)
        if(chch.eq.'-')then
           parity = -1.0
           begin(1) = begin(1) + 1
      else
           parity = 1.0
      endif
c
        decimalptindx = end(1) + 1
        if(numvar.gt.2)then
           write(6,*)'Error parsing argument in getthereal()'
      else
           do i = begin(1),end(1)
                ch = str(i:i)
                ich = ichar(ch)
                if(ich.eq.46)then   ! 46 is ascii code for the decimal point
                   decimalptindx = i
            endif
         enddo
         tempint = 0
           if(decimalptindx.gt.begin(1))then
                ii = decimalptindx - 1
              do i = ii,begin(1),-1
                   power = 10**(ii - i)
                   ch = str(i:i)
                   ich = ichar(ch) - 48
                   tempint = tempint + ich*power
            enddo
              leftofdecimal = float(tempint)
         else
                leftofdecimal = 0.0
         endif
           tempint = 0
           ii = decimalptindx + 1
           if(decimalptindx.lt.end(1))then
            do i = end(1),ii,-1
                   power = 10**(end(1)-i)
                   ch = str(i:i)
                   ich = ichar(ch) - 48
                   tempint = tempint + ich*power
            enddo
                rightofdecimal = float(tempint)/float(10*power)
         else
                rightofdecimal = 0.0
         endif
           rstr = leftofdecimal + rightofdecimal
      endif
        rstr = parity*rstr

      return
        end
