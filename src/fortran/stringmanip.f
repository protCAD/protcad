      subroutine toLowerCase(txt,len)
c
c convert character string to lower case
c
      character*(*) txt
      character*80 save
      character*26 ualpha,lalpha
      integer StringLength
      data ualpha /'ABCDEFGHIJKLMNOPQRSTUVWXYZ'/
      data lalpha /'abcdefghijklmnopqrstuvwxyz'/

      do i=1,StringLength(txt)
        if((txt(i:i).ge.'A').and.(txt(i:i).le.'Z')) then
          match = index(ualpha,txt(i:i))
          save(i:i) = lalpha(match:match)
        else
          save(i:i) = txt(i:i)
        end if
      end do 

      txt = save
      return
      end


c------------------------------------------------------------
      subroutine toUpperCase(txt,len)
c
c convert character string to upper case
c
      character*(*) txt
      character*80 save
      character*26 ualpha,lalpha
      integer StringLength
      data ualpha /'ABCDEFGHIJKLMNOPQRSTUVWXYZ'/
      data lalpha /'abcdefghijklmnopqrstuvwxyz'/

      do i=1,StringLength(txt)
        if((txt(i:i).ge.'a').and.(txt(i:i).le.'z')) then
          match = index(lalpha,txt(i:i))
          save(i:i) = ualpha(match:match)
        else
          save(i:i) = txt(i:i)
        end if
      end do 

      txt = save
      return
      end

c--------------------------------------------------------------
       subroutine searchString (source,string,found,where)

       integer where
       integer slength
       integer srclength
       integer j
       integer StringLength
       character*(*) string
       character*(*) source
       logical found

       srclength = StringLength(source)
       slength = StringLength(string)

       found = .false.
        do j=1,srclength-slength
         if (source(j:j+slength).eq.string) then
           found = .true.
           where = j
           return
         end if
        end do

       where = 999 

999    end
c--------------------------------------------------------------------
       function StringLength (string)

       integer StringLength
       integer maxlength
       integer place
       character*(*) string

       maxlength = len(string)
       
       do place = maxlength,1,-1
        if (string(place:place).ne.'') then
         StringLength = place 
         return
        end if
       end do
       end
c--------------------------------------------------------------------
       subroutine removeTrailingSpaces (string)

       integer LineLength
       integer maxlength
       integer place
       logical isitAlpha
       character*(*) string

       maxlength = len(string)
       
       do place = maxlength,1,-1
        if (string(place:place).eq.' ') then
          string(place:place) = ''
        else if (isitAlpha(string(place:place))) then
         return
        end if
       end do
       end
c---------------------------------------------------------------------
       function nextChar (chr)

       character chr
       character nextChar
       integer tempint

        tempint = ichar(chr)+1
        nextChar = char(tempint)

       end
c--------------------------------------------------------------------
       function incChar (chr,i)
       character incChar
       character chr
       integer i
       integer tempint

        tempint = ichar(chr)
        if (tempint+i.le.126) then
         tempint = tempint + i
         incChar = char(tempint)
         return
        else
         incChar = char(33 + (tempint+i-126))
        end if
       end
c------------------------------------------------------------------
       function isitAlpha(ch)

       character ch
       logical isitAlpha
       integer x

       isitAlpha = .false.
       x = ichar(ch)
       if ((x.ge.48.and.x.le.57).or.
     &     (x.ge.65.and.x.le.90).or.
     &     (x.ge.97.and.x.le.122)) then
         isitAlpha = .true.
         return
        end if

        end

c------------------------------------------------------------------
        subroutine findNextWord(string,eolw,begin,finish,
     &   nextword)

        character *(*) string
        integer eolw
        integer begin
        integer finish
        integer length
        integer StringLength
        logical startfound
        logical finishfound
        logical nextword


        length = StringLength(string)
        startfound = .false.
        finishfound = .false.
        nextword = .false.

        do i=eolw,length
         x = ichar(string(i:i))
         if ((x.ge.33.and.x.le.126).and.(.not.startfound)) then
          startfound = .true.
          begin = i
         end if
         if ((x.le.32.or.x.eq.127).and.(startfound)) then
          finish = i-1
          nextword= .true.
          return
         end if
        end do

        if (startfound) then
         finish = length
         nextword = .true.
        end if
        end

c **************************************************************
c
      subroutine strtoint(str,istr)

      character*(*) str
      character ch
      integer istr,ich,power
      integer begin,finish
      logical exist

      call findNextWord(str,1,begin,finish,exist)

      if (.not.exist) then
       print *,'error! - i cant convert a null string to an int'
       stop 
      end if 

      istr = 0
      do i = finish,begin,-1
       ch = str(i:i)
       ich = ichar(ch) - 48
       power = 10**(finish-i)
       istr = istr + power * ich
      end do

      end

c ****************************************************************
c
      subroutine strtoreal(str,rstr)

      integer begin,finish
      integer ich,power,decimalptindx,tempint
      character*80 str
      character ch,chch
      real parity
      real rstr,leftofdecimal,rightofdecimal
      logical exist

      call findNextWord(str,1,begin,finish,exist)
      if (.not.exist) then
       print *,'error - i cant convert a null string to real!'
       stop
      end if

c     Begin by checking if the argument is negative.
      i = begin
      chch = str(i:i)
      if(chch.eq.'-')then
       parity = -1.0
       begin = begin + 1
      else
       parity = 1.0
      endif

      decimalptindx = finish + 1
      do i = begin,finish
       ch = str(i:i)
       ich = ichar(ch)
       if(ich.eq.46)then   ! 46 is ascii code for the decimal point
        decimalptindx = i
       endif
      enddo

      tempint = 0

      if(decimalptindx.gt.begin)then
       ii = decimalptindx - 1
       do i = ii,begin,-1
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
       if(decimalptindx.lt.finish)then
        do i = finish,ii,-1
         power = 10**(finish-i)
         ch = str(i:i)
         ich = ichar(ch) - 48
         tempint = tempint + ich*power
        enddo
        rightofdecimal = float(tempint)/float(10*power)
       else
        rightofdecimal = 0.0
       endif
       rstr = leftofdecimal + rightofdecimal
       rstr = parity*rstr

       end
c----------------------------------------------------------
c234567
       function compareIgnoringWhitespace(string1, string2)

       logical compareIgnoringWhitespace
       character*(*) string1
       character*(*) string2
       integer begin1
       integer begin2
       integer finish1
       integer finish2
       integer length1
       integer length2
       logical nextword1
       logical nextword2

       call findNextWord(string1,1,begin1,finish1,
     &   nextword1)
       call findNextWord(string2,1,begin2,finish2,
     &   nextword2)

       compareIgnoringWhitespace = .false.
       if (.not.nextword1) then
         return
       endif
       if (.not.nextword2) then
         return
       endif

       if (string1(begin1:finish1).eq.string2(begin2:finish2)) then
          compareIgnoringWhitespace = .true.
	  return
       end if
       return
       end
