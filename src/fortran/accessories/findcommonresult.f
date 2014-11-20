       program findcommonresult

       parameter (bignum = 1000000)

       character*5 junk
       character*11 pdbname(bignum)
       character*11 temppdbname
       character*3  res(bignum)
       character*3  tempres
       character*3  atomname(bignum)
       character*3  tempatomname
       character*6  type(bignum)
       character*6  temptype
       character*20 tempfilename
       character*20 listname
       integer      timesseen(bignum)
       real         energy(bignum)
 
       real tempenergy
       integer attype
       integer tempattype
       integer rank
       integer temprank
       integer coord
       integer tempcoord
       integer resnum(bignum)
       integer tempresnum
       integer numarguments

       integer targettype
       integer targetrank
       character targetjason
       character*80 targetstr
       character*1  tempchain
       character*1 chain(bignum)
       character*1 hit(bignum,20)
       integer      listnum
       
       integer numculled
       logical found
       
10     format (a11,3x,a5,2x,f6.3,1x,i2,i3,i3,1x,a3,1x,i5,1x,a1,3x,a3)
11     format (a3)
12     format (a20)
13     format (i4)
14     format (a1) 
15     format (a11,3x,a3,2x,i4,2x,a1,2x,a3,2x,a6,2x,i3,
     &  a1,a1,a1,a1,a1,a1,a1,a1,a1,a1,a1,a1,a1,
     &  a1,a1,a1,a1,a1,a1,a1,a1) 

      numarguments = iargc()

      if (numarguments.ne.2) then
        print *,'USAGE: findcommonresult filelist atomtype'
        goto 999
      endif

      call getarg(1,listname)
      call getarg(2,targetstr)
      call gettheinteger(targettype,targetstr) 

       targetrank = 20
       numculled = 0
       do i=1,bignum
        timesseen(i) = 0
        do j=1,20
         hit(i,j) = ' '
        end do 
       end do
       listnum = 1 

c0     print *,'Input filenamelist'
c      read (5,12) listname

       open (unit = 2, file = listname,err = 999) 
30     read (2,12,end=900 ) tempfilename

       print *,tempfilename
       open (unit=1,file=tempfilename)

100    read (1,11,end = 800) junk
       if (junk.ne.'pdb') then
        goto 100
       end if

       backspace(1)

200    read (1,10,end = 800) temppdbname,temptype,tempenergy,
     &   tempattype,temprank,
     &   tempcoord,tempres,tempresnum,tempchain,tempatomname

c       print 10, temppdbname,temptype,tempenergy,
c     &   tempattype,temprank,
c     &   tempcoord,tempres,tempresnum,tempchain,tempatomname


       if (tempattype.eq.targettype.and.temprank.le.targetrank) then
c check to see if you've seen this already
c        print *,'found one'
        found = .false.
        if (numculled.ne.0) then
         do i=1,numculled 
          if (tempres.eq.res(i).and.tempresnum.eq.resnum(i).and.
     &        temppdbname.eq.pdbname(i).and.tempatomname.eq.
     &        atomname(i).and.tempchain.eq.chain(i).and.
     &        temptype.eq.type(i)) then
            timesseen(i)  = timesseen(i) + 1
            hit(i,listnum) = 'x'
c            print *,'found a multiple', timesseen(i)
            found = .true.
          end if
         end do
        else
         numculled = numculled + 1
         res(numculled) = tempres
         resnum(numculled) = tempresnum
         pdbname(numculled) = temppdbname
         atomname(numculled)= tempatomname
         chain(numculled) = tempchain
         type(numculled) = temptype
         timesseen(numculled) = 1 
         hit(numculled,listnum) =  'x'
c         print *,'entered first one ',numculled
        end if  
        if (.not.found) then
c         print *,'found a new one',numculled+1
         numculled = numculled + 1
         if (numculled.eq.bignum) then
          print *,'Please increase the value of bignum'
          stop
         end if
         res(numculled) = tempres
         resnum(numculled) = tempresnum
         pdbname(numculled) = temppdbname
         atomname(numculled)= tempatomname
         chain(numculled) = tempchain
         type(numculled) = temptype
         timesseen(numculled) = 1
         hit(numculled,listnum) = 'x'
        end if
       end if

       goto 200

800    continue
       close (1) 
       listnum =listnum+1
       goto 30

900    continue
       close (2)
  
       listnum = listnum-1

       print *,'Atoms which are commonly in Best or Worst Clusters'
       do i=1,numculled
c        if (timesseen(i).gt.3) then
          print 15,pdbname(i),res(i),
     &      resnum(i),chain(i),
     &      atomname(i),type(i),timesseen(i),
     &      ' ','|',(hit(i,j),j=1,listnum)
c        end if
       end do

999    end
c **************************************************************
c
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
