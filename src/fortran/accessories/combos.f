      program combos
      parameter (numtypes = 4)
      parameter (numitems = 5)
      parameter (maxnumunique = 10000)

      integer unique(maxnumunique,numtypes)
      logical first
      logical seen
      logical identical
      integer numunique
      integer icount(numitems)
      integer jcount(numtypes)
      integer i1
      integer i2
      integer icount1,icount2,icount3,icount4,icount5,icount6
      integer icount7,icount8

      do i1 = 1,maxnumunique
       do i2 = 1,numtypes
        unique(i1,i2) = 0
       end do
      end do
  
  
      first = .true.
      do icount1 = 1,numtypes
       icount(1) = icount1  
       do icount2 = 1,numtypes
        icount(2) = icount2
        do icount3 = 1,numtypes
         icount(3) = icount3
         do icount4 = 1,numtypes
          icount(4) = icount4
          do icount5 = 1,numtypes
           icount(5) = icount5


              do i1 = 1,numtypes
               jcount(i1) = 0
              end do
 
              do i1 = 1,numitems
               do j1 = 1,numtypes
                if (icount(i1).eq.j1) then
                 jcount(j1) = jcount(j1) + 1
                end if
               end do
              end do

 
              if (first) then 
               numunique = 1
               do i1 = 1,numtypes
                unique(numunique,i1) = jcount(i1)
               end do
               first = .false.
              else
               seen = .false.
               do i1 = 1,numunique
                identical = .true.
                do i2 = 1,numtypes
                 if (jcount(i2).ne.unique(i1,i2)) then 
                    identical = .false.
                 end if
                end do
                if (identical) then
                  seen = .true.
                end if
               end do
               if (.not.seen) then
                numunique = numunique + 1
                do i1 = 1,numtypes
                  unique(numunique,i1) = jcount(i1)
                end do
               end if
              end if
 
          end do
         end do
        end do
       end do
      end do

      print *,'number of unique combos is: ',numunique
      do i1 = 1,numunique
       print *,(unique(i1,i2),i2=1,numtypes)
      end do
      end

