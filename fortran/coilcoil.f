

c Make a coiled coil from a bundle of straight alpha-helices 
c 
c assume bundle is oriented along the z axis 
c
      parameter (natmx=10000,natmx3=3*natmx,nresmx=600,nresp1=nresmx+1) 
c
      character*4 atnam(natmx)
      character*3 resnam(nresmx)
      integer fatres(nresp1)
      real x(natmx3),chg(natmx),sth,cth,xold,yold,zmin,twopi
      integer ia(natmx3),ib(natmx3)
      character*80 infil,outfil
      real cm(3)
      real*8 v(3,3),scr(36000) 
c
      twopi=2.*3.14159265
 10   write(6,99)
 99   format(' Name of PDB file - ',$)
      read (5,98) infil
 98   format(a80)
      call pdbin1(infil,10,nat,nres,nbond,atnam,resnam,fatres,
     .           x,chg,ia,ib,ierr)
      if(ierr.gt.0.and.ierr.lt.4) then
        write(6,*) ' error ',ierr,' reading file ',infil
        go to 10
      end if
      if(nat.gt.natmx) then
        write(6,*) ' Too many atoms - stop'
        stop
      end if 
c 
c   molecule should be oriented on the z-axis 
c
      xmin=-99.e10
      do i=1,3*nat,3
        if(x(i+2).lt.zmin) zmin=x(i+2)
      end do 
c
 30   write(6,97)
 97   format(' Supercoil pitch - ',$)
      read(5,*) pitch 
c      write(6,96) 
c 96   format(' Supercoil radius - ',$) 
c      read(5,*) radius 
c 
c   Now do the actual twist 
c
      do i=1,3*nat,3 
c   calculate the amount to twist this residue 
C   make theta negative for a left-handed twist 
c
 40     theta=(-(x(i+2)-zmin)/pitch)*twopi
        sth=sin(theta)
        cth=cos(theta) 
c 
c   twist each atom in this residue by theta 
c   about the z-axis and then translate appropriately 
c
        xold=x(i)
        yold=x(i+1)
        x(i)= xold*cth-yold*sth
        x(i+1)=xold*sth+yold*cth
        x(i+2)=x(i+2)
      end do 
c
      write(6,94)
 94   format(' File for output coordinates - ',$)
      read(5,98) outfil
      call pdbot1(outfil,12,nat,nres,0,atnam,resnam,fatres,
     .            x,chg,ia,ib)
      stop
      end
      subroutine pdbin1(fname,lu,nat,nres,nbond,atnam,resnam,fatres,
     .                 x,chg,ia,ib,ierr) 
c 
c   general purpose .pdb file reading program 
c   at the moment, it assumes the residues are in order 
c                  and reads only the first modecule 
c 
c   fname - the name of the .PDB file
      character*(*) fname 
c   nat - number of atoms
c   nres - number of residues
c   atnam(i) - name of the i-th atom
      character*4 atnam(1) 
c   resnam(i) - name of the i-th residue
      character*3 resnam(1) 
c   fatres(i) - atom number for first atom of the i-th residue
      integer fatres(1) 
c   x(i) - coordinates
      real x(1) 
c   chg(i) - charge of the i-th atom
      real chg(1) 
c   ia(i) and ib(i) are the atoms of the i-th bond
      integer ia(1),ib(1) 
c   ierr returned as non-zero on error 
c        = 1 file not found 
c          2 no ATOM cards 
c          3 residues not in order 
c 
c   temporary variables
      character*5 atemp
      character*3 rtemp
      character*6 label
      character*74 string
      integer ibtmp(8) 
c 
c   open the .PDB file
      open(lu,file=fname,status='old') 
c
      iat=0
      kat=-2
      ires=0
      fatres(1)=1
      ksav=0 
c
 10   read(lu,99,end=510) label,string
 99   format(a6,a74)
      if(label.ne.'ATOM ') go to 10 
c 
c   found the ATOM cards 
c
 20   iat=iat+1
      kat=kat+3
      read(string,98) atemp,rtemp,kres,x(kat),x(kat+1),x(kat+2),
     .                chg(iat)
 98   format(6x,a5,a3,1x,i5,4x,3f8.3,f7.3)
      if(atemp(1:1).eq.' ') then
        atnam(iat)=atemp(2:5)
      else
        atnam(iat)=atemp(1:4)
      end if 
c     if(kres.ne.ires.and.kres.ne.ires+1.and.ires.ne.0) go to 520
      if(kres.ne.ires.and.kres.ne.ksav) then
        ksav=kres
        ires=ires+1
        resnam(ires)=rtemp
        fatres(ires)=iat
      end if 
c
      read(lu,99,end=100) label,string
      if(label.eq.'ATOM ') go to 20 
c 
c   finished reading
 100  nat=iat
      nres=ires
      fatres(nres+1)=nat+1 
c 
c   now pick up the connectivity 
c
      nbond=0
 110  read(lu,99,end=530) label,string
      if(label.ne.'CONECT') go to 110 
c
 120  read(string,97) iatmp,(ibtmp(i),i=1,8)
 97   format(10i5)
      if(iatmp.gt.nat) go to 200
      do i=1,8
        if(ibtmp(i).gt.iatmp) then
          nbond=nbond+1
          ia(nbond)=iatmp
          ib(nbond)=ibtmp(i)
        else if (ibtmp(i).eq.0) then
          go to 130
        end if
      end do 
c
 130  read(lu,99,end=200) label,string
      if(label.eq.'CONECT') go to 120 
c
 200  close (lu)
      return 
c 
c   file not found - set ierr to 1
 500  ierr=1
      return 
c 
c   no ATOM cards - set ierr to 2
 510  ierr=2
      close(lu)
      return 
c 
c   residues not in order - set ierr to 3
 520  ierr=3
      close(lu)
      return 
c 
c   no connectivity records - set ierr to 4 
c
 530  ierr=4
      return
      end
      subroutine pdbot1(fname,lu,nat,nres,nbond,atnam,resnam,fatres,
     .                  x,chg,ia,ib) 
c
      character*(*) fname
      integer fatres(1)
      character*3 resnam(1)
      dimension x(1),chg(1)
      dimension ia(1),ib(1)
      character*4 atnam(1)
      dimension list(8) 
c
      open(lu,file=fname,status='unknown')
      do i=1,nres 
c       if(resnam(i).eq.'LIN') write(lu,98)
        do j=fatres(i),fatres(i+1)-1
          jj=3*(j-1)
          write(lu,99) j,atnam(j),resnam(i),i,x(jj+1),x(jj+2),x(jj+3),chg(j)
 99       format('ATOM',2x,i5,2x,a4,a3,2x,i4,4x,3f8.3,f6.3)
        end do
      end do
      write(lu,98)
 98   format('TER')
      if(nbond.gt.0) then
        do i=1,nat
          nlist=0
          do j=1,nbond
            if(ia(j).gt.0) then
              if (ia(j).eq.i) then
                nlist=nlist+1
                list(nlist)=ib(j)
              else if (ib(j).eq.i) then
                nlist=nlist+1
                list(nlist)=ia(j)
                ia(j)=-1 
              end if
            end if
          end do
 10       if(nlist.gt.0) then
            write(lu,97) i,(list(k),k=1,nlist)
 97         format('CONECT',9i5)
          end if
        end do
      end if
      write(lu,96)
 96   format('END')
      return
      end

