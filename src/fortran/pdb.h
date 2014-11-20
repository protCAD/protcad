      character*100 cmnt(maxnumcomments)
      character*6 atheader(maxat)
      character*4 atomname(maxat)
      character*3 res(maxat)
      character*1 subnum(maxat)
      character*1 chain(maxat)
      character*1 multiconf(maxat)
      integer numatm
      integer chainStart
      integer chainEnd
      integer atomnum(maxat)
      integer resnum(maxat)
      integer trueResnum(maxat)
      integer numres
      integer numchains
      integer connectivity(5,maxat)
      integer relevantChain
      integer numcomlines
      real coord(3,maxat)
      real rad(maxat)
      real Bfac(maxat)
      logical unknown(maxat)
      logical reject
      logical hetatm(maxat)
      logical readStdAminos

      character*3   AAname(20)
      character*4   AtName(20,20)
      integer  atomsInAA(20)

      common / pdbblock / cmnt,atheader,atomname,res,subnum,
     &   chain,multiconf,numatm,chainStart,chainEnd,atomnum,
     &   resnum,trueResnum,numres,numchains,connectivity,
     &   relevantChain,numcomlines,coord,rad,Bfac,unknown,reject,
     &   hetatm,readStdAminos,atomsInAA,AAname,AtName
