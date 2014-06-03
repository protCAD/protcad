      integer encadnumentries
      character*1 encadresname(1000)
      character*3 encadatomname(1000)
      character*2 encadatomtype(1000)
      integer     encadatomtypeint(1000)

      integer encadcategorymapsize
      character*2 encadcategorymapchar(40)

      character*20 encadpdbnames(10000)
      integer numencadpdbnames

      common / encadblock / encadresname,encadatomname,encadatomtype,
     &  encadnumentries, encadcategorymapsize, encadcategorymapchar,
     &  encadatomtypeint, encadpdbnames, numencadpdbnames
