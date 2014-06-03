FFLAGS_DEBUG = -C -g2 
FFLAGS_R10 = -O2 -mips4
FFLAGS_OPT = -col72 -O2 -OPT:Olimit=0
FFLAGS_OPT2 = -O2 
FFLAGS_PPC_DEBUG = -Wall -ansi -C -g
FFLAGS_PPC_OPT = -O2 -ansi

F77 = f77
#FORTFLAGS = $(FFLAGS_PPC_DEBUG)
#FORTFLAGS = $(FFLAGS_R10)
FORTFLAGS  = $(FFLAGS_OPT2)
#FORTFLAGS = $(FFLAGS_DEBUG)

BIN = /home/summa/bin

.f.o:
	$(F77)  -c $(FORTFLAGS)  $*.f

con10_pair: con10_pair.o pdb.o significance.o bubbleSort.o \
	pdbselect.o stringmanip.o categorizeatoms.o countcategories.o \
	files.o PMFmodule.o random_subset.o funx.f
	$(F77) $(FORTFLAGS) -o con10_pair con10_pair.o \
	pdb.o significance.o bubbleSort.o pdbselect.o \
	stringmanip.o categorizeatoms.o countcategories.o \
	files.o PMFmodule.o random_subset.o funx.f

	mv con10_pair $(BIN)

con10_std: con10_std.o pdb.o significance.o bubbleSort.o \
	pdbselect.o stringmanip.o categorizeatoms.o countcategories.o \
	files.o PMFmodule.o random_subset.o funx.f
	$(F77) $(FORTFLAGS) -o con10_std con10_std.o \
	pdb.o significance.o bubbleSort.o pdbselect.o \
	stringmanip.o categorizeatoms.o countcategories.o \
	files.o PMFmodule.o random_subset.o funx.f

	mv con10_std $(BIN)

cadis_jae: cadis.o pdb.o funx.o stringmanip.o files.o
	$(F77) $(FORTFLAGS) -o cadis_jae cadis.o pdb.o funx.o stringmanip.o files.o

readresult: readresult.o
	f77 $(FORTFLAGS) -o readresult readresult.o
	mv readresult $(BIN)

findcommonresult: findcommonresult.o
	f77 $(FORTFLAGS) -o findcommonresult findcommonresult.o
	mv findcommonresult $(BIN)

