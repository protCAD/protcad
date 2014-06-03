FFLAGS_DEBUG = -C -mips4 -g2 
FFLAGS_R10 = -O2 -mips4
FFLAGS_OPT = -col72 -O2 -OPT:Olimit=0
FFLAGS_OPT2 = -col120 -O2 
FFLAGS_LINUX_DEBUG = -Wall -ansi -C -g
FFLAGS_LINUX_OPT_386 = -O2 -ansi -malign-double -funroll-loops
FFLAGS_LINUX_OPT_PPC = -O2 -ansi -funroll-loops

F77 = f77
#FORTFLAGS = $(FFLAGS_PPC_DEBUG)
#FORTFLAGS = $(FFLAGS_R10)
#FORTFLAGS  = $(FFLAGS_OPT2)
#FORTFLAGS = $(FFLAGS_LINUX_DEBUG)
#FORTFLAGS = $(FFLAGS_LINUX_OPT_386)
#FORTFLAGS = $(FFLAGS_LINUX_OPT_386)


BIN = /Users/csumma/bin

.f.o:
	$(F77)  -c $(FORTFLAGS)  $*.f

clean: 
	rm ./*.o


blastfilter: blastfilter.o pdbread.o
	$(FCOMP) -o blastfilter blastfilter.o pdbread.o
	mv blastfilter $(BIN)

sumscoredata: sumscoredata.o stringmanip.o cmdline.o
	$(F77) $(FORTFLAGS) -o sumscoredata sumscoredata.o \
	stringmanip.o cmdline.o

	mv sumscoredata $(BIN)

deriveMicroEnv: deriveMicroEnv.o pdb.o significance.o bubbleSort.o \
	pdbselect.o stringmanip.o categorizeatoms.o countcategories.o \
	files.o PMFmodule.o random_subset.o funx.o cmdline.o \
	relevant.o correlate.o makeBinDecision.o readCombos.o
	$(F77) $(FORTFLAGS) -o deriveMicroEnv deriveMicroEnv.o \
	pdb.o significance.o bubbleSort.o pdbselect.o \
	stringmanip.o categorizeatoms.o countcategories.o \
	files.o PMFmodule.o random_subset.o funx.o cmdline.o \
	relevant.o correlate.o makeBinDecision.o readCombos.o

	mv deriveMicroEnv $(BIN)

deriveMicroEnvDiff: deriveMicroEnvDiff.o pdb.o significance.o bubbleSort.o \
	pdbselect.o stringmanip.o categorizeatoms.o countcategories.o \
	files.o funx.o cmdline.o readCombos.o\
	relevant.o correlate.o makeBinDecision.o
	$(F77) $(FORTFLAGS) -o deriveMicroEnvDiff deriveMicroEnvDiff.o \
	pdb.o significance.o bubbleSort.o pdbselect.o \
	stringmanip.o categorizeatoms.o countcategories.o \
	files.o funx.o cmdline.o \
	relevant.o correlate.o makeBinDecision.o readCombos.o

	mv deriveMicroEnvDiff $(BIN)

calculateHiResPMF: calculateHiResPMF.o pdb.o significance.o bubbleSort.o \
	pdbselect.o stringmanip.o categorizeatoms.o countcategories.o \
	files.o random_subset.o funx.o cmdline.o \
	relevant.o correlate.o PMFmodule.o
	$(F77) $(FORTFLAGS) -o calculateHiResPMF calculateHiResPMF.o \
	pdb.o significance.o bubbleSort.o pdbselect.o \
	stringmanip.o categorizeatoms.o countcategories.o \
	files.o random_subset.o funx.o cmdline.o \
	relevant.o correlate.o PMFmodule.o

	mv calculateHiResPMF $(BIN)

calculateLoResPMF: calculateLoResPMF.o pdb.o significance.o bubbleSort.o \
	pdbselect.o stringmanip.o categorizeatoms.o countcategories.o \
	files.o random_subset.o funx.o cmdline.o \
	relevant.o correlate.o PMFmodule.o
	$(F77) $(FORTFLAGS) -o calculateLoResPMF calculateLoResPMF.o \
	pdb.o significance.o bubbleSort.o pdbselect.o \
	stringmanip.o categorizeatoms.o countcategories.o \
	files.o random_subset.o funx.o cmdline.o \
	relevant.o correlate.o PMFmodule.o

	mv calculateLoResPMF $(BIN)

calculateHiResPMF_vol: calculateHiResPMF_vol.o pdb.o significance.o bubbleSort.o \
	pdbselect.o stringmanip.o categorizeatoms.o countcategories.o \
	files.o random_subset.o funx.o cmdline.o \
	relevant.o correlate.o PMFmodule_vol.o
	$(F77) $(FORTFLAGS) -o calculateHiResPMF_vol calculateHiResPMF_vol.o \
	pdb.o significance.o bubbleSort.o pdbselect.o \
	stringmanip.o categorizeatoms.o countcategories.o \
	files.o random_subset.o funx.o cmdline.o \
	relevant.o correlate.o PMFmodule_vol.o

	mv calculateHiResPMF_vol $(BIN)

gen_biomolecule: gen_biomolecule.o pdb.o stringmanip.o files.o cmdline.o
	$(F77) $(FORTFLAGS) -o gen_biomolecule gen_biomolecule.o \
	pdb.o stringmanip.o files.o cmdline.o

	mv gen_biomolecule $(BIN)

cleanpdb: cleanpdb.o pdb.o stringmanip.o files.o cmdline.o
	$(F77) $(FORTFLAGS) -o cleanpdb cleanpdb.o \
	pdb.o stringmanip.o files.o cmdline.o

	mv cleanpdb $(BIN)

commonfrags: commonfrags.o pdb.o \
	stringmanip.o \
	files.o cmdline.o relevant.o
	$(F77) $(FORTFLAGS) -o commonfrags commonfrags.o \
	pdb.o \
	stringmanip.o \
	files.o cmdline.o relevant.o

	mv commonfrags $(BIN)

calculateEncadPMF: calculateEncadPMF.o pdb.o significance.o bubbleSort.o \
	pdbselect.o stringmanip.o categorizeatoms.o countcategories.o \
	files.o random_subset.o funx.o cmdline.o \
	relevant.o correlate.o PMFmodule.o read_encad_namelist.o
	$(F77) $(FORTFLAGS) -o calculateEncadPMF calculateEncadPMF.o \
	pdb.o significance.o bubbleSort.o pdbselect.o \
	stringmanip.o categorizeatoms.o countcategories.o \
	files.o random_subset.o funx.o cmdline.o \
	relevant.o correlate.o PMFmodule.o read_encad_namelist.o

	mv calculateEncadPMF $(BIN)
