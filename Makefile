export PROTCADDIR=$(PWD)
export SRCDIR=$(PROTCADDIR)/src
export TNTINCLUDE=$(SRCDIR)
export BINDIR=$(PROTCADDIR)/bin
export UIDIR=$(PROTCADDIR)/ui
export OBJDIR=$(PROTCADDIR)/obj
export PROJDIR=$(PROTCADDIR)/projects

UNAME := $(shell uname -s)
ifeq ($(UNAME),Linux)
	export $(PATH)=$(PATH):$(PROTCADDIR):$(BINDIR)
endif
ifeq ($(UNAME),Darwin)	
	export $(PATH)=$(PATH):$(PROTCADDIR):$(BINDIR):/usr/local/gfortran
endif


export CXX = g++
export F77 = gfortran
export MAKE = make

SHELL = /bin/sh

TARGETS = protDielectric protEvolver protDihedrals protOligamer protEnergy protFolder protMover protMutator protSequence protInverter protMin protAlign protShaper protBindingEnergy hammingdist

.SUFFIXES: .cc .o .h .a .f

LIB_TARGETS = lib

LIB_CC_OBJECTS = ran.o point.o treeNode.o atom.o atomIterator.o residue.o chain.o residueTemplate.o allowedResidue.o secondaryStructure.o chainPosition.o residueIterator.o chainModBuffer.o molecule.o protein.o ensemble.o CMath.o generalio.o pdbData.o pdbReader.o pdbWriter.o amberVDW.o aaBaseline.o amberElec.o rotamer.o rotamerLib.o PDBAtomRecord.o PDBInterface.o ruler.o line.o lineSegment.o unitSphere.o helixPropensity.o parse.o ramachandranMap.o 

LIB_F77_OBJECTS = bestfit.o

DEFS = -D__STL_USE_EXCEPTIONS

FLAG_OPTMAX = -Wall -oFast -ffast-math -ftree-vectorize -march=native -mtune=native -pipe -msse3 -Wno-deprecated -std=gnu++11

CFLAGS = $(FLAG_OPTMAX) $(DEFS)

FFLAGS = -Wall -g -Wno-tabs -Wno-unused-dummy-argument -Wno-unused-variable

INC_BASE = -I$(SRCDIR)/ensemble -I$(SRCDIR)/io \
-I$(SRCDIR)/math -I$(SRCDIR)/database -I$(TNTINCLUDE)

LIB_BASE = -L$(OBJDIR) -lprotcad  -lc -lm -lstdc++

vpath %.h $(SRCDIR)/ensemble:$(SRCDIR)/database:\
	$(SRCDIR)/ensemble:$(SRCDIR)/io:\
	$(SRCDIR)/math

vpath %.cc $(SRCDIR)/ensemble:$(SRCDIR)/database:\
	$(SRCDIR)/ensemble:$(SRCDIR)/io:\
	$(SRCDIR)/math:$(PROJDIR):
	
vpath %.f $(SRCDIR)/math

vpath %.a $(OBJDIR)

vpath %.o $(OBJDIR)

install : $(LIB_TARGETS) $(TARGETS) protcad
ifeq ($(UNAME),Linux)
	@echo export PROTCADDIR=$(PROTCADDIR) >> ~/.bashrc
	@echo export PATH=$(PATH):$(PROTCADDIR):$(PROTCADDIR)/bin >> ~/.bashrc
endif
ifeq ($(UNAME),Darwin)	
	@echo export PROTCADDIR=$(PROTCADDIR) >> ~/.bash_profile
	@echo export PATH=$(PATH):$(PROTCADDIR):$(PROTCADDIR)/bin >> ~/.bash_profile
endif

all : $(LIB_TARGETS) $(TARGETS) protcad

lib : libprotcad.a

libprotcad.a : $(LIB_CC_OBJECTS) $(LIB_F77_OBJECTS)
	cd $(OBJDIR) && ar rv libprotcad.a $?
	cd $(OBJDIR) && ranlib libprotcad.a

protAlign : libprotcad.a protAlign.cc
	cd $(OBJDIR) && $(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	cd $(OBJDIR) && strip $@ && mv $@ $(BINDIR)

protDielectric : libprotcad.a protDielectric.cc
	cd $(OBJDIR) && $(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	cd $(OBJDIR) && strip $@ && mv $@ $(BINDIR)

protDihedrals : libprotcad.a protDihedrals.cc
	cd $(OBJDIR) && $(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	cd $(OBJDIR) && strip $@ && mv $@ $(BINDIR)

protEnergy : libprotcad.a protEnergy.cc
	cd $(OBJDIR) && $(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	cd $(OBJDIR) && strip $@ && mv $@ $(BINDIR)

protEvolver : libprotcad.a protEvolver.cc
	cd $(OBJDIR) && $(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	cd $(OBJDIR) && strip $@ && mv $@ $(BINDIR)

protFolder : libprotcad.a protFolder.cc
	cd $(OBJDIR) && $(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	cd $(OBJDIR) && strip $@ && mv $@ $(BINDIR)

protInverter : libprotcad.a protInverter.cc
	cd $(OBJDIR) && $(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	cd $(OBJDIR) && strip $@ && mv $@ $(BINDIR)

protMin : libprotcad.a protMin.cc
	cd $(OBJDIR) && $(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	cd $(OBJDIR) && strip $@ && mv $@ $(BINDIR)

protMover : libprotcad.a protMover.cc
	cd $(OBJDIR) && $(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	cd $(OBJDIR) && strip $@ && mv $@ $(BINDIR)

protMutator : libprotcad.a protMutator.cc
	cd $(OBJDIR) && $(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	cd $(OBJDIR) && strip $@ && mv $@ $(BINDIR)

protOligamer : libprotcad.a protOligamer.cc
	cd $(OBJDIR) && $(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	cd $(OBJDIR) && strip $@ && mv $@ $(BINDIR)

protSequence : libprotcad.a protSequence.cc
	cd $(OBJDIR) && $(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	cd $(OBJDIR) && strip $@ && mv $@ $(BINDIR)

protShaper : libprotcad.a protShaper.cc
	cd $(OBJDIR) && $(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	cd $(OBJDIR) && strip $@ && mv $@ $(BINDIR)

protBindingEnergy : libprotcad.a protBindingEnergy.cc
	cd $(OBJDIR) && $(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	cd $(OBJDIR) && strip $@ && mv $@ $(BINDIR)

hammingdist : libprotcad.a hammingdist.cc
	cd $(OBJDIR) && $(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	cd $(OBJDIR) && strip $@ && mv $@ $(BINDIR)

$(LIB_CC_OBJECTS): %.o: %.cc %.h
	$(CXX) -c $(CFLAGS) $(INC_BASE) $< -o $@
	mv $@ $(OBJDIR)

$(LIB_F77_OBJECTS): %.o: %.f
	$(F77) -c $(FFLAGS) $< -o $@
	mv $@ $(OBJDIR)

protcad:
ifeq ($(UNAME),Linux)
	cd $(UIDIR) && qmake protcad.pro && make && strip $@ && mv $@ $(BINDIR)
endif
ifeq ($(UNAME),Darwin)
	cd $(UIDIR) && qmake protcad.pro && make && cp $(BINDIR)/protEvolver protcad.app/Contents/MacOS/
endif

clean: 
	rm -f $(OBJDIR)/*.o 
	rm -f $(OBJDIR)/*.a
	cd $(BINDIR) && rm -f $(TARGETS) && rm -f protcad
	cd $(UIDIR) && if [ -f Makefile ]; then make distclean; fi;
	
