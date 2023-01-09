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

NVCC_RESULT := $(shell which nvcc 2> NULL)
NVCC_TEST := $(notdir $(NVCC_RESULT))

export CXX = g++
export F77 = gfortran
export CU = nvcc
export MAKE = make

SHELL = /bin/sh

TARGETS = protDielectric protEvolver protDihedrals protOligamer protEnergy protMover protMutator protSequence protInverter protMin protAlign protShaper protBindingEnergy hammingdist protSampling protTest

.SUFFIXES: .cc .o .h .a .f .cu

LIB_TARGETS = lib

LIB_CC_OBJECTS = ran.o point.o treeNode.o atom.o atomIterator.o residue.o chain.o residueTemplate.o allowedResidue.o secondaryStructure.o chainPosition.o residueIterator.o chainModBuffer.o molecule.o protein.o ensemble.o CMath.o generalio.o pdbData.o pdbReader.o pdbWriter.o amberVDW.o aaBaseline.o amberElec.o rotamer.o rotamerLib.o PDBAtomRecord.o PDBInterface.o ruler.o line.o lineSegment.o unitSphere.o helixPropensity.o parse.o ramachandranMap.o

LIB_F77_OBJECTS = bestfit.o

LIB_CU_OBJECTS = energy.o

DEFS = -D__STL_USE_EXCEPTIONS

FLAG_OPTMAX = -Wall -oFast -ffast-math -ftree-vectorize -march=native -mtune=native -pipe -msse3 -Wno-deprecated -std=gnu++11

CFLAGS = $(FLAG_OPTMAX) $(DEFS)

FFLAGS = -Wall -g -Wno-tabs -Wno-unused-dummy-argument -Wno-unused-variable

CUFLAGS = -arch=sm_60

INC_BASE = -I$(SRCDIR)/ensemble -I$(SRCDIR)/io \
-I$(SRCDIR)/math -I$(SRCDIR)/database -I$(TNTINCLUDE)

ifeq ($(NVCC_TEST),nvcc)
	LIB_BASE = -L$(OBJDIR) -lprotcad -lc -lm -lstdc++ -lcuda -lcudart
	LIBS = $(LIB_CC_OBJECTS) $(LIB_F77_OBJECTS) $(LIB_CU_OBJECTS)
else
	LIB_BASE = -L$(OBJDIR) -lprotcad -lc -lm -lstdc++
	LIBS = $(LIB_CC_OBJECTS) $(LIB_F77_OBJECTS)
endif

vpath %.h $(SRCDIR)/ensemble:$(SRCDIR)/database:\
	$(SRCDIR)/ensemble:$(SRCDIR)/io:\
	$(SRCDIR)/math

vpath %.cc $(SRCDIR)/ensemble:$(SRCDIR)/database:\
	$(SRCDIR)/ensemble:$(SRCDIR)/io:\
	$(SRCDIR)/math:$(PROJDIR):

vpath %.f $(SRCDIR)/math

vpath %.cu $(SRCDIR)/ensemble

vpath %.a $(OBJDIR)

vpath %.o $(OBJDIR)

install : $(LIB_TARGETS) $(TARGETS)
ifeq ($(UNAME),Linux)
	@echo export PROTCADDIR=$(PROTCADDIR) >> ~/.bashrc
	@echo export PATH=$(PATH):$(PROTCADDIR):$(PROTCADDIR)/bin >> ~/.bashrc
endif
ifeq ($(UNAME),Darwin)
	@echo export PROTCADDIR=$(PROTCADDIR) >> ~/.zshrc
	@echo export PATH=$(PATH):$(PROTCADDIR):$(PROTCADDIR)/bin >> ~/.zshrc
endif

all : $(LIB_TARGETS) $(TARGETS)

lib : libprotcad.a

libprotcad.a : $(LIBS)
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

protSampling : libprotcad.a protSampling.cc
	cd $(OBJDIR) && $(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	cd $(OBJDIR) && strip $@ && mv $@ $(BINDIR)

protTest : libprotcad.a protTest.cc
	cd $(OBJDIR) && $(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	cd $(OBJDIR) && strip $@ && mv $@ $(BINDIR)

$(LIB_CC_OBJECTS): %.o: %.cc %.h
	$(CXX) -c $(CFLAGS) $(INC_BASE) $< -o $@
	mv $@ $(OBJDIR)

$(LIB_F77_OBJECTS): %.o: %.f
	$(F77) -c $(FFLAGS) $< -o $@
	mv $@ $(OBJDIR)

$(LIB_CU_OBJECTS): %.o: %.cu
	$(CU) -c $(CUFLAGS) $< -o $@
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
