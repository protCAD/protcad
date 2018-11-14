export PROTCADDIR=$(PWD)
export SRCDIR=$(PROTCADDIR)/src
export TNTINCLUDE=$(SRCDIR)
export BINDIR=$(PROTCADDIR)/bin
export OBJDIR=$(PROTCADDIR)/obj
export PROJDIR=$(PROTCADDIR)/projects
export $(PATH)=$(PATH):$(PROTCADDIR):$(BINDIR)

export CXX = g++
export MAKE = make

SHELL = /bin/sh

TARGETS = protDielectric protEvolver protMerge protDihedrals protShaper protEnergy protOpt protEvolver protEvolverH protDB protFolder protRandomizer protBindingEnergy triadFinder protMover z_aligner y_aligner protDock protMutator protPointMutator protSequence protInverter protSorter protRotamer protSlipPlane alphaCarbonDihedrals protClashes

.SUFFIXES: .cc .o .h .a

LIB_TARGETS = lib

LIB_CC_OBJECTS = ran1.o ran.o point.o treeNode.o atom.o atomIterator.o residue.o chain.o residueTemplate.o allowedResidue.o secondaryStructure.o chainPosition.o residueIterator.o chainModBuffer.o molecule.o protein.o ensemble.o CMath.o generalio.o pdbData.o pdbReader.o pdbWriter.o amberVDW.o aaBaseline.o amberElec.o rotamer.o rotamerLib.o PDBAtomRecord.o PDBInterface.o ruler.o line.o lineSegment.o unitSphere.o helixPropensity.o parse.o ramachandranMap.o

DEFS = -DHAVE_OPENGL=1 -D__STL_USE_EXCEPTIONS

FLAG_OPTMAX = -Wall -oFast -ffast-math -ftree-vectorize -march=native -mtune=native -pipe -msse3 -Wno-deprecated -std=gnu++11

CFLAGS = $(FLAG_OPTMAX) $(DEFS)

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

install : $(LIB_TARGETS)
	@echo export PROTCADDIR=$(PROTCADDIR) >> ~/.bashrc
	@echo export PATH=$(PATH):$(PROTCADDIR):$(PROTCADDIR)/bin >> ~/.bashrc

all : $(LIB_TARGETS) $(TARGETS)

lib : libprotcad.a

libprotcad.a : $(LIB_CC_OBJECTS)
	cd $(OBJDIR) && ar rv libprotcad.a $?
	cd $(OBJDIR) && ranlib libprotcad.a

protMutator : libprotcad.a protMutator.cc
	cd $(OBJDIR) && $(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	cd $(OBJDIR) && strip $@ && mv $@ $(BINDIR)

protRotamer : libprotcad.a protRotamer.cc
	cd $(OBJDIR) && $(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	cd $(OBJDIR) && strip $@ && mv $@ $(BINDIR)

protPointMutator : libprotcad.a protPointMutator.cc
	cd $(OBJDIR) && $(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	cd $(OBJDIR) && strip $@ && mv $@ $(BINDIR)

protInverter : libprotcad.a protInverter.cc
	cd $(OBJDIR) && $(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	cd $(OBJDIR) && strip $@ && mv $@ $(BINDIR)

protSorter : libprotcad.a protSorter.cc
	cd $(OBJDIR) && $(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	cd $(OBJDIR) && strip $@ && mv $@ $(BINDIR)

protSequence : libprotcad.a protSequence.cc
	cd $(OBJDIR) && $(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	cd $(OBJDIR) && strip $@ && mv $@ $(BINDIR)

protDB : libprotcad.a protDB.cc
	cd $(OBJDIR) && $(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	cd $(OBJDIR) && strip $@ && mv $@ $(BINDIR)

protBindingEnergy : libprotcad.a protBindingEnergy.cc
	cd $(OBJDIR) && $(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	cd $(OBJDIR) && strip $@ && mv $@ $(BINDIR)

protSlipPlane : libprotcad.a protSlipPlane.cc
	cd $(OBJDIR) && $(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	cd $(OBJDIR) && strip $@ && mv $@ $(BINDIR)
	
protDielectric : libprotcad.a protDielectric.cc
	cd $(OBJDIR) && $(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	cd $(OBJDIR) && strip $@ && mv $@ $(BINDIR)

triadFinder : libprotcad.a triadFinder.cc
	cd $(OBJDIR) && $(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	cd $(OBJDIR) && strip $@ && mv $@ $(BINDIR)
	
protEnergy : libprotcad.a protEnergy.cc
	cd $(OBJDIR) && $(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	cd $(OBJDIR) && strip $@ && mv $@ $(BINDIR)
	
protOpt : libprotcad.a protOpt.cc
	cd $(OBJDIR) && $(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	cd $(OBJDIR) && strip $@ && mv $@ $(BINDIR)
	
protFolder : libprotcad.a protFolder.cc
	cd $(OBJDIR) && $(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	cd $(OBJDIR) && strip $@ && mv $@ $(BINDIR)

protShaper : libprotcad.a protShaper.cc
	cd $(OBJDIR) && $(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	cd $(OBJDIR) && strip $@ && mv $@ $(BINDIR)

protMover : libprotcad.a protMover.cc
	cd $(OBJDIR) && $(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	cd $(OBJDIR) && strip $@ && mv $@ $(BINDIR)

protDock : libprotcad.a protDock.cc
	cd $(OBJDIR) && $(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	cd $(OBJDIR) && strip $@ && mv $@ $(BINDIR)

protClashes : libprotcad.a protClashes.cc
	cd $(OBJDIR) && $(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	cd $(OBJDIR) && strip $@ && mv $@ $(BINDIR)

protEvolver : libprotcad.a protEvolver.cc
	cd $(OBJDIR) && $(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	cd $(OBJDIR) && strip $@ && mv $@ $(BINDIR)

protEvolverH : libprotcad.a protEvolverH.cc
	cd $(OBJDIR) && $(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	cd $(OBJDIR) && strip $@ && mv $@ $(BINDIR)

protDihedrals : libprotcad.a protDihedrals.cc
	cd $(OBJDIR) && $(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	cd $(OBJDIR) && strip $@ && mv $@ $(BINDIR)

protMerge : libprotcad.a protMerge.cc
	cd $(OBJDIR) && $(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	cd $(OBJDIR) && strip $@ && mv $@ $(BINDIR)

protRandomizer : libprotcad.a protRandomizer.cc
	cd $(OBJDIR) && $(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	cd $(OBJDIR) && strip $@ && mv $@ $(BINDIR)

alphaCarbonDihedrals : libprotcad.a alphaCarbonDihedrals.cc
	cd $(OBJDIR) && $(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	cd $(OBJDIR) && strip $@ && mv $@ $(BINDIR)

y_aligner : libprotcad.a y_aligner.cc
	cd $(OBJDIR) && $(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	cd $(OBJDIR) && strip $@ && mv $@ $(BINDIR)

z_aligner : libprotcad.a z_aligner.cc
	cd $(OBJDIR) && $(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	cd $(OBJDIR) && strip $@ && mv $@ $(BINDIR)

$(LIB_CC_OBJECTS): %.o: %.cc %.h
	$(CXX) -c $(CFLAGS) $(INC_BASE) $< -o $@
	mv $@ $(OBJDIR)

clean: 
	rm -f $(OBJDIR)/*.o 
	rm -f $(OBJDIR)/*.a
	cd $(BINDIR) && rm -f $(TARGETS)
