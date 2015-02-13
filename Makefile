export TOP=$(PWD)
export SRCDIR=$(TOP)/src
export TNTINCLUDE=$(TOP)/src
export BINDIR=$(TOP)/bin
export OBJDIR=$(TOP)/obj
export PROJDIR=$(TOP)/projects

export CXX = g++
export F77 = g77
export MAKE = make

SHELL = /bin/sh

TARGETS = mutantMaker acidMutator protEvolver mergeComplex structFinder structShaper intraSoluteEnergy protOptSolvent protEvolverBinding database_phipsi protFolder sideChainRandomizer dielectricFit foldingBindingEnergy ligandBindingEnergy bindingEnergy triadFinder protMover z_aligner y_aligner fourEvolver protDock protMutator

.SUFFIXES:
.SUFFIXES: .cc .o .h .f .a

LIB_TARGETS = lib

LIB_CC_OBJECTS = ran1.o ran.o point.o treeNode.o atom.o atomIterator.o residue.o chain.o residueTemplate.o allowedResidue.o secondaryStructure.o chainPosition.o residueIterator.o chainModBuffer.o molecule.o protein.o ensemble.o CMath.o generalio.o ligand.o pdbData.o pdbReader.o pdbWriter.o amberVDW.o aaBaseline.o amberElec.o rotamer.o rotamerLib.o annealer.o PDBAtomRecord.o PDBInterface.o ruler.o line.o lineSegment.o unitSphere.o solvation.o helixPropensity.o ligandTemplate.o parse.o pmf.o microEnvDB.o microEnvironment.o ramachandranMap.o

#DEFS = -DHAVE_OPENGL=1 -D_ALLOWED_RESIDUE_DEBUG
#DEFS = -DHAVE_OPENGL=1 -D__STL_USE_EXCEPTIONS
#DEFS = -DHAVE_OPENGL=1 -D__STL_USE_EXCEPTIONS -DMICROENV_DEBUG_ATOM_TYPES -DATOM_TYPE_DEBUG \
-DMICROENVDB_DEBUG

DEFS = -DHAVE_OPENGL=1 -D__STL_USE_EXCEPTIONS -DMICROENV_DEBUG_ATOM_TYPES -DATOM_TYPE_DEBUG

FLAG_OPT = -Wall -O -g -felide-constructors -Wno-deprecated
FLAG_OPT2 = -Wall -O2 -g -Wno-deprecated
FLAG_OPT3 = -Wall -O3  -g -felide-constructors -Wno-deprecated
FLAG_PROF = -Wall -O3 -felide-constructors -pg -Wno-deprecated
FLAG_DEBUG = -Wall -g2 -felide-constructors -Wno-deprecated
FLAG_DEBUG2 = -Wall -g2 -ansi -pedantic -Wno-deprecated

CFLAGS = $(FLAG_DEBUG2) $(DEFS)
FFLAGS = -Wall -g 

INC_BASE = -I$(SRCDIR)/ensemble -I$(SRCDIR)/io \
-I$(SRCDIR)/math -I$(SRCDIR)/database -I$(SRCDIR)/algorithm \
-I$(TNTINCLUDE) -I$(OBJDIR)

LIB_BASE = -L$(OBJDIR) -lprotcad  -lc -lm -lstdc++

vpath %.h $(SRCDIR)/algorithm:$(SRCDIR)/ensemble:$(SRCDIR)/database:\
	$(SRCDIR)/ensemble:$(SRCDIR)/io:\
	$(SRCDIR)/math

vpath %.cc $(SRCDIR)/algorithm:$(SRCDIR)/ensemble:$(SRCDIR)/database:\
	$(SRCDIR)/ensemble:$(SRCDIR)/io:\
	$(SRCDIR)/math:$(PROJDIR):
	
vpath %.f $(SRCDIR)/math

vpath %.a $(OBJDIR)

vpath %.o $(OBJDIR)

all : $(LIB_TARGETS) $(TARGETS)

lib : $(LIB_CC_OBJECTS) $(LIB_F77_OBJECTS)
		cd $(OBJDIR) && ar rv libprotcad.a $?
		cd $(OBJDIR) && ranlib libprotcad.a
		
mutantMaker : libprotcad.a mutantMaker.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ $(BINDIR)

database_phipsi : libprotcad.a database_phipsi.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ $(BINDIR)

bindingEnergy : libprotcad.a bindingEnergy.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ $(BINDIR)

fourEvolver : libprotcad.a fourEvolver.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ $(BINDIR)
	
dielectricFit : libprotcad.a dielectricFit.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ $(BINDIR)

protEvolverBinding : libprotcad.a protEvolverBinding.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ $(BINDIR)

protMutator : libprotcad.a protMutator.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ $(BINDIR)
	
foldingBindingEnergy : libprotcad.a foldingBindingEnergy.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ $(BINDIR)

triadFinder : libprotcad.a triadFinder.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ $(BINDIR)
	
ligandBindingEnergy : libprotcad.a ligandBindingEnergy.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ $(BINDIR)
	
intraSoluteEnergy : libprotcad.a intraSoluteEnergy.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ $(BINDIR)
	
protOptSolvent : libprotcad.a protOptSolvent.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ $(BINDIR)
	
protFolder : libprotcad.a protFolder.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ $(BINDIR)
	
foldingC : libprotcad.a foldingC.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ $(BINDIR)
	
protEvolverHomo : libprotcad.a protEvolverHomo.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ $(BINDIR)

structShaper : libprotcad.a structShaper.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ $(BINDIR)

protMover : libprotcad.a protMover.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ $(BINDIR)

protDock : libprotcad.a protDock.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ $(BINDIR)

protEvolver : libprotcad.a protEvolver.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ $(BINDIR)

structFinder : libprotcad.a structFinder.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ $(BINDIR)	

mergeComplex : libprotcad.a mergeComplex.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ $(BINDIR)

dFinder : libprotcad.a dFinder.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ $(BINDIR)

acidMutator : libprotcad.a acidMutator.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ $(BINDIR)

sideChainRandomizer : libprotcad.a sideChainRandomizer.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ $(BINDIR)

y_aligner : libprotcad.a y_aligner.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ $(BINDIR)

z_aligner : libprotcad.a z_aligner.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ $(BINDIR)

$(LIB_F77_OBJECTS): %.o: %.f
	$(F77) -c $(FFLAGS) $^ -o $@
	mv $@ $(OBJDIR)

$(LIB_CC_OBJECTS): %.o: %.cc %.h
	$(CXX) -c $(CFLAGS) $(INC_BASE) $< -o $@
	mv $@ $(OBJDIR)

clean: 
	rm -f $(OBJDIR)/*.o 
	rm -f $(OBJDIR)/*.a
	rm -f $(OBJDIR)/Makefile
	cd $(BINDIR) && rm -f $(TARGETS)
