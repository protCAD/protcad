SHELL = /bin/sh

.SUFFIXES:
.SUFFIXES: .cc .o .h .f .a

LIB_CC_OBJECTS = ran1.o ran.o point.o treeNode.o atom.o atomIterator.o residue.o chain.o residueTemplate.o allowedResidue.o secondaryStructure.o chainPosition.o residueIterator.o chainModBuffer.o molecule.o protein.o ensemble.o CMath.o generalio.o ligand.o pdbData.o pdbReader.o pdbWriter.o amberVDW.o aaBaseline.o amberElec.o rotamer.o rotamerLib.o annealer.o PDBAtomRecord.o PDBInterface.o ruler.o line.o lineSegment.o unitSphere.o solvation.o helixPropensity.o ligandTemplate.o parse.o pmf.o microEnvDB.o microEnvironment.o ramachandranMap.o

TARGETS = libprotcad y_aligner acidMutator protEvolver dFinder mergeComplex structFinder structShaper intraSoluteEnergy protOptSolvent protEvolverBinding database_phipsi protFolder sideChainRandomizer dielectricFit foldingBindingEnergy ligandBindingEnergy bindingEnergy triadFinder protMover fourEvolver protDock protMutator

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

INC_BASE = -I$(srcdir)/ensemble -I$(srcdir)/io \
-I$(srcdir)/math -I$(srcdir)/database -I$(srcdir)/algorithm \
-I$(TNTINCLUDE)

INC_GRAPHICS = -I$(srcdir)/graphics -I/usr/X11R6/include \

LIB_BASE = -L$(LIBDIR) -lprotcad  -lc -lm -lstdc++

LIB_GRAPHICS = -L/usr/X11R6/lib -lFOX \
-lglut -lpthread -lGL -lGLU -lXi -lXmu \
-lX11 -lga

vpath %.h $(srcdir)/algorithm:$(srcdir)/ensemble:$(srcdir)/database:\
	$(srcdir)/ensemble:$(srcdir)/graphics:$(srcdir)/io:\
	$(srcdir)/math

vpath %.cc $(srcdir)/algorithm:$(srcdir)/ensemble:$(srcdir)/database:\
	$(srcdir)/ensemble:$(srcdir)/graphics:$(srcdir)/io:\
	$(srcdir)/math:$(TOP)/projects:$(projdir)/McLendon/src:
	
vpath %.f $(srcdir)/math

vpath %.a $(LIBDIR)

vpath %.o $(TOP)/obj/

libprotcad.a: $(LIB_CC_OBJECTS) $(LIB_F77_OBJECTS)
		ar rv libprotcad.a $?
		ranlib libprotcad.a
		cp libprotcad.a $(TOP)/lib/

.cc.a:

	
mutantMaker : libprotcad.a mutantMaker.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/

database_phipsi : libprotcad.a database_phipsi.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/

bindingEnergy : libprotcad.a bindingEnergy.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/

fourEvolver : libprotcad.a fourEvolver.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/
	
dielectricFit : libprotcad.a dielectricFit.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/

protEvolverBinding : libprotcad.a protEvolverBinding.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/

protMutator : libprotcad.a protMutator.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/
	
foldingBindingEnergy : libprotcad.a foldingBindingEnergy.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/

triadFinder : libprotcad.a triadFinder.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/
	
ligandBindingEnergy : libprotcad.a ligandBindingEnergy.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/
	
intraSoluteEnergy : libprotcad.a intraSoluteEnergy.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/
	
protOptSolvent : libprotcad.a protOptSolvent.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/
	
protFolder : libprotcad.a protFolder.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/
	
foldingC : libprotcad.a foldingC.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/
	
protEvolverHomo : libprotcad.a protEvolverHomo.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/

structShaper : libprotcad.a structShaper.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/

protMover : libprotcad.a protMover.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/

protDock : libprotcad.a protDock.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/

protEvolver : libprotcad.a protEvolver.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/

structFinder : libprotcad.a structFinder.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/	

mergeComplex : libprotcad.a mergeComplex.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/

dFinder : libprotcad.a dFinder.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/

acidMutator : libprotcad.a acidMutator.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/

sideChainRandomizer : libprotcad.a sideChainRandomizer.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/

y_aligner : libprotcad.a y_aligner.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/

z_aligner : libprotcad.a z_aligner.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/

$(LIB_F77_OBJECTS): %.o: %.f
	$(F77) -c $(FFLAGS) $^ -o $@

$(LIB_CC_OBJECTS): %.o: %.cc %.h
	$(CXX) -c $(CFLAGS) $(INC_BASE) $< -o $@

$(LIB_GRAF_OBJECTS): %.o: %.cc %.h
	$(CXX) -c $(CFLAGS) $(INC_BASE) $(INC_GRAPHICS) $< -o $@

cleanlibs: 
	rm -f *.o 
	rm -f *.a
	rm -f $(LIBDIR)/*.a
	rm -f $(LIBDIR)/*.o
