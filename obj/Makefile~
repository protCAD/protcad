SHELL = /bin/sh

.SUFFIXES:
.SUFFIXES: .cc .o .h .f .a

LIB_CC_OBJECTS = ran1.o ran.o point.o treeNode.o atom.o atomIterator.o residue.o chain.o residueTemplate.o allowedResidue.o secondaryStructure.o chainPosition.o residueIterator.o chainModBuffer.o molecule.o protein.o ensemble.o CMath.o generalio.o ligand.o pdbData.o pdbReader.o pdbWriter.o amberVDW.o aaBaseline.o amberElec.o rotamer.o rotamerLib.o annealer.o PDBAtomRecord.o PDBInterface.o ruler.o line.o lineSegment.o unitSphere.o solvation.o helixPropensity.o ligandTemplate.o parse.o

TARGETS = libprotcad y_aligner acidMutator protEvolver dFinder protDock mergeComplex structFinder structShaper dielectric intraSoluteEnergy chainBindingEnergy protOptSolvent aggreSim amberAnalyzer protEvolverBinding database_phipsi protFolder sideChainRandomizer simpleFolder dielectricFit foldingBindingEnergy ligandBindingEnergy bindingEnergy triadFinder

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
-I$(GALIBINCLUDE) -I$(TNTINCLUDE)

INC_GRAPHICS = -I$(srcdir)/graphics -I/usr/X11R6/include \
iI$(GALIBINCLUDE) -I$(FOXINCLUDE)

LIB_BASE = -L$(LIBDIR) -lprotcad  -lc -lm -lstdc++

LIB_GRAPHICS = -L/usr/X11R6/lib -lFOX \
-lglut -lpthread -lGL -lGLU -lXi -lXmu \
-lX11 -lga

vpath %.h $(srcdir)/algorithm:$(srcdir)/ensemble:$(srcdir)/database:\
	$(srcdir)/ensemble:$(srcdir)/graphics:$(srcdir)/io:\
	$(srcdir)/math

vpath %.cc $(srcdir)/algorithm:$(srcdir)/ensemble:$(srcdir)/database:\
	$(srcdir)/ensemble:$(srcdir)/graphics:$(srcdir)/io:\
	$(srcdir)/math:$(srcdir)/test:$(projdir)/McLendon/src:\
	$(projdir)/glycophorin/src:$(projdir)/GXXG/src:\
	$(projdir)/Synaptobrevin/src:$(projdir)/kchannel/src: \
	$(projdir)/liang/src:$(projdir)/biphenyl/src:$(projdir)/jeff/src:\
	$(projdir)/plaw/src:$(projdir)/alessandrotest

vpath %.f $(srcdir)/math

vpath %.a $(LIBDIR)

vpath %.o $(TOP)/obj/$(ARCH)

libprotcad.a: $(LIB_CC_OBJECTS) $(LIB_F77_OBJECTS)
		cd $(ARCH) && ar rv libprotcad.a $?
		ranlib $(ARCH)/libprotcad.a
		cp $(ARCH)/libprotcad.a $(TOP)/lib/$(ARCH)

.cc.a:

	
mutantMaker : libprotcad.a mutantMaker.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/$(ARCH)

database_phipsi : libprotcad.a database_phipsi.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/$(ARCH)

bindingEnergy : libprotcad.a bindingEnergy.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/$(ARCH)
	
dielectricFit : libprotcad.a dielectricFit.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/$(ARCH)

protEvolverBinding : libprotcad.a protEvolverBinding.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/$(ARCH)
	
foldingBindingEnergy : libprotcad.a foldingBindingEnergy.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/$(ARCH)

triadFinder : libprotcad.a triadFinder.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/$(ARCH)
	
simpleFolder : libprotcad.a simpleFolder.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/$(ARCH)
	
ligandBindingEnergy : libprotcad.a ligandBindingEnergy.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/$(ARCH)
	
intraSoluteEnergy : libprotcad.a intraSoluteEnergy.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/$(ARCH)
	
amberAnalyzer : libprotcad.a amberAnalyzer.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/$(ARCH)
	
protOptSolvent : libprotcad.a protOptSolvent.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/$(ARCH)
	
protFolder : libprotcad.a protFolder.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/$(ARCH)
	
foldingC : libprotcad.a foldingC.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/$(ARCH)
	
chainBindingEnergy : libprotcad.a chainBindingEnergy.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/$(ARCH)

protEvolverHomo : libprotcad.a protEvolverHomo.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/$(ARCH)
	
aggreSim : libprotcad.a aggreSim.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/$(ARCH)

structShaper : libprotcad.a structShaper.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/$(ARCH)

protMover : libprotcad.a protMover.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/$(ARCH)

dielectric : libprotcad.a dielectric.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/$(ARCH)

protEvolver : libprotcad.a protEvolver.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/$(ARCH)

structFinder : libprotcad.a structFinder.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/$(ARCH)

protDock : libprotcad.a protDock.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/$(ARCH)	

mergeComplex : libprotcad.a mergeComplex.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/$(ARCH)

dFinder : libprotcad.a dFinder.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/$(ARCH)

acidMutator : libprotcad.a acidMutator.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/$(ARCH)

sideChainRandomizer : libprotcad.a sideChainRandomizer.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/$(ARCH)

y_aligner : libprotcad.a y_aligner.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/$(ARCH)

z_aligner : libprotcad.a z_aligner.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/$(ARCH)

$(LIB_F77_OBJECTS): %.o: %.f
	$(F77) -c $(FFLAGS) $^ -o $(ARCH)/$@

$(LIB_CC_OBJECTS): %.o: %.cc %.h
	$(CXX) -c $(CFLAGS) $(INC_BASE) $< -o $(ARCH)/$@

$(LIB_GRAF_OBJECTS): %.o: %.cc %.h
	$(CXX) -c $(CFLAGS) $(INC_BASE) $(INC_GRAPHICS) $< -o $(ARCH)/$@

zzyxyz: 
	rm -f $(ARCH)/*.o 
	rm -f $(ARCH)/*.a
	rm -f $(LIBDIR)/*.a
