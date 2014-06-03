SHELL = /bin/sh

.SUFFIXES:
.SUFFIXES: .cc .o .h .f .a

LIB_CC_OBJECTS = ran1.o ran.o point.o treeNode.o atom.o \
	atomIterator.o residue.o chain.o residueTemplate.o \
	allowedResidue.o secondaryStructure.o chainPosition.o \
	ramachandranMap.o residueIterator.o chainModBuffer.o \
	molecule.o protein.o ensemble.o CMath.o generalio.o \
	ligand.o pdbData.o pdbReader.o pdbWriter.o amberVDW.o pmf.o \
	microEnvDB.o microEnvironment.o aaBaseline.o amberElec.o rotamer.o \
	rotamerLib.o annealer.o  \
	deadEndEliminator.o PDBAtomRecord.o PDBInterface.o ruler.o line.o \
	lineSegment.o unitSphere.o solvation.o helixPropensity.o ligandTemplate.o parse.o

#LIB_F77_OBJECTS = bestfit.o

#LIB_GRAF_OBJECTS = position.o lineObject.o myFXGLLine.o glviewWindow.o \
	myFXGLSphere.o FXGLProtViewer.o FXRepackDialog.o \
	gaDriver.o FXEditDialog.o
	
TARGETS = libprotcad hemePacker hetatmTester hetatmProtTest hemeBundlePacker \
	biphenPacker nativeDesigner glycoHBondEnergy2 primalSoup chiralTransfer\
	primalSheet chiralTest primalHelixLoopHelix chiralInteraction z_aligner\
	dimerBuilder trimerBuilder originator repack zHelixVectorAligner y_aligner \
	nitrileProjection downHillSimplex mutationDocker thrRotamer gcn4mutator \
	serAntiDimer serDimer homoOligomer helphipsiscan turnBuilder DFscpacker\
	helixangle deberdesign helixEval fiber ramachandran patchCalc curveHelix \
	pentamerBuilder backboneMerge coords evenMC mapPhiPsi mapSequence m2Modeler\
	calcHelixAngle m2Fit changeLoop skeleton mutateGFPSurface coilBuilder Energy\
	lowestRotCUDA rotOpt rotOpthelix helixMeasure electroMeasure indigestible \
	repulseEnergy attractEnergy digestible gastricAnalyzer counter positionAnalyzer \
	denovoBuilder denovoBuilderSym test savingStateTest indigestibleSym digestibleSym \
	pkaAve digestIndex resdigestIndex patternFinder chiralFlipper bbOpt condOpt fiveElements \
	protOpt ionicEnergy acidMutator folding protEvolver dFinder \
	protDock mergeComplex oligOpt structFinder structShaper \
	pointMutator_ins impSolvent dielectric burialSurface carbonDensity intraSoluteEnergy \
	chainBindingEnergy protOptSolvent foldingC aggreSim amberAnalyzer protEvolverBinding \
	peptideBondSASA traj_deformation database_phipsi protFolder sideChainRandomizer \
	simpleFolder protRandomizer dielectricFit foldingBindingEnergy ligandBindingEnergy \
	solvationPDB bindingEnergy triadFinder

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

#INC_BASE = -I$(SVMTINCLUDE) -I$(srcdir)/ensemble -I$(srcdir)/io \
#-I$(srcdir)/math -I$(srcdir)/database -I$(srcdir)/algorithm \
#-I$(GALIBINCLUDE) -I$(TNTINCLUDE)

INC_BASE = -I$(srcdir)/ensemble -I$(srcdir)/io \
-I$(srcdir)/math -I$(srcdir)/database -I$(srcdir)/algorithm \
-I$(GALIBINCLUDE) -I$(TNTINCLUDE)

INC_GRAPHICS = -I$(srcdir)/graphics -I/usr/X11R6/include \
iI$(GALIBINCLUDE) -I$(FOXINCLUDE)

#LIB_BASE = -L$(LIBDIR) \
#-lprotcad -ldoubleComplexTensor -ldoubleTensor -lboolTensor \
#-ldoubleComplexSquare -ldoubleSquare \
#-ldoubleComplexMatrix -ldoubleMatrix -lboolMatrix \
#-ldoubleComplexVector -ldoubleVector -lboolVector -loffsetVector \
#-ldoubleComplex -lsvmt_error -lm -lga -lstdc++

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

libprotcadgraphics.a : $(LIB_GRAF_OBJECTS)
		cd $(ARCH) && ar rv libprotcadgraphics.a $?
		ranlib $(ARCH)/libprotcadgraphics.a
		cp $(ARCH)/libprotcadgraphics.a $(TOP)/lib/$(ARCH)

.cc.a:



mutateGFPSurface : libprotcad.a mutateGFPSurface.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/$(ARCH)

solvationPDB : libprotcad.a solvationPDB.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/$(ARCH)

coilBuilder : libprotcad.a coilBuilder.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/$(ARCH)
	
mutantMaker : libprotcad.a mutantMaker.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/$(ARCH)

functions : libprotcad.a functions.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/$(ARCH)
	
traj_deformation : libprotcad.a traj_deformation.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/$(ARCH)
	
database_phipsi : libprotcad.a database_phipsi.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/$(ARCH)

burialSurface : libprotcad.a burialSurface.cc
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
	
peptideBondSASA : libprotcad.a peptideBondSASA.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/$(ARCH)
	
protRandomizer : libprotcad.a protRandomizer.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/$(ARCH)

pointMutator : libprotcad.a pointMutator.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/$(ARCH)
	
carbonDensity : libprotcad.a carbonDensity.cc
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

pointMutator_ins : libprotcad.a pointMutator_ins.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/$(ARCH)

structShaper : libprotcad.a structShaper.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/$(ARCH)

protMover : libprotcad.a protMover.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/$(ARCH)

impSolvent : libprotcad.a impSolvent.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/$(ARCH)

dielectric : libprotcad.a dielectric.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/$(ARCH)

fiveElements : libprotcad.a fiveElements.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/$(ARCH)

protEvolver : libprotcad.a protEvolver.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/$(ARCH)

oligOpt : libprotcad.a oligOpt.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/$(ARCH)

structFinder : libprotcad.a structFinder.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/$(ARCH)

smartMutator : libprotcad.a smartMutator.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/$(ARCH)

folding : libprotcad.a folding.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/$(ARCH)

protDock : libprotcad.a protDock.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/$(ARCH)	

mergeComplex : libprotcad.a mergeComplex.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/$(ARCH)

rotOpt : libprotcad.a rotOpt.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/$(ARCH)

protOpt : libprotcad.a protOpt.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/$(ARCH)

protOpt0 : libprotcad.a protOpt0.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/$(ARCH)

dFinder : libprotcad.a dFinder.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/$(ARCH)

acidMutator : libprotcad.a acidMutator.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/$(ARCH)

wtOptimizer : libprotcad.a wtOptimizer.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/$(ARCH)

atomMap : libprotcad.a atomMap.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/$(ARCH)
	
sideChainRandomizer : libprotcad.a sideChainRandomizer.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/$(ARCH)

ionicEnergy : libprotcad.a ionicEnergy.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/$(ARCH)

condOpt : libprotcad.a condOpt.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/$(ARCH)

bbOpt : libprotcad.a bbOpt.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/$(ARCH)

pkaAve : libprotcad.a pkaAve.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/$(ARCH)

attractEnergy : libprotcad.a attractEnergy.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/$(ARCH)

chiralFlipper : libprotcad.a chiralFlipper.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/$(ARCH)

savingStateTest : libprotcad.a savingStateTest.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/$(ARCH)

repulseEnergy : libprotcad.a repulseEnergy.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/$(ARCH)

digestIndex : libprotcad.a digestIndex.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/$(ARCH)

resdigestIndex : libprotcad.a resdigestIndex.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/$(ARCH)

denovoBuilder : libprotcad.a denovoBuilder.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/$(ARCH)

denovoBuilderSym : libprotcad.a denovoBuilderSym.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/$(ARCH)

patternFinder : libprotcad.a patternFinder.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/$(ARCH)

indigestible : libprotcad.a indigestible.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/$(ARCH)

indigestibleEvo : libprotcad.a indigestibleEvo.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/$(ARCH)

indigestibleSym : libprotcad.a indigestibleSym.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/$(ARCH)

positionAnalyzer : libprotcad.a positionAnalyzer.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/$(ARCH)

gastricAnalyzer : libprotcad.a gastricAnalyzer.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/$(ARCH)

digestible : libprotcad.a digestible.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/$(ARCH)

digestibleSym : libprotcad.a digestibleSym.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/$(ARCH)

test : libprotcad.a test.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/$(ARCH)
	
rotOpthelix : libprotcad.a rotOpthelix.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/$(ARCH)

counter : libprotcad.a counter.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/$(ARCH)

helixMeasure : libprotcad.a helixMeasure.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/$(ARCH)
	
skeleton : libprotcad.a skeleton.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/$(ARCH)

changeLoop : libprotcad.a changeLoop.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/$(ARCH)

calcHelixAngle : libprotcad.a calcHelixAngle.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/$(ARCH)

m2Fit : libprotcad.a m2Fit.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/$(ARCH)

m2Modeler : libprotcad.a m2Modeler.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/$(ARCH)

Energy : libprotcad.a Energy.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/$(ARCH)

lowestRotCUDA : libprotcad.a lowestRotCUDA.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/$(ARCH)

mapSequence : libprotcad.a mapSequence.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/$(ARCH)

mapPhiPsi : libprotcad.a mapPhiPsi.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/$(ARCH)

evenMC : libprotcad.a evenMC.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/$(ARCH)

coords: libprotcad.a coords.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/$(ARCH)

backboneMerge : libprotcad.a backboneMerge.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/$(ARCH)

pentamerBuilder : libprotcad.a pentamerBuilder.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/$(ARCH)

curveHelix : libprotcad.a curveHelix.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/$(ARCH)

patchCalc : libprotcad.a patchCalc.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/$(ARCH)

ramachandran : libprotcad.a ramachandran.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/$(ARCH)

fiber : libprotcad.a fiber.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/$(ARCH)

helixEval : libprotcad.a helixEval.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/$(ARCH)

deberdesign : libprotcad.a deberdesign.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/$(ARCH)

helixangle : libprotcad.a helixangle.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/$(ARCH)

deber : libprotcad.a deber.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/$(ARCH)

ms1Dimer : libprotcad.a ms1Dimer.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/$(ARCH)

serAntiDimer : libprotcad.a serAntiDimer.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/$(ARCH)

serDimer : libprotcad.a serDimer.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/$(ARCH)

DFscpacker : libprotcad.a DFscpacker.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/$(ARCH)

turnBuilder : libprotcad.a turnBuilder.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/$(ARCH)

helphipsiscan : libprotcad.a helphipsiscan.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/$(ARCH)

gcn4mutator : libprotcad.a gcn4mutator.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/$(ARCH)

homoOligomer : libprotcad.a homoOligomer.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/$(ARCH)


thrRotamer : libprotcad.a thrRotamer.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/$(ARCH)

mutationDocker : libprotcad.a mutationDocker.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/$(ARCH)

downHillSimplex : libprotcad.a downHillSimplex.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/$(ARCH)

zHelixVectorAligner : libprotcad.a zHelixVectorAligner.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/$(ARCH)

repack : libprotcad.a repack.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/$(ARCH)

originator : libprotcad.a originator.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/$(ARCH)

trimerBuilder : libprotcad.a trimerBuilder.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/$(ARCH)

dimerBuilder : libprotcad.a dimerBuilder.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/$(ARCH)

nitrileProjection : libprotcad.a nitrileProjection.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/$(ARCH)

y_aligner : libprotcad.a y_aligner.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/$(ARCH)

z_aligner : libprotcad.a z_aligner.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/$(ARCH)

chiralInteraction : libprotcad.a chiralInteraction.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/$(ARCH)

primalHelixLoopHelix : libprotcad.a primalHelixLoopHelix.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/$(ARCH)

glycoDockerGridRef : libprotcad.a glycoDockerGridRef.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/$(ARCH)

chiralTest : libprotcad.a chiralTest.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/$(ARCH)

primalSheet : libprotcad.a primalSheet.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/$(ARCH)

chiralTransfer : libprotcad.a chiralTransfer.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/$(ARCH)

primalSoup : libprotcad.a primalSoup.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/$(ARCH)

glycoHBondEnergy2 : libprotcad.a glycoHBondEnergy2.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/$(ARCH)

nativeDesigner : libprotcad.a nativeDesigner.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	mv $@ ../bin/$(ARCH)

trypsin : libprotcad.a trypsin.cc
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
