export TOP=$(PWD)
export srcdir=$(TOP)/src
export TNTINCLUDE=$(TOP)/src
export BINDIR=$(TOP)/bin
export LIBDIR=$(TOP)/lib

export CXX = g++
export F77 = g77
export MAKE = make

LIB_TARGETS = lib
EXEC_TARGETS = mutantMaker acidMutator protEvolver mergeComplex structFinder structShaper intraSoluteEnergy protOptSolvent protEvolverBinding database_phipsi protFolder sideChainRandomizer dielectricFit foldingBindingEnergy ligandBindingEnergy bindingEnergy triadFinder protMover z_aligner y_aligner fourEvolver protDock protMutator

all : $(LIB_TARGETS) $(EXEC_TARGETS)

lib : FORCE
	cd obj && $(MAKE) libprotcad.a

fourEvolver : FORCE
	cd obj && $(MAKE) fourEvolver

sideChainRandomizer : FORCE
	cd obj && $(MAKE) sideChainRandomizer

protDock : FORCE
	cd obj && $(MAKE) protDock

protMutator : FORCE
	cd obj && $(MAKE) protMutator
acidMutator : FORCE
	cd obj && $(MAKE) acidMutator
	
foldingBindingEnergy : FORCE
	cd obj && $(MAKE) foldingBindingEnergy
	
database_phipsi : FORCE
	cd obj && $(MAKE) database_phipsi
	
intraSoluteEnergy : FORCE
	cd obj && $(MAKE) intraSoluteEnergy

protMover : FORCE
	cd obj && $(MAKE) protMover
	
bindingEnergy : FORCE
	cd obj && $(MAKE) bindingEnergy
	
protFolder : FORCE
	cd obj && $(MAKE) protFolder

triadFinder : FORCE
	cd obj && $(MAKE) triadFinder

dielectricFit : FORCE
	cd obj && $(MAKE) dielectricFit
	
mergeComplex : FORCE
	cd obj && $(MAKE) mergeComplex
	
ligandBindingEnergy : FORCE
	cd obj && $(MAKE) ligandBindingEnergy

protEvolverBinding : FORCE
	cd obj && $(MAKE) protEvolverBinding

structFinder : FORCE
	cd obj && $(MAKE) structFinder

protOptSolvent : FORCE
	cd obj && $(MAKE) protOptSolvent

protEvolver : FORCE
	cd obj && $(MAKE) protEvolver

structShaper : FORCE
	cd obj && $(MAKE) structShaper

mutateGFPSurface : FORCE
	cd obj && $(MAKE) mutateGFPSurface

mutantMaker : FORCE
	cd obj && $(MAKE) mutantMaker

dFinder : FORCE
	cd obj && $(MAKE) dFinder

y_aligner : FORCE
	cd obj && $(MAKE) y_aligner

z_aligner : FORCE
	cd obj && $(MAKE) z_aligner

clean: FORCE
	cd bin && rm $(EXEC_TARGETS)
	cd obj && $(MAKE) cleanlibs

FORCE:
