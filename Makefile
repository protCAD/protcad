export TOP=$(PWD)
export srcdir=$(TOP)/src
export GALIBINCLUDE=$(TOP)/galib245
export TNTINCLUDE=$(TOP)

export ARCH=$(MACHTYPE)

export LIBDIR=$(TOP)/lib/$(ARCH)

export CXX = g++
export F77 = g77
export MAKE = make

LIB_TARGETS = lib
EXEC_TARGETS = mutantMaker acidMutator protEvolver protDock mergeComplex structFinder structShaper dielectric intraSoluteEnergy chainBindingEnergy protOptSolvent aggreSim amberAnalyzer protEvolverBinding database_phipsi protFolder sideChainRandomizer simpleFolder dielectricFit foldingBindingEnergy ligandBindingEnergy solvationPDB bindingEnergy triadFinder protMover z_aligner y_aligner

all : $(LIB_TARGETS) $(EXEC_TARGETS)

lib : FORCE
	cd obj && $(MAKE) libprotcad.a

sideChainRandomizer : FORCE
	cd obj && $(MAKE) sideChainRandomizer

acidMutator : FORCE
	cd obj && $(MAKE) acidMutator
	
chainBindingEnergy : FORCE
	cd obj && $(MAKE) chainBindingEnergy
	
foldingBindingEnergy : FORCE
	cd obj && $(MAKE) foldingBindingEnergy
	
aggreSim : FORCE
	cd obj && $(MAKE) aggreSim
	
database_phipsi : FORCE
	cd obj && $(MAKE) database_phipsi
	
intraSoluteEnergy : FORCE
	cd obj && $(MAKE) intraSoluteEnergy
	
bindingEnergy : FORCE
	cd obj && $(MAKE) bindingEnergy
	
protFolder : FORCE
	cd obj && $(MAKE) protFolder
	
simpleFolder : FORCE
	cd obj && $(MAKE) simpleFolder

triadFinder : FORCE
	cd obj && $(MAKE) triadFinder

dielectricFit : FORCE
	cd obj && $(MAKE) dielectricFit
	
protDock : FORCE
	cd obj && $(MAKE) protDock
	
mergeComplex : FORCE
	cd obj && $(MAKE) mergeComplex
	
ligandBindingEnergy : FORCE
	cd obj && $(MAKE) ligandBindingEnergy

protEvolverBinding : FORCE
	cd obj && $(MAKE) protEvolverBinding

amberAnalyzer : FORCE
	cd obj && $(MAKE) amberAnalyzer

structFinder : FORCE
	cd obj && $(MAKE) structFinder

dielectric : FORCE
	cd obj && $(MAKE) dielectric

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
	rm -f bin/$(ARCH)/$(EXEC_TARGETS)
	cd obj && $(MAKE) zzyxyz

realclean: FORCE
	$(MAKE) clean
	rm -f $(LIBDIR)/*.a 
FORCE:
