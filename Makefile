export TOP=$(PWD)
export srcdir=$(TOP)/src
export projdir=$(TOP)/projects
export SVMTROOT=$(TOP)/svmt
export SVMTINCLUDE=$(SVMTROOT)/include
export GALIBROOT=$(TOP)/galib245
export GALIBINCLUDE=$(GALIBROOT)
export FOXROOT=$(TOP)/fox
export FOXINCLUDE=$(FOXROOT)/include
export TNTINCLUDE=$(TOP)

export ARCH=$(MACHTYPE)

export LIBDIR=$(TOP)/lib/$(ARCH)

export CXX = g++
export F77 = g77
export MAKE = make

LIB_TARGETS = galib fox lib lib_graphics
EXEC_TARGETS = chiralTest primalSheet chiralTransfer primalSoup nativeDesigner \
	trypsin glycoHBondEnergy2 primalHelixLoopHelix chiralInteraction z_aligner \
	dimerBuilder trimerBuilder originator y_aligner mutationDocker thrRotamer\
	serAntiDimer serDimer gcn4mutator homoOligomer helphipsiscan turnBuilder DFscpacker \
	ms1Dimer helixangle deberdesign helixEval fiber ramachandran patchCalc \
	curveHelix pentamerBuilder backboneMerge coords evenMC mapSequence mapPhiPsi \
	m2Modeler calcHelixAngle m2Fit changeLoop skeleton mutateGFPSurface rotOpthelix rotOpt \
	mutantMaker helixMeasure attractEnergy repulseEnergy indigestible digestible \
	gastricAnalyzer counter positionAnalyzer denovoBuilder denovoBuilderSym test Energy \
	savingStateTest indigestibleSym digestibleSym pkaAve digestIndex resdigestIndex \
	indigestibleEvo patternFinder chiralFlipper bbOpt protOpt condOpt fiveElements \
	ionicEnergy acidMutator folding protEvolver dFinder protDock \
	mergeComplex oligOpt structFinder structShaper pointMutator_ins \
	impSolvent dielectric burialSurface carbonDensity intraSoluteEnergy chainBindingEnergy \
	protOptSolvent foldingC aggreSim amberAnalyzer protEvolverBinding peptideBondSASA \
	traj_deformation database_phipsi protFolder sideChainRandomizer simpleFolder protRandomizer \
	dielectricFit foldingBindingEnergy ligandBindingEnergy solvationPDB bindingEnergy triadFinder

all : $(LIB_TARGETS) $(EXEC_TARGETS)

svmt : FORCE  
	cd $(SVMTROOT)/src && $(MAKE)
	ranlib $(SVMTROOT)/lib/*.a
	cp $(SVMTROOT)/lib/*.a $(LIBDIR)

galib : FORCE
	cd $(GALIBROOT) && $(MAKE) lib
	cp $(GALIBROOT)/ga/libga.a $(LIBDIR)

fox : FORCE
	cd $(FOXROOT) && ./configure && make

lib : FORCE
	cd obj && $(MAKE) libprotcad.a

lib_graphics : FORCE
	cd obj && $(MAKE) libprotcadgraphics.a

skeleton : FORCE
	cd obj && $(MAKE) skeleton

sideChainRandomizer : FORCE
	cd obj && $(MAKE) sideChainRandomizer
	
carbonDensity : FORCE
	cd obj && $(MAKE) carbonDensity
	
chainBindingEnergy : FORCE
	cd obj && $(MAKE) chainBindingEnergy
	
foldingBindingEnergy : FORCE
	cd obj && $(MAKE) foldingBindingEnergy
	
aggreSim : FORCE
	cd obj && $(MAKE) aggreSim
	
traj_deformation : FORCE
	cd obj && $(MAKE) traj_deformation
	
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

solvationPDB : FORCE
	cd obj && $(MAKE) solvationPDB

triadFinder : FORCE
	cd obj && $(MAKE) triadFinder

wtOptimizer : FORCE
	cd obj && $(MAKE) wtOptimizer
	
protRandomizer : FORCE
	cd obj && $(MAKE) protRandomizer
	
dielectricFit : FORCE
	cd obj && $(MAKE) dielectricFit
	
protDock : FORCE
	cd obj && $(MAKE) protDock
	
peptideBondSASA : FORCE
	cd obj && $(MAKE) peptideBondSASA
	
foldingC : FORCE
	cd obj && $(MAKE) foldingC

mergeComplex : FORCE
	cd obj && $(MAKE) mergeComplex
	
ligandBindingEnergy : FORCE
	cd obj && $(MAKE) ligandBindingEnergy

protEvolverBinding : FORCE
	cd obj && $(MAKE) protEvolverBinding

protEvolverHomo : FORCE
	cd obj && $(MAKE) protEvolverHomo
	
amberAnalyzer : FORCE
	cd obj && $(MAKE) amberAnalyzer

atomMap : FORCE
	cd obj && $(MAKE) atomMap

functions : FORCE
	cd obj && $(MAKE) functions

protMover : FORCE
	cd obj && $(MAKE) protMover

structFinder : FORCE
	cd obj && $(MAKE) structFinder

impSolvent : FORCE
	cd obj && $(MAKE) impSolvent

burialSurface : FORCE
	cd obj && $(MAKE) burialSurface

dielectric : FORCE
	cd obj && $(MAKE) dielectric

oligOpt : FORCE
	cd obj && $(MAKE) oligOpt
	
protOptSolvent : FORCE
	cd obj && $(MAKE) protOptSolvent

protEvolver : FORCE
	cd obj && $(MAKE) protEvolver

smartMutator : FORCE
	cd obj && $(MAKE) smartMutator

folding : FORCE
	cd obj && $(MAKE) folding

structShaper : FORCE
	cd obj && $(MAKE) structShaper

fiveElements : FORCE
	cd obj && $(MAKE) fiveElements

mutateGFPSurface : FORCE
	cd obj && $(MAKE) mutateGFPSurface

coilBuilder : FORCE
	cd obj && $(MAKE) coilBuilder
	
mutantMaker : FORCE
	cd obj && $(MAKE) mutantMaker

positionAnalyzer : FORCE
	cd obj && $(MAKE) positionAnalyzer

denovoBuilder : FORCE
	cd obj && $(MAKE) denovoBuilder

denovoBuilderSym : FORCE
	cd obj && $(MAKE) denovoBuilderSym

counter : FORCE
	cd obj && $(MAKE) counter

protOpt : FORCE
	cd obj && $(MAKE) protOpt

protOpt0 : FORCE
	cd obj && $(MAKE) protOpt0

dFinder : FORCE
	cd obj && $(MAKE) dFinder

ionicEnergy : FORCE
	cd obj && $(MAKE) ionicEnergy

condOpt : FORCE
	cd obj && $(MAKE) condOpt

helixMeasure : FORCE
	cd obj && $(MAKE) helixMeasure

digestIndex : FORCE
	cd obj && $(MAKE) digestIndex

resdigestIndex : FORCE
	cd obj && $(MAKE) resdigestIndex

patternFinder : FORCE
	cd obj && $(MAKE) patternFinder

chiralFlipper : FORCE
	cd obj && $(MAKE) chiralFlipper

bbOpt : FORCE
	cd obj && $(MAKE) bbOpt

rotOpthelix : FORCE
	cd obj && $(MAKE) rotOpthelix

pkaAve : FORCE
	cd obj && $(MAKE) pkaAve

savingStateTest : FORCE
	cd obj && $(MAKE) savingStateTest
	
test : FORCE
	cd obj && $(MAKE) test

rotOpt : FORCE
	cd obj && $(MAKE) rotOpt

indigestible : FORCE
	cd obj && $(MAKE) indigestible
	
indigestibleEvo : FORCE
	cd obj && $(MAKE) indigestibleEvo

indigestibleSym : FORCE
	cd obj && $(MAKE) indigestibleSym

gastricAnalyzer : FORCE
	cd obj && $(MAKE) gastricAnalyzer

digestible : FORCE
	cd obj && $(MAKE) digestible

digestibleSym : FORCE
	cd obj && $(MAKE) digestibleSym

attractEnergy : FORCE
	cd obj && $(MAKE) attractEnergy

repulseEnergy : FORCE
	cd obj && $(MAKE) repulseEnergy

changeLoop : FORCE
	cd obj && $(MAKE) changeLoop

calcHelixAngle : FORCE
	cd obj && $(MAKE) calcHelixAngle

m2Fit : FORCE
	cd obj && $(MAKE) m2Fit

m2Modeler : FORCE
	cd obj && $(MAKE) m2Modeler

Energy : FORCE
	cd obj && $(MAKE) Energy

lowestRotCUDA : FORCE
	cd obj && $(MAKE) lowestRotCUDA

iQuickpick : FORCE
	cd obj && $(MAKE) iQuickpick

mapSequence : FORCE
	cd obj && $(MAKE) mapSequence

mapPhiPsi : FORCE
	cd obj && $(MAKE) mapPhiPsi

coords : FORCE
	cd obj && $(MAKE) coords

evenMC : FORCE
	cd obj && $(MAKE) evenMC

backboneMerge : FORCE
	cd obj && $(MAKE) backboneMerge

pentamerBuilder : FORCE
	cd obj && $(MAKE) pentamerBuilder

curveHelix : FORCE
	cd obj && $(MAKE) curveHelix

patchCalc : FORCE
	cd obj && $(MAKE) patchCalc

fiber : FORCE
	cd obj && $(MAKE) fiber

helixEval : FORCE
	cd obj && $(MAKE) helixEval

helixangle : FORCE
	cd obj && $(MAKE) helixangle

deberdesign : FORCE
	cd obj && $(MAKE) deberdesign

deber : FORCE
	cd obj && $(MAKE) deber

ms1Dimer : FORCE
	cd obj && $(MAKE) ms1Dimer

serAntiDimer : FORCE
	cd obj && $(MAKE) serAntiDimer

serDimer : FORCE
	cd obj && $(MAKE) serDimer

DFscpacker : FORCE
	cd obj && $(MAKE) DFscpacker

turnBuilder : FORCE
	cd obj && $(MAKE) turnBuilder

helphipsiscan : FORCE
	cd obj && $(MAKE) helphipsiscan

gcn4mutator : FORCE
	cd obj && $(MAKE) gcn4mutator

homoOligomer : FORCE
	cd obj && $(MAKE) homoOligomer

thrRotamer : FORCE
	cd obj && $(MAKE) thrRotamer

mutationDocker : FORCE
	cd obj && $(MAKE) mutationDocker

downHillSimplex : FORCE
	cd obj && $(MAKE) downHillSimplex

originator : FORCE
	cd obj && $(MAKE) originator

trimerBuilder : FORCE
	cd obj && $(MAKE) trimerBuilder

dimerBuilder : FORCE
	cd obj && $(MAKE) dimerBuilder

nitrileProjection : FORCE
	cd obj && $(MAKE) nitrileProjection

y_aligner : FORCE
	cd obj && $(MAKE) y_aligner

z_aligner : FORCE
	cd obj && $(MAKE) z_aligner

zHelixVectorAligner : FORCE
	cd obj && $(MAKE) zHelixVectorAligner

repack : FORCE
	cd obj && $(MAKE) repack

nativeDesigner : FORCE
	cd obj && $(MAKE) nativeDesigner

glycoHBondEnergy2 : FORCE
	cd obj && $(MAKE) glycoHBondEnergy2

primalSoup : FORCE
	cd obj && $(MAKE) primalSoup

chiralTransfer : FORCE
	cd obj && $(MAKE) chiralTransfer

primalSheet : FORCE
	cd obj && $(MAKE) primalSheet

chiralTest : FORCE
	cd obj && $(MAKE) chiralTest

primalHelixLoopHelix : FORCE
	cd obj && $(MAKE) primalHelixLoopHelix

pointMutator : FORCE
	cd obj && $(MAKE) pointMutator

pointMutator_ins : FORCE
	cd obj && $(MAKE) pointMutator_ins

acidMutator : FORCE
	cd obj && $(MAKE) acidMutator

chiralInteraction : FORCE
	cd obj && $(MAKE) chiralInteraction

glycoDockerGridRef : FORCE
	cd obj && $(MAKE) glycoDockerGridRef

trypsin : FORCE
	cd obj && $(MAKE) trypsin

ramachandran : FORCE
	cd obj && $(MAKE) ramachandran

clean: FORCE
	rm -f bin/$(ARCH)/$(EXEC_TARGETS)
	cd obj && $(MAKE) zzyxyz

realclean: FORCE
	$(MAKE) clean
	rm -f $(LIBDIR)/*.a 
FORCE:
