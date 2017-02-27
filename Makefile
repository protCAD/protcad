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

TARGETS = example

.SUFFIXES: .cc .o .h .a

LIB_TARGETS = lib

LIB_CC_OBJECTS = ran1.o ran.o point.o treeNode.o atom.o atomIterator.o residue.o chain.o residueTemplate.o allowedResidue.o secondaryStructure.o chainPosition.o residueIterator.o chainModBuffer.o molecule.o protein.o ensemble.o CMath.o generalio.o pdbData.o pdbReader.o pdbWriter.o amberVDW.o aaBaseline.o amberElec.o rotamer.o rotamerLib.o annealer.o PDBAtomRecord.o PDBInterface.o ruler.o line.o lineSegment.o unitSphere.o helixPropensity.o parse.o ramachandranMap.o

DEFS = -DHAVE_OPENGL=1 -D__STL_USE_EXCEPTIONS

FLAG_OPTMAX = -Wall -oFast -ffast-math -ftree-vectorize -march=native -mtune=native -pipe -msse3 -Wno-deprecated -std=gnu++11

CFLAGS = $(FLAG_OPTMAX) $(DEFS)

INC_BASE = -I$(SRCDIR)/ensemble -I$(SRCDIR)/io \
-I$(SRCDIR)/math -I$(SRCDIR)/database -I$(SRCDIR)/algorithm \
-I$(TNTINCLUDE)

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

install : $(LIB_TARGETS)
	@echo export PROTCADDIR=$(PROTCADDIR) >> ~/.bashrc
	@echo export PATH=$(PATH):$(PROTCADDIR):$(PROTCADDIR)/bin >> ~/.bashrc

lib : libprotcad.a

libprotcad.a : $(LIB_CC_OBJECTS)
	cd $(OBJDIR) && ar rv libprotcad.a $?
	cd $(OBJDIR) && ranlib libprotcad.a

example : libprotcad.a example.cc
	cd $(OBJDIR) && $(CXX) $(CFLAGS) $^ -o $@ $(INC_BASE) $(LIB_BASE)
	cd $(OBJDIR) && mv $@ $(BINDIR)

$(LIB_CC_OBJECTS): %.o: %.cc %.h
	$(CXX) -c $(CFLAGS) $(INC_BASE) $< -o $@
	mv $@ $(OBJDIR)

clean: 
	rm -f $(OBJDIR)/*.o 
	rm -f $(OBJDIR)/*.a
	cd $(BINDIR) && rm -f $(TARGETS)
