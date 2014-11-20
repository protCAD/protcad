export TOP=$(PWD)
export srcdir=$(TOP)/src
export TNTINCLUDE=$(TOP)/src
export BINDIR=$(TOP)/bin
export LIBDIR=$(TOP)/lib

export CXX = g++
export F77 = g77
export MAKE = make

LIB_TARGETS = lib
EXEC_TARGETS = example

all : $(LIB_TARGETS) $(EXEC_TARGETS)

lib : FORCE
	cd obj && $(MAKE) libprotcad.a

example : FORCE
	cd obj && $(MAKE) example

clean: FORCE
	cd bin && rm $(EXEC_TARGETS)
	cd obj && $(MAKE) cleanlibs

FORCE:
