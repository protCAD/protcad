#include <string.h>

#ifndef PDB_DATA_H
#include "pdbData.h"
#endif

#ifndef MOLECULE_H
#include "molecule.h"
#endif

#ifndef PDB_WRITER_H
#define PDB_WRITER_H

#ifndef PROTEIN_H
class protein;
#endif

#ifndef ATOM_ITERATOR_H
#include "atomIterator.h"
#endif

#ifndef RESIDUEITERATOR_H
#include "residueIterator.h"
#endif

unsigned int pdbWriter(protein* _pProtein, const string& _filename);

void renumberAtoms(protein* _pProtein);
#endif
