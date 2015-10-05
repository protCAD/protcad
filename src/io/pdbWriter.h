#include <string.h>
//#include "svmt.h"

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

#ifndef LIGAND_H
class ligand;
#endif

#ifndef ATOM_ITERATOR_H
#include "atomIterator.h"
#endif

#ifndef RESIDUEITERATOR_H
#include "residueIterator.h"
#endif

unsigned int pdbWriter(vector<protein*> _protVec, vector<ligand*> _ligVec, const string& _pdbFile);
unsigned int pdbWriter(protein* _pProtein, const string& _filename);
unsigned int pdbWriter(ligand* _pLigand, const string& _filename);

void renumberAtoms(protein* _pProtein);
void renumberAtoms(ligand* _pLigand);
#endif
