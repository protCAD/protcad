// filename: atomIterator.h

#include "assert.h"
#include <iostream>
#include "typedef.h"
#include "protein.h"

#ifndef ATOMITERATOR_H
#define ATOMITERATOR_H
//#warning "atomIterator.h read in"
class atomIterator
{
public:
	atomIterator(protein* _pProtein);
	~atomIterator();

	atomIterator& operator++ (int _x);	
	void initializeLigand();
	void initialize();

private:
	void increment();
	bool updateCurrentAtomIndex();
	bool updateCurrentResidueIndex();
	bool updateCurrentChainIndex();
public:
	UInt last() { return failFlag; }
	void setHeavyAtomsOnly() {heavyOnly = true;}
	void setAllAtoms() {heavyOnly = false;}

	protein* getProteinPointer() {return pInputProtein;}
	chain* getChainPointer() {return pItsCurrentChain;}
	residue* getResiduePointer() {return pItsCurrentResidue;}
	atom* getAtomPointer() {return pItsCurrentAtom;}
	UInt getChainIndex() {return itsCurrentChainIndex;}
	UInt getResidueIndex() {return itsCurrentResidueIndex;}
	UInt getAtomIndex() {return itsCurrentAtomIndex;}
	
private:
	protein* pInputProtein;	
	chain*   pItsCurrentChain;
	residue* pItsCurrentResidue;
	atom*    pItsCurrentAtom;

	UInt itsCurrentChainIndex;
	UInt itsCurrentResidueIndex;
	UInt itsCurrentAtomIndex;

	UInt numChainsInCurrent;
	UInt numResiduesInCurrent;
	UInt numAtomsInCurrent;
	UInt failFlag;

	bool hetatmFlag;
	
	bool heavyOnly;
};

#endif
