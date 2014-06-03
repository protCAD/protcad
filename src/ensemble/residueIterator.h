// filename: residueIterator.h

#include "assert.h"
#include <iostream>
#include "typedef.h"
#include "protein.h"

#ifndef RESIDUEITERATOR_H
#define RESIDUEITERATOR_H

class residueIterator
{
public:
	residueIterator(protein* _pProtein);
	~residueIterator();

	residueIterator& operator++ (int _x);	
	void initialize();
	void increment();
	bool updateCurrentResidueIndex();
	bool updateCurrentChainIndex();
	UInt last() { return failFlag; }

	protein* getProteinPointer() {return pInputProtein;}
	chain* getChainPointer() {return pItsCurrentChain;}
	residue* getResiduePointer() {return pItsCurrentResidue;}
	UInt getChainIndex() {return itsCurrentChainIndex;}
	UInt getResidueIndex() {return itsCurrentResidueIndex;}
	
private:
	protein* pInputProtein;	
	chain*   pItsCurrentChain;
	residue* pItsCurrentResidue;

	UInt itsCurrentChainIndex;
	UInt itsCurrentResidueIndex;

	UInt numChainsInCurrent;
	UInt numResiduesInCurrent;
	UInt failFlag;
	
	bool heavyOnly;
};

#endif
