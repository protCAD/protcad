// filename: superRotamerIterator.h

#include "assert.h"
#include <iostream.h>
#include "typedef.h"
#include "protein.h"

#ifndef SUPERROTAMERITERATOR_H
#define SUPERROTAMERITERATOR_H

class superRotamerIterator
{
public:
	superRotamerIterator(chainPosition* _pChainPosition);
	~superRotamerIterator();

	superRotamerIterator& operator++ (int _x);	
	
private:
	void initialize();
	void increment();
	bool updateCurrentAllowedResidueIndex();
	bool updateCurrentAllowedRotamerIndex();

public:
	UInt last() { return failFlag; }
	allowedResidue* getAllowedResiduePointer() {return pItsCurrentAllowedResidue;}
	rotamer* getAllowedRotamerPointer() {return pItsCurrentAllowedRotamer;}
	UInt getAllowedResidueIndex() {return itsCurrentAllowedResidueIndex;}
	UInt getAllowedRotamerIndex() {return itsCurrentAllowedRotamerIndex;}
	
private:
	chainPosition* pInputChainPosition;	
	allowedResidue* pItsCurrentAllowedResidue;
	rotamer* pItsCurrentAllowedRotamer;

	UInt itsCurrentAllowedResidueIndex;
	UInt itsCurrentAllowedRotamerIndex;

	UInt numAllowedResiduesInCurrent;
	UInt numAllowedRotamersInCurrent;
	UInt numBranchpointsInCurrent;
	UInt failFlag;
};

#endif
