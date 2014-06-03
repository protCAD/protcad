// filename: deadEndIterator.h

#include "assert.h"
#include <iostream.h>
#include "typedef.h"
#include "protein.h"

#ifndef DEAD_END_ITERATOR_H
#define DEAD_END_ITERATOR_H

class deadEndIterator
{
public:
	deadEndIterator(protein* _pProtein);
	~deadEndIterator();

	deadEndIterator& operator++ (int _x);	
private:
	void initialize(UInt _chain, UInt _residue);
	void increment();
	bool updateCurrentResidueIndex();
	bool updateCurrentAllowedResidueIndex();
	bool updateCurrentBranchpointIndex();
	bool updateCurrentRotamerIndex();
public:
	UInt last() { return failFlag; }

	protein* getProteinPointer()			{return pInputProtein;}
	chain* getChainPointer()			{return pItsCurrentChain;}
	chainPosition* getChainPositionPointer()	{return pItsCurrentChainPosition;}
	residue* getResiduePointer()			{return pItsCurrentResidue;}
	allowedResidue* getAllowedResiduePointer()	{return pItsCurrentAllowedResidue;}

	UInt getChainIndex()		{return itsCurrentChainIndex;}
	UInt getResidueIndex()		{return itsCurrentResidueIndex;}
	UInt getAllowedResidueIndex()	{return itsCurrentAllowedResidue;}
	UInt getBranchpointIndex()	{return itsBranchpointIndex;}
	UInt getRotamerIndex()		{return itsRotamerIndex;}
	
private:
	protein* pInputProtein;	
	chain*   pItsCurrentChain;
	residue* pItsCurrentResidue;
	chainPosition* pItsCurrentChainPosition;
	allowedResidue* pItsCurrentAllowedResidue;

	UInt itsCurrentChainIndex;
	UInt itsCurrentResidueIndex;
	UInt itsCurrentAllowedResidue;
	UInt itsCurrentBranchpointIndex;
	UInt itsCurrentRotamerIndex;

	UInt numChainsInCurrent;
	UInt numResiduesInCurrent;
	UInt numAllowedResiduesInCurrent;
	UInt numBranchpointsInCurrent;
	UInt numRotamersInCurrent;
	UInt failFlag;
};

#endif
