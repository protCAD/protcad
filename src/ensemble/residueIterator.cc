#include "residueIterator.h"

residueIterator::residueIterator(protein* _pProtein)
{	
	pInputProtein = _pProtein;
	initialize();
}

residueIterator::~residueIterator()
{
}

void residueIterator::initialize()
{
	failFlag = 0;
	itsCurrentChainIndex = 0;
	itsCurrentResidueIndex = 0;
	
// need to add some bounds checking here to make sure that each of these is
// fully defined upon instantiation...
// Note: Initialization in conditions of "if" loop are intended in this case....

	if( (numChainsInCurrent = pInputProtein->itsChains.size()) )
	{	pItsCurrentChain   = (pInputProtein->itsChains)[itsCurrentChainIndex];
		if( (numResiduesInCurrent = pItsCurrentChain->itsResidues.size()) )
		{	pItsCurrentResidue = (pItsCurrentChain->itsResidues)[itsCurrentResidueIndex];
		}
		else
		{	pItsCurrentResidue = 0;
		}
	}
	else
	{	pItsCurrentChain = 0;
		pItsCurrentResidue = 0;
	}
}

residueIterator& residueIterator::operator++ (int _x)
{	increment();
//	cout << "Incrementing..";
//	cout << " failFlag = " << failFlag << endl;
	return *this;
}

void residueIterator::increment()
{
	if ( updateCurrentResidueIndex() )
	{	return;
	}
	else if ( updateCurrentChainIndex() )
	{	return;
	}
	else
	{	failFlag = 1;
		return;
	}
}

bool residueIterator::updateCurrentResidueIndex()
{	if (itsCurrentResidueIndex < numResiduesInCurrent-1)
        {       itsCurrentResidueIndex++;
                pItsCurrentResidue = (pItsCurrentChain->itsResidues)[itsCurrentResidueIndex];
                return true;
        }
	else return false;
}

bool residueIterator::updateCurrentChainIndex()
{	if (itsCurrentChainIndex < numChainsInCurrent-1)
	{	itsCurrentChainIndex++;
		pItsCurrentChain = (pInputProtein->itsChains)[itsCurrentChainIndex];
		numResiduesInCurrent = pItsCurrentChain->itsResidues.size();
		itsCurrentResidueIndex = 0;
		pItsCurrentResidue = (pItsCurrentChain->itsResidues)[itsCurrentResidueIndex];
		return true;
	}
	else return false;
}
