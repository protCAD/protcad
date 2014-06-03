#include "superRotamerIterator.h"

superRotamerIterator::superRotamerIterator(protein* _pProtein)
{	
	pInputProtein = _pProtein;
	initialize();
}

superRotamerIterator::~superRotamerIterator()
{
}

void superRotamerIterator::initialize()
{
	failFlag = 0;
	itsCurrentAllowedResidueIndex = 0;
	
// need to add some bounds checking here to make sure that each of these is
// fully defined upon instantiation...
// Note: Initialization in conditions of "if" loop are intended in this case....

	if( (numAllowedResiduesInCurrent = pInputChainPosition->itsAllowedResidues.size()) )
	{	pItsCurrentAllowedResidue   = &(pInputChainPosition->itsAllowedResidues)[itsCurrentAllowedResidueIndex];
		if( (numAllowedRotamersInCurrent = pItsCurrentAllowedResidue->itsAllowedRotamers.size()) )
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

superRotamerIterator& superRotamerIterator::operator++ (int _x)
{	increment();
	return *this;
}

void superRotamerIterator::increment()
{
	if ( updateCurrentAllowedResidueIndex() )
	{	return;
	}
	else if ( updateCurrentAllowedRotamerIndex() )
	{	return;
	}
	else
	{	failFlag = 1;
		return;
	}
}

bool superRotamerIterator::updateCurrentAllowedResidueIndex()
{	if (itsCurrentAllowedResidueIndex < numAllowedResiduesInCurrent-1)
        {       itsCurrentAllowedResidueIndex++;
                itsCurrentAllowedResidueIdentity = (pItsCurrentChainPosition->itsResidues)[itsCurrentResidueIndex];
                return true;
        }
	else return false;
}

bool superRotamerIterator::updateCurrentChainIndex()
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
