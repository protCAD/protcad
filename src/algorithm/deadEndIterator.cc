#include "deadEndIterator.h"

deadEndIterator::deadEndIterator(protein* _pProtein)
{//	cout << "deadEndIterator constructor called:"
//	     << "deadEndIterator::deadEndIterator()" << endl;
	pItsProtein = _pProtein;
	initialize();
}

deadEndIterator::~deadEndIterator()
{//	cout<< "deadEndIterator destructor called " << endl;
}

deadEndIterator& deadEndIterator::operator++ (int _x)
{	increment();
	return *this;	
}

void deadEndIterator::initialize()
{
	failFlag = 0;
	itsCurrentChainIndex = 0;
	itsCurrentResidueIndex = 0;
	itsCurrentBranchpointIndex = 0;
	itsCurrentRotamerIndex = 0;

	if ( (numChainsInCurrent = pInputProtein->itsChains.size()) )
	{	pItsCurrentChain = (pInputProtein->itsChains)[itsCurrentChainIndex];
