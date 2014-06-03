#include "secondaryStructure.h"

ramachandranMap secondaryStructure::theMap;

secondaryStructure::secondaryStructure()
{	itsSecondaryStructureIndex = -1;
}

secondaryStructure::secondaryStructure(int _indexInChain)
{
	itsIndexInChain = _indexInChain;
}

//deep copy
secondaryStructure::secondaryStructure(const secondaryStructure& _rhs)
{
	itsStartingResidueIndex = _rhs.itsStartingResidueIndex;
	itsEndingResidueIndex = _rhs.itsEndingResidueIndex;
	itsSecondaryStructureIndex = _rhs.itsSecondaryStructureIndex;
	itsIndexInChain = _rhs.itsIndexInChain;
}

secondaryStructure::~secondaryStructure()
{
}

void secondaryStructure::assign(chain* _pTheChain)
{
	assignViaPhiPsi(_pTheChain);
}

void secondaryStructure::assignViaPhiPsi(chain* _pTheChain)
{
//	cout << " assignViaPhiPsi ";
	if (itsIndexInChain == 0 || UInt(itsIndexInChain + 1) == _pTheChain->getNumResidues())
	{
//		cout << itsIndexInChain; 
		itsSecondaryStructureIndex = 4;
	}
	else
	{
//		cout << itsIndexInChain; 
		double phi = _pTheChain->getPhi(itsIndexInChain);
		double psi = _pTheChain->getPsi(itsIndexInChain);
		int tempint = theMap.getSecondaryStructureIndex(phi,psi);
//		cout << " phi " << phi << " psi " << psi << " index " << tempint << endl;
		itsSecondaryStructureIndex = tempint;
	}
}
