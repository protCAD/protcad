// filename: protein.cc
// contents: class protein implementation

#include "protein.h"
bool protein::messagesActive = false;
typedef vector<UInt>::iterator iterUINT;
typedef vector<int>::iterator iterINT;
typedef vector<vector<int> >:: iterator iterINTVEC;

bool protein::calcSelfEnergy = true;
UInt protein::howMany = 0;
UInt protein::itsSolvationParam = 0;

protein::protein() : molecule()
{
#ifdef PROTEIN_DEBUG
	cout<< "default protein constructor called" << endl;
	cout << itsName << endl;
#endif
	itsChains.resize(0);
	setMoleculeType(1);
	resetAllBuffers();
	itsLastModifiedChain = -1;
	initializeModificationMethods();
	howMany++;
}

protein::protein(const string& _name) : molecule(_name)
{
#ifdef PROTEIN_DEBUG
	cout<< "protein constructor for " << _name << " called" << endl;
#endif
	itsChains.resize(0);
    itsIndependentChainsMap.resize(0);
	itsChainLinkageMap.resize(0);
	itsLastModifiedChain = -1;
	setMoleculeType(1);
	resetAllBuffers();
	initializeModificationMethods();
	howMany++;
}

protein::protein(const protein& _rhs)
{	for (UInt i=0; i<_rhs.itsChains.size(); i++)
	{	add(new chain( *(_rhs.itsChains[i])));
	}
	itsName = _rhs.itsName;
	itsLastModifiedChain = _rhs.itsLastModifiedChain;
	itsLastModificationMethod = _rhs.itsLastModificationMethod;
	initializeModificationMethods();
	setMoleculeType(1);
	itsChainLinkageMap = _rhs.itsChainLinkageMap;
	itsIndependentChainsMap = _rhs.itsIndependentChainsMap;
}

protein::~protein()
{
#ifdef PROTEIN_DEBUG
	cout << "protein destructor called " << endl;
#endif
	for(UInt i=0; i<itsChains.size(); i++)
	{	delete itsChains[i];
	}
    howMany--;
}

void protein::resetAllBuffers()
{
	itsLastModifiedChain = -1;
	itsLastModificationMethod = -1;
}

/*chain& protein::getChain(UInt _chainIndex)
{
	if (_chainIndex < itsChains.size())
		return itsChains(_chainIndex);
	else
	{
		cout << "ERROR:  chain " << _chainIndex << " not found." << endl;
		return itsChains(0);
	}
}*/

void protein::add(chain* _pChain)
{
	itsChains.push_back(_pChain);
	UInt index = itsChains.size() - 1;
	itsIndependentChainsMap.push_back(index);
	vector<int> tempIntVec;
	tempIntVec.resize(0);
	tempIntVec.push_back(-1);
	itsChainLinkageMap.push_back(tempIntVec);
}

void protein::initializeModificationMethods()
{
	bool (protein::* pFunc)(ran&);
	pFunc = &protein::performRandomMutation;
	itsModificationMethods[0] = pFunc;
	pFunc = &protein::performRandomRotamerChange;
	itsModificationMethods[1] = pFunc;
	pFunc = &protein::performRandomRotamerRotation;
	itsModificationMethods[2] = pFunc;
	itsModificationMethods[3] = 0;
	itsModificationMethods[4] = 0;
}
void protein::removeChain(UInt _chainIndex)
{
	
	for (UInt i = 0; i < itsChains.size(); i++)
	{
		if (i > _chainIndex)
		{
			cout << "map " << itsIndependentChainsMap[i] << endl;
			itsIndependentChainsMap[i] = itsIndependentChainsMap[i]-1;
		}
	}
	delete itsChains[_chainIndex];
	itsChains.resize(itsChains.size()-1);
	
}

void protein::symmetryLinkChainAtoB(UInt _aIndex, UInt _bIndex)
{
	//cout << "Before Linkage" << endl;
	//printAllLinkageInfo();
	// first, remove index of chain A from itsIndependentChainsMap
	UInt Aindex=0;
	UInt Bindex=0;
	bool AfoundInIndependentList = false;
	bool BfoundInIndependentList = false;
	// find chain A in ChainsMap
	for (UInt i=0; i<itsIndependentChainsMap.size(); i++)
	{	if ( _aIndex == itsIndependentChainsMap[i])
		{	Aindex = i;
			AfoundInIndependentList = true;
		}
	}

	if (AfoundInIndependentList)
	{
		// remove the cell from itsIndependentChainsMap
		iterUINT firstUINT;
		firstUINT = itsIndependentChainsMap.begin();
		if ( Aindex < itsIndependentChainsMap.size())
		{       itsIndependentChainsMap.erase(firstUINT + Aindex);
		}

		// also remove the corresponding cell from itsChainLinkageMap, if it exists,
		// and retain the data
		vector<int> subLinkages = itsChainLinkageMap[Aindex];
		iterINTVEC firstINTVEC;
		firstINTVEC = itsChainLinkageMap.begin();
		if ( Aindex < itsChainLinkageMap.size())
		{	itsChainLinkageMap.erase(firstINTVEC + Aindex);
		}

		if (subLinkages.size() == 1 && subLinkages[0] == -1)
		{	// its not linked to anthing -- ignore
		}
		else
		{
		// we need to so something about sublinkages here.... recursion?
			for (UInt i=0; i<subLinkages.size(); i++)
			{
				symmetryLinkChainAtoB(subLinkages[i], _bIndex);
			}
		}
	}

	// now, find out where chain B resides in the ChainsMap
	for (UInt i=0; i<itsIndependentChainsMap.size(); i++)
	{	if ( _bIndex == itsIndependentChainsMap[i])
		{	Bindex = i;
			BfoundInIndependentList = true;
		}
	}
	if (BfoundInIndependentList)
	{
		// is this the first one?
		if ( itsChainLinkageMap[Bindex].size() == 1 && itsChainLinkageMap[Bindex][0] == -1)
		{	itsChainLinkageMap[Bindex][0] = _aIndex;
		}
		else
		{	itsChainLinkageMap[Bindex].push_back(_aIndex);
		}
	}
	else
	{	cout << "Error in  protein::symmetryLinkChainAtoB" << endl;
		cout << "Cannot find the independent chain " << _bIndex << " specified" << endl;
	}
	// now we need to copy active residue information from chain B to chain A

	//cout << "After Linkage" << endl;
	//cout << endl;
	if (messagesActive)	cout << "Chain " << _aIndex << " linked to chain " << _bIndex << endl;
	printAllLinkageInfo();
}

void protein::printAllLinkageInfo()
{
	if (messagesActive) cout << "Linkage Info:" << endl << "--------------" << endl;
	UInt numIndependentChains = itsIndependentChainsMap.size();
	if (messagesActive) cout << "Independent Chains" << endl;
	for (UInt i=0; i<numIndependentChains; i++)
	{	if (messagesActive) cout << itsIndependentChainsMap[i] << " " ;
	}
	if (messagesActive) cout << endl;

	UInt sizeOfChainLinkageMap = itsChainLinkageMap.size();

	if (messagesActive) cout << "Linked Chains" << endl;
	UInt largest = 0;
	for (UInt i=0; i<sizeOfChainLinkageMap; i++)
	{	if (itsChainLinkageMap[i].size() > largest)
		{	largest = itsChainLinkageMap[i].size();
		}
	}

	UInt counter = 0;
	while (counter <= largest)
	{
		for (UInt i=0; i<sizeOfChainLinkageMap; i++)
		{
			if (counter < itsChainLinkageMap[i].size())
			{	if (itsChainLinkageMap[i][counter] != -1)
				{	if (messagesActive) cout << itsChainLinkageMap[i][counter];
					if (messagesActive) cout << " ";
				}
				else
				{	if (messagesActive) cout << "- ";
				}
			}
		}
		if (messagesActive) cout << endl;
		counter++;
	}
	if (messagesActive) cout << "Done with linkage info" << endl;
}

void protein::activateForRepacking(const UInt _chainIndex, const UInt _residueIndex)
{
	if (_chainIndex < itsChains.size())
	{	if( !(itsChains[_chainIndex]->activateForRepacking(_residueIndex)))
		{	//cout << "Activation Failure on chain " << _chainIndex;
			//cout << ", residue " << _residueIndex << endl;
		}
	}
	else
	{	cout << "Error from protein::activateForRepacking" << endl;
		cout << "chain index out of bounds :" << _chainIndex << endl;
	}

	// find chain in independent chain list
	UInt indChainIndex = 0;
	for (UInt i=0; i<itsIndependentChainsMap.size(); i++)
	{	if (itsIndependentChainsMap[i] == _chainIndex)
		{	indChainIndex = i;
		}
	}

	UInt numSymLinkedChains = itsChainLinkageMap[indChainIndex].size();
	if ( /*numSymLinkedChains != 1 &&*/ itsChainLinkageMap[indChainIndex][0] != -1)
	{       for (UInt i=0; i<numSymLinkedChains;i++)
		{
			itsChains[itsChainLinkageMap[indChainIndex][i]]->activateForRepacking(_residueIndex);
		}
	}
}

void protein::activateForRepacking(const UInt _chainIndex, const UInt _start, const UInt _end)
{
	for (UInt i=_start; i<=_end; i++)
	{	activateForRepacking(_chainIndex,i);
	}
}

void protein::activateAllForRepacking(const UInt _chainIndex)
{
	if (_chainIndex < itsChains.size())
	{	UInt tempint = itsChains[_chainIndex]->getNumResidues();
		for (UInt i=0; i<tempint; i++)
		{	activateForRepacking(_chainIndex,i);
		}
	}
	else
	{	cout << "Error from protein::activateAllForRepacking" << endl;
		cout << "chain index out of bounds :" << _chainIndex << endl;
	}
}

UIntVec protein::getResAllowed (const UInt _chainIndex, const UInt _residueIndex)
{
	UIntVec allowedResidues;
	if (_chainIndex < itsChains.size() && _chainIndex >= 0)
	{
		allowedResidues = itsChains[_chainIndex]->getResAllowed(_residueIndex);
		return allowedResidues;
	}
	else
	{
		cout << "Error from protein::getResAllowed ... chain value passed is " << _chainIndex << endl;
		allowedResidues.resize(0);
		return allowedResidues; // return null list
	}
}


void protein::setResNotAllowed(const UInt _chainIndex, const UInt _residueIndex, const UInt _residueType)
{
	if (_chainIndex < itsChains.size())
	{	itsChains[_chainIndex]->setResNotAllowed(_residueIndex, _residueType);
	}
	else
	{	cout << "Error from protein::setResNotAllowed" << endl;
		cout << "chain index out of bounds :" << _chainIndex << endl;
	}

	// find chain in independent chain list
	UInt indChainIndex = 0;
	for (UInt i=0; i<itsIndependentChainsMap.size(); i++)
	{	if (itsIndependentChainsMap[i] == _chainIndex)
		{	indChainIndex = i;
		}
	}

	UInt numSymLinkedChains = itsChainLinkageMap[indChainIndex].size();
	if (/* numSymLinkedChains != 1 &&*/ itsChainLinkageMap[indChainIndex][0] != -1)
	{       for (UInt i=0; i<numSymLinkedChains;i++)
		{
			itsChains[itsChainLinkageMap[indChainIndex][i]]->setResNotAllowed(_residueIndex, _residueType);
		}
	}
}

void protein::setResAllowed(const UInt _chainIndex, const UInt _residueIndex, const UInt _residueType)
{
	if (_chainIndex < itsChains.size())
	{	itsChains[_chainIndex]->setResAllowed(_residueIndex, _residueType);
	}
	else
	{	cout << "Error from protein::setResAllowed" << endl;
		cout << "chain index out of bounds :" << _chainIndex << endl;
	}
	// find chain in independent chain list
	UInt indChainIndex = 0;
	for (UInt i=0; i<itsIndependentChainsMap.size(); i++)
	{	if (itsIndependentChainsMap[i] == _chainIndex)
		{	indChainIndex = i;
		}
	}

	UInt numSymLinkedChains = itsChainLinkageMap[indChainIndex].size();
	if ( /*numSymLinkedChains != 1 &&*/ itsChainLinkageMap[indChainIndex][0] != -1)
	{       for (UInt i=0; i<numSymLinkedChains;i++)
		{
			itsChains[itsChainLinkageMap[indChainIndex][i]]->setResAllowed(_residueIndex, _residueType);
		}
	}
}

void protein::setOnlyCharged(const UInt _chainIndex, const UInt _residueIndex)
{
	// only allow LYS, ARG, GLU, HIS and ASP

        setResNotAllowed(_chainIndex,_residueIndex,0); //ALA
        setResNotAllowed(_chainIndex,_residueIndex,2); //ASN
        setResNotAllowed(_chainIndex,_residueIndex,4); //CYS
        setResNotAllowed(_chainIndex,_residueIndex,5); //GLN
        setResNotAllowed(_chainIndex,_residueIndex,7); //GLY
        setResNotAllowed(_chainIndex,_residueIndex,9); //ILE
        setResNotAllowed(_chainIndex,_residueIndex,10); //LEU
        setResNotAllowed(_chainIndex,_residueIndex,12); //MET
        setResNotAllowed(_chainIndex,_residueIndex,13); //PHE
        setResNotAllowed(_chainIndex,_residueIndex,14); //PRO
        setResNotAllowed(_chainIndex,_residueIndex,15); //SER
        setResNotAllowed(_chainIndex,_residueIndex,16); //THR
        setResNotAllowed(_chainIndex,_residueIndex,17); //TRP
        setResNotAllowed(_chainIndex,_residueIndex,18); //TYR
        setResNotAllowed(_chainIndex,_residueIndex,19); //VAL
		return;
}

void protein::setListNotAllowed(const UInt _chainIndex, const UInt _residueIndex, const vector <UInt> _typeIndexVector)
{
	for (UInt i = 0; i < _typeIndexVector.size(); i++)
		setResNotAllowed(_chainIndex,_residueIndex,_typeIndexVector[i]);
	return;
}

void protein::setOnlyHydrophilic(const UInt _chainIndex, const UInt _residueIndex)
{
		setResNotAllowed(_chainIndex,_residueIndex,0); //ALA
		setResNotAllowed(_chainIndex,_residueIndex,4); //CYS
		setResNotAllowed(_chainIndex,_residueIndex,7); //GLY
		setResNotAllowed(_chainIndex,_residueIndex,9); //ILE
		setResNotAllowed(_chainIndex,_residueIndex,10); //LEU
		setResNotAllowed(_chainIndex,_residueIndex,12); //MET
		setResNotAllowed(_chainIndex,_residueIndex,13); //PHE
		setResNotAllowed(_chainIndex,_residueIndex,14); //PRO
		setResNotAllowed(_chainIndex,_residueIndex,15); //SER
		setResNotAllowed(_chainIndex,_residueIndex,16); //THR
		setResNotAllowed(_chainIndex,_residueIndex,17); //TRP
		setResNotAllowed(_chainIndex,_residueIndex,18); //TYR
		setResNotAllowed(_chainIndex,_residueIndex,19); //VAL
}

void protein::setOnlyROCHydrophobic(const UInt _chainIndex, const UInt _residueIndex)
{
		setResNotAllowed(_chainIndex,_residueIndex,0); //ALA
		setResNotAllowed(_chainIndex,_residueIndex,1); //ARG
		setResNotAllowed(_chainIndex,_residueIndex,2); //ASN
		setResNotAllowed(_chainIndex,_residueIndex,3); //ASP
		setResNotAllowed(_chainIndex,_residueIndex,4); //CYS
		setResNotAllowed(_chainIndex,_residueIndex,5); //GLN
		setResNotAllowed(_chainIndex,_residueIndex,6); //GLU
		setResNotAllowed(_chainIndex,_residueIndex,7); //GLY
		setResNotAllowed(_chainIndex,_residueIndex,8); //HIS
		setResNotAllowed(_chainIndex,_residueIndex,11); //LYS
		setResNotAllowed(_chainIndex,_residueIndex,12); //MET
		setResNotAllowed(_chainIndex,_residueIndex,14); //PRO
		setResNotAllowed(_chainIndex,_residueIndex,15); //SER
		setResNotAllowed(_chainIndex,_residueIndex,16); //THR
		setResNotAllowed(_chainIndex,_residueIndex,17); //TRP
		setResNotAllowed(_chainIndex,_residueIndex,18); //TYR
}

void protein::setOnlyNativeIdentity(const UInt _chainIndex, const UInt _residueIndex)
{
	itsChains[_chainIndex]->setOnlyNativeIdentity(_residueIndex);
	// find chain in independent chain list
	UInt indChainIndex = 0;
	for (UInt i=0; i<itsIndependentChainsMap.size(); i++)
	{	if (itsIndependentChainsMap[i] == _chainIndex)
		{	indChainIndex = i;
		}
	}

	UInt numSymLinkedChains = itsChainLinkageMap[indChainIndex].size();
	if ( /*numSymLinkedChains != 1 &&*/ itsChainLinkageMap[indChainIndex][0] != -1)
	{       for (UInt i=0; i<numSymLinkedChains;i++)
		{
			itsChains[itsChainLinkageMap[indChainIndex][i]]->setOnlyNativeIdentity(_residueIndex);
		}
	}
}

void protein::setAllHydrogensOn(const bool _hydrogensOn)
{
	for(UInt i=0; i<itsChains.size();i++)
	{
		itsChains[i]->setAllHydrogensOn(_hydrogensOn);
	}
}

void protein::setAllPolarHydrogensOn(const bool _polarHydrogensOn)
{
	// note: if true, (all) hydrogensOn should be false
	for(UInt i=0;i<itsChains.size();i++)
	{
		itsChains[i]->setAllPolarHydrogensOn(_polarHydrogensOn);
	}
}

double protein::netCharge()
{
    double nCharge = 0.0;
    for(UInt i=0;i<itsChains.size();i++)
    {
        nCharge += itsChains[i]->netCharge();
    }
    return nCharge;
}


void protein::stripToGlycine()
{
	for (UInt theChain =0; theChain<itsChains.size(); theChain++)
	{	for (UInt res=0; res< itsChains[theChain]->getNumResidues(); res++)
		{	mutateWBC(theChain,res,7);
		}
	}
}

int protein::mutate(vector <int> _position, UInt _resType)
{
	if (_position[1] >=0 && _position[1] < (int)itsChains.size())
	{
		mutateWBC(_position[1], _position[2], _resType);
		return 1;
	}
	else cout << "ERROR in protein::mutate(vector<int> _position, UInt _resType) ... chain position passed to function is out of range." << endl;
	return -1;
}

void protein::mutateWBC(const UInt _chainIndex, const UInt _resIndex, const UInt _aaIndex)
{	if(_chainIndex < itsChains.size())
	{	itsChains[_chainIndex]->mutate(_resIndex,_aaIndex);
		itsChains[_chainIndex]->commitLastMutation();
		// find chain in independent chain list
		UInt indChainIndex = 0;
		for (UInt i=0; i<itsIndependentChainsMap.size(); i++)
		{	if (itsIndependentChainsMap[i] == _chainIndex)
			{	indChainIndex = i;
			}
		}

		UInt numSymLinkedChains = itsChainLinkageMap[indChainIndex].size();
		if ( /*numSymLinkedChains != 1 &&*/ itsChainLinkageMap[indChainIndex][0] != -1)
		{       for (UInt i=0; i<numSymLinkedChains;i++)
			{
				itsChains[itsChainLinkageMap[indChainIndex][i]]->mutate(_resIndex,_aaIndex);
				itsChains[itsChainLinkageMap[indChainIndex][i]]->commitLastMutation();
			}
		}
	}
}

void protein::mutate(const UInt _chainIndex, const UInt _resIndex, const UInt _aaIndex)
{	itsChains[_chainIndex]->mutate(_resIndex,_aaIndex);
	itsChains[_chainIndex]->commitLastMutation();
	// find chain in independent chain list
	UInt indChainIndex = 0;
	for (UInt i=0; i<itsIndependentChainsMap.size(); i++)
	{	if (itsIndependentChainsMap[i] == _chainIndex)
		{	indChainIndex = i;
		}
	}

	UInt numSymLinkedChains = itsChainLinkageMap[indChainIndex].size();
	if (itsChainLinkageMap[indChainIndex][0] != -1)
	{       for (UInt i=0; i<numSymLinkedChains;i++)
		{
			itsChains[itsChainLinkageMap[indChainIndex][i]]->mutate(_resIndex,_aaIndex);
			itsChains[itsChainLinkageMap[indChainIndex][i]]->commitLastMutation();
		}
	}
}

void protein::makeAtomSilent(const UInt _chainIndex, const UInt _resIndex, const UInt _atomIndex)
{
	if ( _chainIndex >=0 && _chainIndex < itsChains.size())
	{
		itsChains[_chainIndex]->makeAtomSilent(_resIndex, _atomIndex);
	}
	else
		cout << "ERROR in protein::makeAtomSilent(...) ...\n\t" << _chainIndex << " value for chain index is out of range." << endl;
	return;
}

void protein::makeResidueSilent(const UInt _chainIndex, const UInt _resIndex)
{
	if ( _chainIndex >=0 && _chainIndex < itsChains.size())
	{
		itsChains[_chainIndex]->makeResidueSilent(_resIndex);
	}
	else
		cout << "ERROR in protein::makeAtomSilent(...) ...\n\t" << _chainIndex << " value for chain index is out of range." << endl;
	return;
}

vector <int> protein::getLastModification()
{
	vector <int> position;
	position.resize(0);
	position.push_back(itsLastModifiedChain);
	position.push_back(itsChains[itsLastModifiedChain]->getLastModificationPosition());
	return position;
}

int protein::modify(ran& _ran, vector <int> _position)
{
    if (chainPosition::getHowMany() == 0)
    {
        cout << "Error reported from protein::modify()" << endl;
        cout << "No chain positions have been activated for modification!" << endl;
        return -2;
    }

    int modificationMethod = chooseModificationMethod(_ran);
    if (modificationMethod >= 0)
    {
		itsLastModificationMethod = modificationMethod;
		switch (modificationMethod)
		{
			case 0:
				if(performRandomMutation(_ran, _position)) return 1;
				break;
			case 1:
				if(performRandomRotamerChange(_ran, _position)) return 1;
				break;
			case 2:
				if(performRandomRotamerRotation(_ran, _position)) return 1;
				break;
			default:
				cout << "ERROR or ABORT in protein::modify(_ran, _position)" << endl;
				return -1;
		}
    }
    return -1;
}

int protein::modify(ran& _ran)
{
    if (chainPosition::getHowMany() == 0)
    {
		cout << "Error reported from protein::modify()" << endl;
		cout << "No chain positions have been activated for modification!" << endl;
        return -2;
    }
	int modificationMethod = chooseModificationMethod(_ran);
	if (modificationMethod >= 0)
	{
        if( (this->*itsModificationMethods[modificationMethod])(_ran) )
		{
			itsLastModificationMethod = modificationMethod;
			return 1;
		}
	}
	return -1;
}

bool protein::performRandomMutation(ran& _ran)
{	int chainToModify = chooseTargetChain(_ran);
	if (messagesActive) cout << "CHAIN " << chainToModify << " ";
	if (messagesActive) cout << "MUT ";
	if (chainToModify >= 0)
	{
		vector<chainModBuffer> theBuffers;
		theBuffers = itsChains[chainToModify]->performRandomMutation(_ran);
		/*
		cout << "ChainModBuffers for chain " << chainToModify << endl;
		theBuffers[0].printAll();
		theBuffers[1].printAll();
		*/
		if ( theBuffers[0].containsData() )
		{	itsLastModifiedChain = chainToModify;
			itsLastModificationMethod = 0;
			// find chain in independent chain list
			UInt indChainIndex = 0;
			for (UInt i=0; i<itsIndependentChainsMap.size(); i++)
			{	if (int(itsIndependentChainsMap[i]) == chainToModify)
				{	indChainIndex = i;
				}
			}

			UInt numSymLinkedChains = itsChainLinkageMap[indChainIndex].size();
			if (/* numSymLinkedChains != 1 &&*/ itsChainLinkageMap[indChainIndex][0] != -1)
			{       for (UInt i=0; i<numSymLinkedChains;i++)
				{
					itsChains[itsChainLinkageMap[indChainIndex][i]]->repeatModification(theBuffers[1]);
				}
			}
			return true;
		}
		else
		{
		//	cout << "Protein level has detected abort..." << endl;
		}

	}
	else
	{	cout << "Error from protein::performRandomMutation" << endl;
		cout << "Protein contains no chains" << endl;
	}
	resetAllBuffers();
	return false;
}

bool protein::performRandomMutation(ran& _ran, vector <int> _position)
{
	int chainToModify = _position[1];
	if (messagesActive) cout << "CHAIN " << chainToModify << " ";
	if (messagesActive) cout << "MUT ";
	if (chainToModify >= 0)
	{
		vector<chainModBuffer> theBuffers;
		theBuffers = itsChains[chainToModify]->performRandomMutation(_ran, _position);
		if (theBuffers[0].containsData())
		{
			itsLastModifiedChain = chainToModify;
			itsLastModificationMethod = 0;
			UInt indChainIndex = 0;
            for (UInt i=0; i<itsIndependentChainsMap.size(); i++)
            {   if (int(itsIndependentChainsMap[i]) == chainToModify)
                {   indChainIndex = i;
                }
            }

            UInt numSymLinkedChains = itsChainLinkageMap[indChainIndex].size();
            if (/* numSymLinkedChains != 1 &&*/ itsChainLinkageMap[indChainIndex][0] != -1)
            {       for (UInt i=0; i<numSymLinkedChains;i++)
                {
                    itsChains[itsChainLinkageMap[indChainIndex][i]]->repeatModification(theBuffers[1]);
                }
            }
            return true;
        }
        else
        {
        //  cout << "Protein level has detected abort..." << endl;
        }

    }
    else
    {   cout << "Error from protein::performRandomMutation" << endl;
        cout << "Protein contains no chains" << endl;
    }
    resetAllBuffers();
    return false;
}

void protein::setupSystem(ran& _ran)
{
	// try mutation first
	// if mutation fails, try rotamer changes
	vector<UInt> targetID;
	for (UInt i=0; i < itsIndependentChainsMap.size(); i++)
	{
        UInt iChain = itsIndependentChainsMap[i];
        // get number of active residues in chain
        vector<chainPosition*> cpVector;
        vector<UInt> activeMap;
        UIntVec allowedRes;
        cpVector = itsChains[iChain]->getChainPositionVector();
        activeMap = itsChains[iChain]->getRepackActivePositionMap();
        UInt numActive = activeMap.size();
        for (UInt j=0; j< numActive; j++)
        {   if (cpVector[activeMap[j]])
            {   allowedRes = cpVector[activeMap[j]]->getResAllowed();
                UInt tempTargetID = allowedRes[int(_ran.getNext() * (allowedRes.size()-1))];
                mutate(iChain,activeMap[j],tempTargetID);
            }
        }
	}
/*
           UInt numSymLinkedChains = itsChainLinkageMap[iChain].size();
           if ( numSymLinkedChains != 1 && itsChainLinkageMap[iChain][0] != -1)
           {       for (UInt i=0; i<numSymLinkedChains;i++)
               {
                   itsChains[itsChainLinkageMap[iChain][i]]->repeatModification(theBuffers[1]);
               }
           }
*/
}

int protein::chooseModificationMethod(ran& _ran)
{
	// in this 'bare bones' version, every modification
	// method is given equal weight.  Not necessarily
	// the best implementation, but useful for testing
	int methodSize = 3;
	if (methodSize != 0)
	{
		int i =  int(_ran.getNext() * methodSize);
		return i;
	}
	return -1;
}

void protein::acceptModification()
{	if (itsLastModifiedChain >=0)
	{	switch (itsLastModificationMethod)
		{
			case 0:	commitLastMutation();
					break;
			case 1:	commitLastRotamerChange();
					break;
			case 2:	commitLastRotamerRotation();
					break;
			default:	cout << "oops! hit default value in protein::acceptModification";
					cout << endl << "something is definitely wrong here." << endl;
					break;
		}
	}
}

void protein::rejectModification()
{	if (itsLastModifiedChain >=0)
	{	switch (itsLastModificationMethod)
		{
			case 0:	undoLastMutation();
					break;
			case 1:	undoLastRotamerChange();
					break;
			case 2:	undoLastRotamerRotation();
					break;
			default:	cout <<"oops! hit default value in protein::rejectModification";
					cout << endl << "something is definitely wrong here." << endl;
					break;
		}
	}
}

void protein::finishProteinBuild()
{
	for (unsigned int i=0; i< itsChains.size(); i++)
	{
		itsChains[i]->rebuildResiduesInChain();
		itsChains[i]->finishChainBuild();
	}
/*
	vector < UIntVec > artificialPos;
	artificialPos.resize(0);
	for (UInt i = 0; i < itsChains.size(); i++)
	{
		UIntVec artificialChainAndRes;
		for (UInt j = 0; j < itsChains[i]->getNumResidues(); j++)
		{
			if (itsChains[i]->isArtificiallyBuilt(j))
			{
				artificialChainAndRes.push_back(i);
				artificialChainAndRes.push_back(j);
				artificialPos.push_back(artificialChainAndRes);
				artificialChainAndRes.resize(0);
				//cout << "Flag squirted " << i << " "<< j<< endl;
			}
		}
	}
    if (artificialPos.size() <= 10 && artificialPos.size() != 0)
    {
        optimizeRotamers(artificialPos);
    }
    else
    {
        if (artificialPos.size() != 0) cout << "Too many residues to optimize with this algorithm.\n OPTIMIZATION ABORTED." << endl;
    }
 */
    return;
}

void protein::listSecondaryStructure()
{
	if (messagesActive) cout << "SECONDARY STRUCTURE" << endl;
	for (unsigned int i=0; i< itsChains.size(); i++)
	{	if (messagesActive) cout << "CHAIN " << i << endl;
		itsChains[i]->listSecondaryStructure();
	}
}

void protein::listDihedrals()
{
	if (messagesActive) cout << "DIHEDRAL ANGLES" << endl;
	for (unsigned int i=0; i< itsChains.size(); i++)
	{	if (messagesActive) cout << "CHAIN " << i << endl;
		itsChains[i]->listDihedrals();
	}
}

vector <int> protein::chooseNextTargetPosition(ran& _ran)
{
	vector <int> position;
	position.resize(0);

	int chainPos = chooseTargetChain(_ran);
	int resPos = itsChains[chainPos]->chooseNextTargetPosition(_ran);

	position.push_back(chainPos);
	position.push_back(resPos);

	return position;
}

UInt protein::chooseNextMutationIdentity(ran& _ran, vector <int> _position)
{
	if (_position[1] >=0 && _position[1] < (int)itsChains.size())
	{
		return itsChains[_position[1]]->chooseNextMutationIdentity(_ran, _position);
	}
	else cout << "ERROR in chooseNextMutationIdentity at protein level ... chain ID is out of range." << endl;
	return 0;
}

int protein::chooseTargetChain(ran& _ran)
{	int chainSize = itsIndependentChainsMap.size();
	if (chainSize != 0)
	{	return int(itsIndependentChainsMap[int(_ran.getNext() * chainSize)]);
	}
	return -1;
}

void protein::commitLastMutation()
{	if (itsLastModifiedChain >=0 && itsLastModificationMethod == 0)
	{	itsChains[itsLastModifiedChain]->commitLastMutation();
		// find chain in independent chain list
		UInt indChainIndex = 0;
		for (UInt i=0; i<itsIndependentChainsMap.size(); i++)
		{	if (int(itsIndependentChainsMap[i]) == itsLastModifiedChain)
			{	indChainIndex = i;
			}
		}

		UInt numSymLinkedChains = itsChainLinkageMap[indChainIndex].size();
		if (/* numSymLinkedChains != 1 &&*/ itsChainLinkageMap[indChainIndex][0] != -1)
		{       for (UInt i=0; i<numSymLinkedChains;i++)
			{
				itsChains[itsChainLinkageMap[indChainIndex][i]]->commitLastMutation();
			}
		}
		resetAllBuffers();
	}
	else
	{	cout << "Error reported by protein::commitLastMutation()";
		cout << endl << "No last modified chain" << endl;
	}
	if (messagesActive) cout << "ACCEPTED " << endl;
}

void protein::undoLastMutation()
{	if (itsLastModifiedChain >=0 && itsLastModificationMethod == 0)
	{	itsChains[itsLastModifiedChain]->undoLastMutation();
		// find chain in independent chain list
		UInt indChainIndex = 0;
		for (UInt i=0; i<itsIndependentChainsMap.size(); i++)
		{	if (int(itsIndependentChainsMap[i]) == itsLastModifiedChain)
			{	indChainIndex = i;
			}
		}

		UInt numSymLinkedChains = itsChainLinkageMap[indChainIndex].size();
		if (/* numSymLinkedChains != 1 &&*/ itsChainLinkageMap[indChainIndex][0] != -1)
		{       for (UInt i=0; i<numSymLinkedChains;i++)
			{
				itsChains[itsChainLinkageMap[indChainIndex][i]]->undoLastMutation();
			}
		}
		resetAllBuffers();
	}
	else
	{	cout << "Error reported by protein::undoLastMutation()";
		cout << endl << "No last modified chain" << endl;
	}
	if (messagesActive) cout << "REJECTED " << endl;
}

vector <chainPosition*> protein::getChainPositionVector(const UInt _chain)
{	if (itsChains[_chain])
	{
		return itsChains[_chain]->getChainPositionVector();
	}
	else
	{	cout << "Error reported by protein::getChainPositionVector()";
		cout << endl << "chain " << _chain << " is invalid specifier" << endl;
		cout << "Returning null vector - further behavior unpredictable!!!" << endl;
	}
	chainPosition* nullCP = 0;
	vector <chainPosition*> nullVec;
	nullVec.push_back(nullCP);
	return nullVec;
}

double protein::getVolume(UInt _method)
{
	double itsVolume = 0.0;
	for (UInt i = 0; i < itsChains.size(); i++)
	{
		itsVolume += itsChains[i]->getVolume(_method);
	}
	return itsVolume;
}

// get the interaction energy of a particular modified position with the rest of the system, taking into account positions
// symmetry linked to the passed position
double protein::getPositionEnergy(UIntVec _pos)
{
	vector <int> pos;
	pos.resize(0);
	for (UInt i = 0; i < _pos.size(); i ++)
	{
		pos.push_back((int)_pos[i]);
	}
	return getPositionEnergy(pos);
}

double protein::getPositionEnergy(vector <int> _position)
{
   // bool calcSurfEnergy = false;
   // if (solvation::getItsScaleFactor() != 0.0) calcSurfEnergy = true;
   // if (calcSurfEnergy) getSurfaceArea();
	if (_position.size() == 2) // ensemble number not in
	{

		vector <int> temp;
		temp.resize(0);
		temp.push_back(0);
		temp.push_back(_position[0]);
		temp.push_back(_position[1]);
		_position = temp;
	}

    double intraEnergy = 0.0;
    if (_position[1] >=0 && _position[1] < (int)itsChains.size())
    {

        double tempE = itsChains[_position[1]]->getPositionIntraEnergy(_position);
        //if (calcSelfEnergy)
        //{
        //    double selfEnergy = getSelfEnergy(_position[1], _position[2]);
        //    tempE -= selfEnergy;
        //}
        //if (calcSurfEnergy) tempE += getSolvationEnergyNoSASA(_position[1], _position[2]);
        if (messagesActive) cout << " intraEnergy of chain " << _position[1] << ": " << tempE << endl;
        intraEnergy += tempE;
        for (UInt i = 0; i < itsChains.size(); i++)
        {
            if ((int)i != _position[1])
            {
                tempE = itsChains[_position[1]]->getPositionInterEnergy(_position, itsChains[i]);
                if (messagesActive) cout << " interEnergy of chains " << _position[1] << " and " << i << ": " << tempE << endl;
                intraEnergy += tempE;
            }
        }
    }

/*	for (UInt i = 0; i < itsChains.size(); i ++)
	{
		for (UInt j = 0; j < getNumResidues(i); j ++)
		{
			if (i == _position[1] && j == _position[2])
			{ // skip
			}
			else
			{
				intraEnergy += getHBondEnergy(_position[1], _position[2], i, j);
			}
		}
	}
*/
    return intraEnergy;
}

int protein::setPhi(const UInt _chain, const UInt _res, double _phi)
{
	if (_chain < itsChains.size())
	{
		return itsChains[_chain]->setPhi(_res, _phi);
	}
	else
	{
		cout << "chain index out of range" << endl;
		return -1;
	}
}

int protein::setPsi(const UInt _chain, const UInt _res, double _psi)
{
	if (_chain < itsChains.size())
	{
		return itsChains[_chain]->setPsi(_res, _psi);
	}
	else
	{
		cout << "chain index out of range" << endl;
		return -1;
	}
}

int protein::setDihedral(const UInt _chainIndex, const UInt _resIndex, double _dihedral, UInt _angleType, UInt _direction)
{
	if (_chainIndex < itsChains.size())
	{
		return itsChains[_chainIndex]->setDihedral(_resIndex, _dihedral, _angleType, _direction);
	}
	else
	{
		cout << "chain index out of range" << endl;
		return -1;
	}
}

double protein::getPositionEnergy(UInt _chain, UInt _residue)
{
	vector <int> position;
	position.resize(0);
	position.push_back(-1);
	position.push_back((int)_chain);
	position.push_back((int)_residue);
	double energy=getPositionEnergy(position);

	return energy;
}

double protein::getIntraEnergy(const UInt _chain1, const UInt _residue1, const UInt _atom1, const UInt _chain2, const UInt _residue2, const UInt _atom2)
{
	if (_chain1 >=0 && _chain1 < itsChains.size() && _chain2 >= 0 && _chain2 < itsChains.size())
	{
		return itsChains[_chain1]->getInterEnergy(_residue1, _atom1, itsChains[_chain2], _residue2, _atom2);
	}
	else
	{
		cout << "ERROR in getIntraEnergy(chain1,res1,atom1,chain2,res2,atom2) ..." << endl;
		cout << "chain index of 1 or 2 out of range." << endl;
	}
	exit(1);
}

double protein::intraEnergy()
{	double intraEnergy = 0.0;
	for(UInt i=0; i<itsChains.size(); i++)
	{
#ifdef _PROTEIN_DEBUG
		if (messagesActive) cout << "intraEnergy being calculated! ";
#endif
		intraEnergy += itsChains[i]->intraEnergy();
		if (messagesActive) cout << "\n" << "intraEnergy - chain " << i << " is: " << itsChains[i]->intraEnergy() << endl;
		for(UInt j=i+1; j<itsChains.size(); j++)
		{
			intraEnergy += itsChains[i]->interEnergy(itsChains[j]);
			if (messagesActive) cout << "\n" << "interEnergy - chain " << i << "," << j << " is: " << itsChains[i]->interEnergy(itsChains[j]) << endl;
		}
	}
    //if (solvation::getItsScaleFactor() != 0.0)
    //{
     //   intraEnergy += tabulateSolvationEnergy();
    //}
	return intraEnergy;
}


//**************Default non-Redundant optimized Energy Function***************************************
double protein::protEnergy()
{
	updateEnergy();
	double Energy = 0.0;
	for(UInt i=0; i<itsChains.size(); i++)
	{
		Energy += itsChains[i]->getEnergy();
	}
	return Energy;
}

double protein::protEnergy(UInt chainIndex) //Energy of chain alone
{
	setMoved(true,0);
	updateEnergy(chainIndex);
	double Energy = itsChains[chainIndex]->getEnergy();
	return Energy;
}

void protein::updateEnergy()
{
	updateMovedDependence(0);
	updateDielectrics();
	for(UInt i=0; i<itsChains.size(); i++)
	{
		itsChains[i]->updateEnergy();
		for(UInt j=i+1; j<itsChains.size(); j++)
		{
			itsChains[i]->updateEnergy(itsChains[j]);
		}
	}
	setMoved(false,0);
}

void protein::updateEnergy(UInt chainIndex)
{
	updateMovedDependence(0);
	updateDielectrics(chainIndex);
	itsChains[chainIndex]->updateEnergy();
	setMoved(false,0);
}

void protein::updateDielectrics()
{
	for(UInt i=0; i<itsChains.size(); i++)
	{
		itsChains[i]->polarizability();
		for(UInt j=i+1; j<itsChains.size(); j++)
		{
			itsChains[i]->polarizability(itsChains[j]);
		}
	}
	calculateDielectrics();
}

void protein::updateDielectrics(UInt chainIndex)
{
	itsChains[chainIndex]->polarizability();
	calculateDielectrics(chainIndex);
}

void protein::calculateDielectrics()
{
	for(UInt i=0; i<itsChains.size(); i++)
	{
		itsChains[i]->calculateDielectrics();
	}
}

void protein::calculateDielectrics(UInt chainIndex)
{
	itsChains[chainIndex]->calculateDielectrics();
}

void protein::updateMovedDependence(UInt _EorC)
{
	for(UInt i=0; i<itsChains.size(); i++)
	{
		itsChains[i]->updateMovedDependence(_EorC);
		for(UInt j=i+1; j<itsChains.size(); j++)
		{
			itsChains[i]->updateMovedDependence(itsChains[j], _EorC);
		}
	}
}

void protein::setMoved(bool _moved, UInt _EorC)
{
	for(UInt i=0; i<itsChains.size(); i++)
	{
		itsChains[i]->setMoved(_moved,_EorC);
	}
}

double protein::protEnergy(UInt chainIndex, UInt resIndex)
{
	if (getMoved(chainIndex, resIndex, 0)){
		updateEnergy();
	}
	return itsChains[chainIndex]->getEnergy(resIndex);
}

double protein::getMedianResidueEnergy()
{
	updateEnergy();
	double median, resE;
	vector <double> resEnergies;
	for (UInt i = 0; i < itsChains.size(); i++)
	{
		for (UInt j = 0; j < itsChains[i]->itsResidues.size(); j++)
		{
			resE = itsChains[i]->getEnergy(j);
			resEnergies.push_back(resE);
		}
	}
	
	size_t size = resEnergies.size();
	sort(resEnergies.begin(), resEnergies.end());
	if (size % 2 == 0)
	{
		median = (resEnergies[size / 2 - 1] + resEnergies[size / 2]) / 2;
	}
	else
	{
		median = resEnergies[size / 2];
	}
	return median;
}

double protein::getMedianResidueEnergy(UIntVec _activeChains)
{
	updateEnergy();
	double median, resE;
	vector <double> resEnergies;
	for (UInt i = 0; i < _activeChains.size(); i++)
	{
		for (UInt j = 0; j < itsChains[_activeChains[i]]->itsResidues.size(); j++)
		{
			resE = itsChains[_activeChains[i]]->getEnergy(j);
			resEnergies.push_back(resE);
		}
	}
	
	size_t size = resEnergies.size();
	sort(resEnergies.begin(), resEnergies.end());
	if (size % 2 == 0)
	{
		median = (resEnergies[size / 2 - 1] + resEnergies[size / 2]) / 2;
	}
	else
	{
		median = resEnergies[size / 2];
	}
	return median;
}

bool protein::boltzmannEnergyCriteria(double _deltaEnergy) //calculate boltzmann probability of an energy to determine acceptance criteria
{
	bool acceptance = false;
	double Entropy = 1000000/((rand() % 1000000)+1); //generate high precision random probability as entropy
	double PiPj = pow(EU,(_deltaEnergy/residue::getKT()));
	if (PiPj < Entropy){acceptance = true;}
	return acceptance;
}

double protein::boltzmannProbabilityToEnergy(double Pi, double Pj) //calculate boltzmann Energy from a probability (Pi) compared to a reference (Pj) probability to determine acceptance criteria
{
	double KT = residue::getTemperature()*KB;
	double Energy = -KT*log(Pi/Pj);
	return Energy;
}

void protein::updateBackboneClashes()
{
	updateMovedDependence(2);
	for(UInt i=0; i<itsChains.size(); i++)
	{
		itsChains[i]->updateBackboneClashes();
		for(UInt j=i+1; j<itsChains.size(); j++)
		{
			itsChains[i]->updateBackboneClashes(itsChains[j]);
		}
	}
	setMoved(false, 2);
}

UInt protein::getNumHardBackboneClashes()
{
	updateBackboneClashes();
	UInt clashes = 0;
	for(UInt i=0; i<itsChains.size(); i++)
	{
		clashes += itsChains[i]->getBackboneClashes();
	}
	return clashes;
}

void protein::updateClashes()
{
	updateMovedDependence(1);
	for(UInt i=0; i<itsChains.size(); i++)
	{
		itsChains[i]->updateClashes();
		for(UInt j=i+1; j<itsChains.size(); j++)
		{
			itsChains[i]->updateClashes(itsChains[j]);
		}
	}
	setMoved(false, 1);
}

UInt protein::getNumHardClashes()
{
	updateClashes();
	UInt clashes = 0;
	for(UInt i=0; i<itsChains.size(); i++)
	{
		clashes += itsChains[i]->getClashes();
	}
	return clashes;
}

UInt protein::getNumHardClashes(UInt chainIndex, UInt resIndex)
{
	if (getMoved(chainIndex, resIndex,1)){
		updateClashes();
	}
	return itsChains[chainIndex]->getClashes(resIndex);
}

UInt protein::getMedianResidueNumHardClashes()
{
	updateClashes();
	UInt median, clashes;
	vector <UInt> resClashes;
	for (UInt i = 0; i < itsChains.size(); i++)
	{
		for (UInt j = 0; j < itsChains[i]->itsResidues.size(); j++)
		{
			clashes = itsChains[i]->getClashes(j);
			resClashes.push_back(clashes);
		}
	}
	
	size_t size = resClashes.size();
	sort(resClashes.begin(), resClashes.end());
	if (size % 2 == 0)
	{
		median = (resClashes[size / 2 - 1] + resClashes[size / 2]) / 2;
	}
	else
	{
		median = resClashes[size / 2];
	}
	return median;
}

double protein::getMedianResidueEnergy(UIntVec _activeChains, UIntVec _activeResidues)
{
	updateEnergy();
	double median, resE;
	vector <double> resEnergies;
	for (UInt i = 0; i < _activeChains.size(); i++)
	{
		for (UInt j = 0; j < _activeResidues.size(); j++)
		{
			resE = itsChains[_activeChains[i]]->getEnergy(_activeResidues[j]);
			resEnergies.push_back(resE);
		}
	}
	size_t size = resEnergies.size();

	sort(resEnergies.begin(), resEnergies.end());

	if (size  % 2 == 0)
	{
	  median = (resEnergies[size / 2 - 1] + resEnergies[size / 2]) / 2;
	}
	else
	{
	  median = resEnergies[size / 2];
	}
	return median;
}

void protein::updateResiduesPerTurnType()
{
	for(UInt i=0; i<itsChains.size(); i++)
	{
		itsChains[i]->updateResiduesPerTurnType();
	}
}

void protein::updateTotalNumResidues()
{
	UInt numResidues = 0;
	for(UInt i=0; i<itsChains.size(); i++)
	{
		numResidues += getNumResidues(i);
	}
	itsNumResidues = numResidues;
}

double protein::intraEnergy(UInt _chain1, UInt _chain2)
{
    if (_chain1 >=0 && _chain2 >=0 && _chain1 < itsChains.size() && _chain2 < itsChains.size() )
    {
        return itsChains[_chain1]->interEnergy(itsChains[_chain2]);
    }
    else
    {
        cout << "ERROR in intraEnergy(chain, chain) chain index out of range." << endl;
        exit(1);
    }
}

double protein::intraEnergy(const UInt _chain)
{
	if((_chain>=0) && (_chain<itsChains.size()))
	{
		return itsChains[_chain]->intraEnergy();
	}
	else
	{
		cout << "ERROR in protein::intraEnergy(_chain), exiting\n";
		exit(1);
	}
}

UInt protein::getNumResidues(UInt _chainIndex) const
{	if (_chainIndex < itsChains.size())
	{	return itsChains[_chainIndex]->getNumResidues();
	}
	else
	{	cout << "Error reparoted from protein::getNumResidues" << endl;
		cout << "chainIndex out of bounds" << endl;
	}
	return 0;
}

int protein::getIndexFromResNum(UInt _chainIndex, UInt _resnum)
{
	// is chain valid?
	if (_chainIndex < itsChains.size())
	{
		// first, take a guess and find direction to search
		int tempint = itsChains[_chainIndex]->mapResNumToChainPosition(_resnum);
		return tempint;
	}
	else
	{
		cout << "Invalid Chain Specifier: " << _chainIndex << endl;
	}
	return -1;
}

void protein::setChi(const UInt _chainIndex, const UInt _resIndex, const UInt _bpt, const UInt _chi, const double _angle)
{
	if (_chainIndex >= 0 && _chainIndex < itsChains.size())
	{
		itsChains[_chainIndex]->setChi(_resIndex, _bpt, _chi, _angle);
	}
	else
	{
		cout << "ERROR in protein::setChi(...)\n\tchain index is out of bounds ..." << endl;
		return;
	}
	    // make same change to symmetry linked chains
    UInt independentChainIndex = 0;
    for (UInt i = 0; i < itsIndependentChainsMap.size(); i++)
    {
        if (itsIndependentChainsMap[i] == _chainIndex) independentChainIndex = i;
    }
    UInt numSymmetryLinkedChains = itsChainLinkageMap[independentChainIndex].size();
    if (itsChainLinkageMap[independentChainIndex][0] != -1)
    {
        for (UInt i = 0; i < numSymmetryLinkedChains; i ++)
        {
            itsChains[itsChainLinkageMap[independentChainIndex][i]]->setChi(_resIndex, _bpt, _chi, _angle);
        }
    }

    return;
}

void protein::rotateChain(UInt _chain, const axis _axis, const double _theta)
{
	if (_chain >= 0 && _chain < itsChains.size())
	{
		itsChains[_chain]->rotate(_axis, _theta);
	}
	else
	{
		cout << "ERROR in protein::rotateChain(...)\n\tchain index is out of bounds ..." << endl;
		return;
	}
	return;
}

void protein::rotateChainRelative(UInt _chain, const axis _axis, const double _theta)
{
	if (_chain >= 0 && _chain < itsChains.size())
	{
		itsChains[_chain]->rotateRelative(_axis, _theta);
	}
	else
	{
		cout << "ERROR in protein::rotateChainRelative(...)\n\tchain index is out of bounds ..." << endl;
		return;
	}
	return;
}

void protein::translateChain(UInt _chain, const double _x, const double _y, const double _z)
{
	if (_chain >= 0 && _chain < itsChains.size())
	{
		itsChains[_chain]->translate( _x, _y, _z);
	}
	else
	{
		cout << "ERROR in protein::translateChain(...)\n\tchain index is out of bounds ..." << endl;
		return;
	}
	return;
}

void protein::translateChainR(UInt _chain, const double _x, const double _y, const double _z)
{
	if (_chain >= 0 && _chain < itsChains.size())
	{
		itsChains[_chain]->translateR( _x, _y, _z);
	}
	else
	{
		cout << "ERROR in protein::translateChain(...)\n\tchain index is out of bounds ..." << endl;
		return;
	}
	return;
}

void protein::setSidechainDihedralAngles(UInt _chainIndex, UInt _indexInChain, vector <vector <double> > Angles)
{
	if (_chainIndex >= 0 && _chainIndex < itsChains.size())
	{
		itsChains[_chainIndex]->setSidechainDihedralAngles(_indexInChain, Angles);
	}
	else
	{
		cout << "ERROR in protein::setSidechainDihedralAngles(...)\n\tchain index is out of bounds ..." << endl;
		return;
	}
	return;
}

void protein::setRelativeChi(const UInt _chainIndex, const UInt _resIndex, const UInt _bpt, const UInt _chi, const double _angle)
{
	if (_chainIndex >= 0 && _chainIndex < itsChains.size())
	{
		itsChains[_chainIndex]->setRelativeChi(_resIndex, _bpt, _chi, _angle);
	}
	else
	{
		cout << "ERROR in protein::setRelativeChi(...)\n\tchain index is out of bounds ..." << endl;
		return;
	}

	// make same change to symmetry linked chains
	UInt independentChainIndex = 0;
	for (UInt i = 0; i < itsIndependentChainsMap.size(); i++)
    {
        if (itsIndependentChainsMap[i] == _chainIndex) independentChainIndex = i;
    }
	UInt numSymmetryLinkedChains = itsChainLinkageMap[independentChainIndex].size();
    if (itsChainLinkageMap[independentChainIndex][0] != -1)
    {
        for (UInt i = 0; i < numSymmetryLinkedChains; i ++)
        {
            itsChains[itsChainLinkageMap[independentChainIndex][i]]->setRelativeChi(_resIndex, _bpt, _chi, _angle);
        }
    }

    return;
}

void protein::setRotamerNotAllowed(const UInt _chainIndex, const UInt _residueIndex, const UInt _resType, const UInt _bpt, const UInt _rotamer )
{
    if (_chainIndex >= 0 && _chainIndex < itsChains.size())
    {
        itsChains[_chainIndex]->setRotamerNotAllowed(_residueIndex, _resType, _bpt, _rotamer);
    }
    else
    {
        cout << "Error from protein::setRotamerNotAllowed" << endl;
        cout << " chain index out of bounds : " << _chainIndex << endl;
		return;
    }

	// symmetry mapping code:
	UInt independentChainIndex = 0;
    for (UInt i = 0; i < itsIndependentChainsMap.size(); i++)
    {
        if (itsIndependentChainsMap[i] == _chainIndex) independentChainIndex = i;
    }

    UInt numSymmetryLinkedChains = itsChainLinkageMap[independentChainIndex].size();
    if (itsChainLinkageMap[independentChainIndex][0] != -1)
    {
        for (UInt i = 0; i < numSymmetryLinkedChains; i ++)
        {
            itsChains[itsChainLinkageMap[independentChainIndex][i]]->setRotamerNotAllowed( _residueIndex, _resType, _bpt, _rotamer );
        }
    }

    return;
}

void protein::setCanonicalHelixRotamersOnly(const UInt _chainIndex)
{
	UIntVec activeResidues = itsChains[_chainIndex]->getActiveResidues();
	if (activeResidues.size() > 0)
	{
		for (UInt i = 0; i < activeResidues.size(); i++)
		{
			setCanonicalHelixRotamersOnly(_chainIndex, activeResidues[i]);
		}
	}
	else
		cout << "No active positions on chain " << _chainIndex << endl;

	return;
}

void protein::setCanonicalHelixRotamersOnly(const UInt _chainIndex, const UInt _resIndex)
{
	UIntVec activeResidues = itsChains[_chainIndex]->getActiveResidues();
	bool isActive = false;
	for (UInt i = 0; i < activeResidues.size(); i++) if (activeResidues[i] == _resIndex) isActive = true;
	if (isActive)
	{
		UIntVec allowedResidues = getResAllowed(_chainIndex, _resIndex);
		for (UInt i = 0; i < allowedResidues.size(); i++)
		{
			setCanonicalHelixRotamersOnly( _chainIndex, _resIndex, allowedResidues[i]);
		}
	}
	return;
}

void protein::setCanonicalHelixRotamersOnly( const UInt _chainIndex, const UInt _resIndex, const UInt _resType )
{
	UIntVec allowedResidues = getResAllowed(_chainIndex, _resIndex);
	bool resAllowed = false;
	for (UInt i = 0; i < allowedResidues.size(); i ++)
	{
			if (_resType == allowedResidues[i]) resAllowed = true;
	}
	if (resAllowed)
	{
		if (messagesActive) cout << "Setting canonical helix rotamers ONLY for chain " << _chainIndex << " pos " << _resIndex << " for " << residue::getDataBaseItem(_resType) << endl;
		for (UInt bpt = 0; bpt < residue::getNumBpt(_resType); bpt++)
		{
			UIntVec allowedRotamers;
			allowedRotamers  = itsChains[_chainIndex]->getAllowedRotamers(_resIndex,  _resType, bpt);
			for (UInt j = 0; j < allowedRotamers.size(); j ++)
			{
				if (isValidHelixRotamer(_resType, bpt, allowedRotamers[j]) )
				{
					if (messagesActive) cout << " keeping rotamer - " << allowedRotamers[j] << " ";
				}
				else
				{
					setRotamerNotAllowed ( _chainIndex, _resIndex, _resType, bpt, allowedRotamers[j]);
				}
			}
			if (messagesActive) cout << endl;
		}
	}
	else
	{
		cout << "ERROR in protein::setCanonicalHelixRotamersOnly." << endl;
		cout << " " << residue::getDataBaseItem(_resType) << " not allowed at position " << _chainIndex << " " << _resIndex << endl;
		cout << " rotamer modification aborted for this residue type." << endl;
	}
	return;
}

bool protein::isValidHelixRotamer(UInt _resType, UInt _bpt, UInt _rotamer)
{
	if (_resType == 0) return true;
	if (_resType == 1) { if (_rotamer == 37 || _rotamer == 40 || _rotamer == 43) return true;}
	if (_resType == 2) { if (_rotamer == 7 || _rotamer == 8) return true;}
	if (_resType == 3) { if (_rotamer == 7 || _rotamer == 8) return true;}
	if (_resType == 4) { if (_rotamer == 2) return true;}
	if (_resType == 5) { if (_rotamer == 9 || _rotamer == 21 || _rotamer == 22 || _rotamer == 23 || _rotamer == 28) return true;}
	if (_resType == 6) { if (_rotamer == 13 || _rotamer == 22) return true;}
	if (_resType == 7) return true;
	if (_resType == 8) { if (_rotamer == 3 || _rotamer == 6) return true;}
	if (_resType == 9) { if (_rotamer == 7 || _rotamer == 8) return true;}
	if (_resType == 10){ if (_rotamer == 3 || _rotamer == 7) return true;}
	if (_resType == 11){ if (_rotamer == 40 || _rotamer == 67) return true;}
	if (_resType == 12){ if (_rotamer == 21 || _rotamer == 23) return true;}
	if (_resType == 13){ if (_rotamer == 3 || _rotamer == 6) return true;}
	if (_resType == 14){ if (_rotamer == 1) return true;}
	if (_resType == 15){ if (_rotamer == 0 || _rotamer == 2) return true;}
	if (_resType == 16){ if (_rotamer == 0 || _rotamer == 2) return true;}
	if (_resType == 17){ if (_rotamer == 3 || _rotamer == 5 || _rotamer == 8) return true;}
	if (_resType == 18){ if (_rotamer == 3 || _rotamer == 6) return true;}
	if (_resType == 19){ if (_rotamer == 1 ) return true;}

	return false;
}

void protein::setRotamerWBC(const UInt _chainIndex, const UInt _resIndex, const UInt _bpt, const UInt _rotamer)
{
    //cout << "Setting rotamer res:"<<_resIndex<<" bpt:"<<_bpt<< " rotamer:"<<_rotamer<< endl;
    if(_chainIndex < itsChains.size())
    {   itsChains[_chainIndex]->setRotamerWithoutBuffering(_resIndex,_bpt,_rotamer);
        // find chain in independent chain list
        UInt indChainIndex = 0;
        for (UInt i=0; i<itsIndependentChainsMap.size(); i++)
        {   if (itsIndependentChainsMap[i] == _chainIndex)
            {   indChainIndex = i;
            }
        }

        UInt numSymLinkedChains = itsChainLinkageMap[indChainIndex].size();
        if (/* numSymLinkedChains != 1 &&*/ itsChainLinkageMap[indChainIndex][0] != -1)
        {       for (UInt i=0; i<numSymLinkedChains;i++)
            {
                itsChains[itsChainLinkageMap[indChainIndex][i]]->setRotamerWithoutBuffering(_resIndex,_bpt,_rotamer);
            }
        }
    }
}

void protein::setRotamer(const UInt _chainIndex, const UInt _resIndex, const UInt _bpt, const UInt _rotamer)
{
    itsChains[_chainIndex]->setRotamerWithoutBuffering(_resIndex,_bpt,_rotamer);
    // find chain in independent chain list
    UInt indChainIndex = 0;
    for (UInt i=0; i<itsIndependentChainsMap.size(); i++)
    {   if (itsIndependentChainsMap[i] == _chainIndex)
        {   indChainIndex = i;
        }
    }

    UInt numSymLinkedChains = itsChainLinkageMap[indChainIndex].size();
    if (/* numSymLinkedChains != 1 &&*/ itsChainLinkageMap[indChainIndex][0] != -1)
    {       for (UInt i=0; i<numSymLinkedChains;i++)
        {
            itsChains[itsChainLinkageMap[indChainIndex][i]]->setRotamerWithoutBuffering(_resIndex,_bpt,_rotamer);
        }
    }
}

void protein::setPolarHRotamer(const UInt _chainIndex, const UInt _resIndex, const UInt _rotamer)
{
    itsChains[_chainIndex]->setPolarHRotamerWithoutBuffering(_resIndex,_rotamer);
    // find chain in independent chain list
    UInt indChainIndex = 0;
    for (UInt i=0; i<itsIndependentChainsMap.size(); i++)
    {   if (itsIndependentChainsMap[i] == _chainIndex)
        {   indChainIndex = i;
        }
    }

    UInt numSymLinkedChains = itsChainLinkageMap[indChainIndex].size();
    if (/* numSymLinkedChains != 1 &&*/ itsChainLinkageMap[indChainIndex][0] != -1)
    {       for (UInt i=0; i<numSymLinkedChains;i++)
        {
            itsChains[itsChainLinkageMap[indChainIndex][i]]->setPolarHRotamerWithoutBuffering(_resIndex,_rotamer);
        }
    }
}

bool protein::performRandomRotamerChange(ran& _ran)
{   int chainToModify = chooseTargetChain(_ran);
    cout << "CHAIN " << chainToModify << " ";
    cout << "ROTC ";
    if (chainToModify >= 0)
    {
        vector<chainModBuffer> theBuffers;
        theBuffers = itsChains[chainToModify]->performRandomRotamerChange(_ran);
        if ( theBuffers[0].containsData() )
        {
            itsLastModifiedChain = chainToModify;
            itsLastModificationMethod = 1;
            // find chain in independent chain list
            UInt indChainIndex = 0;
            for (UInt i=0; i<itsIndependentChainsMap.size(); i++)
            {   if (int(itsIndependentChainsMap[i]) == chainToModify)
                {   indChainIndex = i;
                }
            }

            UInt numSymLinkedChains = itsChainLinkageMap[indChainIndex].size();
            if (/* numSymLinkedChains != 1 &&*/ itsChainLinkageMap[indChainIndex][0] != -1)
            {       for (UInt i=0; i<numSymLinkedChains;i++)
                {
                    itsChains[itsChainLinkageMap[indChainIndex][i]]->repeatModification(theBuffers[1]);
                }
            }
            return true;
        }
        else
        {
            //cout << "Protein level has detected abort..." << endl;
        }
    }
    else
    {   cout << "Error from protein::performRandomRotamerChange" << endl;
        cout << "Protein contains no chains" << endl;
    }
    resetAllBuffers();
    return false;
}

bool protein::performRandomRotamerChange(ran& _ran, vector <int> _position)
{   int chainToModify = _position[1];
    cout << "CHAIN " << chainToModify << " ";
    cout << "ROTC ";
    if (chainToModify >= 0)
    {
        vector<chainModBuffer> theBuffers;
        theBuffers = itsChains[chainToModify]->performRandomRotamerChange(_ran, _position);
        if ( theBuffers[0].containsData() )
        {
            itsLastModifiedChain = chainToModify;
            itsLastModificationMethod = 1;
            // find chain in independent chain list
            UInt indChainIndex = 0;
            for (UInt i=0; i<itsIndependentChainsMap.size(); i++)
            {   if (int(itsIndependentChainsMap[i]) == chainToModify)
                {   indChainIndex = i;
                }
            }

            UInt numSymLinkedChains = itsChainLinkageMap[indChainIndex].size();
            if (/* numSymLinkedChains != 1 &&*/ itsChainLinkageMap[indChainIndex][0] != -1)
            {       for (UInt i=0; i<numSymLinkedChains;i++)
                {
                    itsChains[itsChainLinkageMap[indChainIndex][i]]->repeatModification(theBuffers[1]);
                }
            }
            return true;
        }
        else
        {
            cout << "Protein level has detected abort..." << endl;
        }
    }
    else
    {   cout << "Error from protein::performRandomRotamerChange" << endl;
        cout << "Protein contains no chains" << endl;
    }
    resetAllBuffers();
    return false;
}

bool protein::performRandomRotamerRotation(ran& _ran)
{   int chainToModify = chooseTargetChain(_ran);
    cout << "CHAIN " << chainToModify << " ";
    cout << "RR ";
    if (chainToModify >= 0)
    {
        vector<chainModBuffer> theBuffers;
        theBuffers = itsChains[chainToModify]->performRandomRotamerRotation(_ran);
        if ( theBuffers[0].containsData() )
        {   itsLastModifiedChain = chainToModify;
            itsLastModificationMethod = 2;
            // find chain in independent chain list
            UInt indChainIndex = 0;
            for (UInt i=0; i<itsIndependentChainsMap.size(); i++)
            {   if (int(itsIndependentChainsMap[i]) == chainToModify)
                {   indChainIndex = i;
                }
            }

            UInt numSymLinkedChains = itsChainLinkageMap[indChainIndex].size();
            if (/* numSymLinkedChains != 1 &&*/ itsChainLinkageMap[indChainIndex][0] != -1)
            {       for (UInt i=0; i<numSymLinkedChains;i++)
                {
                    itsChains[itsChainLinkageMap[indChainIndex][i]]->repeatModification(theBuffers[1]);
                }
            }
            return true;
        }
        else
        {
            //cout << "Protein level has detected abort..." << endl;
        }
    }
    else
    {   cout << "Error from protein::performRandomRotamerRotation" << endl;
        cout << "Protein contains no chains" << endl;
    }
    resetAllBuffers();
    return false;
}

bool protein::performRandomRotamerRotation(ran& _ran, vector <int> _position)
{   int chainToModify = _position[1];
    cout << "CHAIN " << chainToModify << " ";
    cout << "RR ";
    if (chainToModify >= 0)
    {
        vector<chainModBuffer> theBuffers;
        theBuffers = itsChains[chainToModify]->performRandomRotamerRotation(_ran, _position);
        if ( theBuffers[0].containsData() )
        {   itsLastModifiedChain = chainToModify;
            itsLastModificationMethod = 2;
            // find chain in independent chain list
            UInt indChainIndex = 0;
            for (UInt i=0; i<itsIndependentChainsMap.size(); i++)
            {   if (int(itsIndependentChainsMap[i]) == chainToModify)
                {   indChainIndex = i;
                }
            }

            UInt numSymLinkedChains = itsChainLinkageMap[indChainIndex].size();
            if (/* numSymLinkedChains != 1 &&*/ itsChainLinkageMap[indChainIndex][0] != -1)
            {       for (UInt i=0; i<numSymLinkedChains;i++)
                {
                    itsChains[itsChainLinkageMap[indChainIndex][i]]->repeatModification(theBuffers[1]);
                }
            }
            return true;
        }
        else
        {
            //cout << "Protein level has detected abort..." << endl;
        }
    }
    else
    {   cout << "Error from protein::performRandomRotamerRotation" << endl;
        cout << "Protein contains no chains" << endl;
    }
    resetAllBuffers();
    return false;
}

void protein::commitLastRotamerChange()
{   if(itsLastModifiedChain >= 0 && itsLastModificationMethod == 1)
    {   itsChains[itsLastModifiedChain]->commitLastRotamerChange();
        // find chain in independent chain list
        UInt indChainIndex = 0;
        for (UInt i=0; i<itsIndependentChainsMap.size(); i++)
        {   if (int(itsIndependentChainsMap[i]) == itsLastModifiedChain)
            {   indChainIndex = i;
            }
        }

        UInt numSymLinkedChains = itsChainLinkageMap[indChainIndex].size();
        if (/* numSymLinkedChains != 1 &&*/ itsChainLinkageMap[indChainIndex][0] != -1)
        {       for (UInt i=0; i<numSymLinkedChains;i++)
            {
                itsChains[itsChainLinkageMap[indChainIndex][i]]->commitLastRotamerChange();
            }
        }
        resetAllBuffers();
    }
    else
    {   cout << "Error reported by protein::commitLastRotamerChange()";
        cout << endl << "No last modified chain" << endl;
    }
    cout << "ACCEPTED " << endl;
}

void protein::saveCurrentState()
{
	for (UInt i = 0; i < itsChains.size(); i ++)
	{
		itsChains[i]->saveCurrentState();
	}
	return;
}

void protein::undoState()
{
	for (UInt i = 0; i < itsChains.size(); i ++)
	{
		itsChains[i]->undoState();
	}
	saveCurrentState();
	return;
}

void protein::commitState()
{
	for (UInt i = 0; i < itsChains.size(); i ++)
	{
		itsChains[i]->commitState();
	}
	return;
}

void protein::undoLastRotamerChange()
{   if (itsLastModifiedChain >=0 && itsLastModificationMethod == 1)
    {   itsChains[itsLastModifiedChain]->undoLastRotamerChange();
        // find chain in independent chain list
        UInt indChainIndex = 0;
        for (UInt i=0; i<itsIndependentChainsMap.size(); i++)
        {   if (int(itsIndependentChainsMap[i]) == itsLastModifiedChain)
            {   indChainIndex = i;
            }
        }

        UInt numSymLinkedChains = itsChainLinkageMap[indChainIndex].size();
        if (/* numSymLinkedChains != 1 &&*/ itsChainLinkageMap[indChainIndex][0] != -1)
        {       for (UInt i=0; i<numSymLinkedChains;i++)
            {
                itsChains[itsChainLinkageMap[indChainIndex][i]]->undoLastRotamerChange();
            }
        }
        resetAllBuffers();
    }
    else
    {   cout << "Error reported by protein::undoLastRotamerChange()";
        cout << endl << "No last modified chain" << endl;
    }
    cout << "REJECTED " << endl;
}

void protein::commitLastRotamerRotation()
{   if(itsLastModifiedChain >= 0 && itsLastModificationMethod == 2)
    {   itsChains[itsLastModifiedChain]->commitLastRotamerRotation();
        // find chain in independent chain list
        UInt indChainIndex = 0;
        for (UInt i=0; i<itsIndependentChainsMap.size(); i++)
        {   if (int(itsIndependentChainsMap[i]) == itsLastModifiedChain)
            {   indChainIndex = i;
            }
        }

        UInt numSymLinkedChains = itsChainLinkageMap[indChainIndex].size();
        if (/* numSymLinkedChains != 1 &&*/ itsChainLinkageMap[indChainIndex][0] != -1)
        {       for (UInt i=0; i<numSymLinkedChains;i++)
            {
                itsChains[itsChainLinkageMap[indChainIndex][i]]->commitLastRotamerRotation();
            }
        }
        resetAllBuffers();
    }
    else
    {   cout << "Error reported by protein::commitLastRotamerRotation()";
        cout << endl << "No last modified chain" << endl;
    }
    cout << "ACCEPTED " << endl;
}

void protein::undoLastRotamerRotation()
{   if (itsLastModifiedChain >=0 && itsLastModificationMethod == 2)
    {   itsChains[itsLastModifiedChain]->undoLastRotamerRotation();
        // find chain in independent chain list
        UInt indChainIndex = 0;
        for (UInt i=0; i<itsIndependentChainsMap.size(); i++)
        {   if (int(itsIndependentChainsMap[i]) == itsLastModifiedChain)
            {   indChainIndex = i;
            }
        }

        UInt numSymLinkedChains = itsChainLinkageMap[indChainIndex].size();
        if ( /*numSymLinkedChains != 1 &&*/ itsChainLinkageMap[indChainIndex][0] != -1)
        {       for (UInt i=0; i<numSymLinkedChains;i++)
            {
                itsChains[itsChainLinkageMap[indChainIndex][i]]->undoLastRotamerRotation();
            }
        }
        resetAllBuffers();
    }
    else
    {   cout << "Error reported by protein::undoLastRotamerRotation()";
        cout << endl << "No last modified chain" << endl;
    }
    cout << "REJECTED " << endl;
}

void protein::listAllowedRotamers(UInt _chain, UInt _resIndex)
{   if (_chain < itsChains.size())
    {   itsChains[_chain]->listAllowedRotamers(_resIndex);
    }
}

double protein::getRMSD(protein* _other)
{
    vector<dblVec> coord1;
    vector<dblVec> coord2;
    atomIterator theIter1(this);
    atomIterator theIter2(_other);
    atom* pAtom;
    bool first = true;
    
    // Load backbone atoms into vector for fit and alignment
    for (;!(theIter1.last());theIter1++)
    {
       pAtom = theIter1.getAtomPointer(); 
       if(pAtom->getName() == "N" || pAtom->getName() == "CA" || pAtom->getName() == "C" || pAtom->getName() == "O"){
            coord1.push_back(pAtom->getCoords());
       }
    }
    for (;!(theIter2.last());theIter2++)
    {
       pAtom = theIter2.getAtomPointer(); 
       if(pAtom->getName() == "N" || pAtom->getName() == "CA" || pAtom->getName() == "C" || pAtom->getName() == "O"){
            coord2.push_back(pAtom->getCoords());
       }
    }
    int diff = 0;
    if(coord1.size() != coord2.size()){
		if (coord2.size() < coord1.size()){ diff = coord1.size()-coord2.size(); first = true;}
		else{diff = coord2.size()-coord1.size(); first = false;}
    }
    int maxsize;
    if (first){maxsize = coord2.size();}else{maxsize = coord1.size();}
	double rotmat[9]; double centroid1[3]; double centroid2[3]; double rmsd = 0; int ierr = 0;
	int list1[maxsize]; int list2[maxsize]; int trials = 1; double bestRMSD= 1E10;
	double newCoord1[maxsize*3]; double newCoord2[maxsize*3]; double newCoord3[maxsize*3];
	
	if (diff != 0){trials = diff;}
	for (int h = 0; h < trials; h++)
	{
		for (int i=0; i<maxsize; i++)
		{	
			for (int j=0; j<3; j++)
			{
				if(first){
					newCoord1[ (i*3) + j] = coord1[i+h][j];
					newCoord2[ (i*3) + j] = coord2[i+h][j];
				}
				else{
					newCoord1[ (i*3) + j] = coord2[i+h][j];
					newCoord2[ (i*3) + j] = coord1[i+h][j];
				}
			}
			list1[i] = i+1;
			list2[i] = i+1;
		}
		
		// Calculate best fit of backbone atoms, rotation matrix and rmsd using fortran algorithm based on Machlachlan
		bestfit_(newCoord1, &maxsize, newCoord2, &maxsize, &maxsize, newCoord3, list1, list2, &rmsd, &ierr, rotmat, centroid1, centroid2);
		if (rmsd < bestRMSD){
			bestRMSD = rmsd;
		}
	}
	return bestRMSD;
}

void protein::alignToAxis(const axis _axis)
{
    vector<dblVec> coord1;
    vector<dblVec> coord2;
    dblVec coords(3);
    atomIterator theIter(this);
    atom* pAtom;
    
    // Load backbone atoms into vector for alignment to axis
    for (;!(theIter.last());theIter++)
    {
       pAtom = theIter.getAtomPointer(); 
       if (pAtom->getName() == "N" || pAtom->getName() == "CA" || pAtom->getName() == "C" || pAtom->getName() == "O" || pAtom->getName() == "CB"){
            coord2.push_back(pAtom->getCoords());
       }
    }
    
    // Find longest axis of protein
    double maxdist = -1E10, dist;
	for (UInt i = 0; i < coord2.size(); i++)
	{
		for (UInt j = i+1; j < coord2.size(); j++)
		{
			dist = CMath::distance(coord2[i], coord2[j]);
			if (dist > maxdist){maxdist = dist;}
		}
	}
	
	// Create set of points in space equally spaced of the size of the backbone on the axis
	double increment = maxdist/coord2.size();
	double coord = (maxdist/2)*-1;
	for (UInt i = 0; i < coord2.size(); i++)
	{
		if (_axis == X_axis){
			coords[0] = coord; coords[1] = 0.0; coords[2] = 0.0;
			coord1.push_back(coords);
		}
		if (_axis == Y_axis){
			coords[0] = 0.0; coords[1] = coord; coords[2] = 0.0;
			coord1.push_back(coords);
		}
		if (_axis == Z_axis){
			coords[0] = 0.0; coords[1] = 0.0; coords[2] = coord;
			coord1.push_back(coords);
		}
		coord += increment;
	}
	
    int maxsize = coord2.size();
	double rotmat[9]; double centroid1[3]; double centroid2[3]; double rmsd = 0; int ierr = 0;
	int list1[maxsize]; int list2[maxsize];
	double newCoord1[maxsize*3]; double newCoord2[maxsize*3]; double newCoord3[maxsize*3];
	for (int i=0; i<maxsize; i++)
	{	
		for (int j=0; j<3; j++)
		{
			newCoord1[ (i*3) + j] = coord1[i][j];
			newCoord2[ (i*3) + j] = coord2[i][j];
		}
		list1[i] = i+1;
		list2[i] = i+1;
	}
	bestfit_(newCoord1, &maxsize, newCoord2, &maxsize, &maxsize, newCoord3, list1, list2, &rmsd, &ierr, rotmat, centroid1, centroid2);

	// Load rotation vector into rotation matrix
	dblMat rotMat(3,3,3);
    for (UInt i=0; i<3; i++)
    {	for (UInt j=0; j<3; j++)
		{
			rotMat[i][j] = rotmat[(j*3) + i];
		}
    }
    
    // Translate protein centroid to zero, rotate and translate centroid to final axis
	for (UInt i = 0; i < getNumChains(); i++)
	{
		translateChain(i, -centroid2[0], -centroid2[1], -centroid2[2]);
		transform(i,rotMat);
		translateChain(i,centroid1[0], centroid1[1], centroid1[2]);
	}
}

void protein::translate(const UInt _index, const dblVec& _dblVec)
{
	if(_index < itsChains.size())
    {
		itsChains[_index]->translate(_dblVec);
    }
}

void protein::translate(const dblVec& _dblVec)
{
	// translate every chain
	for (UInt i = 0; i < itsChains.size(); i ++)
	{
		itsChains[i]->translate(_dblVec);
	}
}

void protein::eulerRotate(UInt _chain, const double _phi, const double _theta, const double _psi)
{
	// calculate rotation matrix
	double a11 = cos(_psi)*cos(_phi) - cos(_theta)*sin(_phi)*sin(_psi);
	double a12 = cos(_psi)*sin(_phi) + cos(_theta)*cos(_phi)*sin(_psi);
	double a13 = sin(_psi)*sin(_theta);
	double a21 = -1*sin(_psi)*cos(_phi) - cos(_theta)*sin(_phi)*cos(_psi);
	double a22 = -1*sin(_psi)*sin(_phi) + cos(_theta)*cos(_phi)*cos(_psi);
	double a23 = cos(_psi)*sin(_theta);
	double a31 = sin(_theta)*sin(_phi);
	double a32 = -1*sin(_theta)*cos(_phi);
	double a33 = cos(_theta);

	dblMat rotMatrix(3,3,0.0);
	rotMatrix[0][0] = a11;
	rotMatrix[1][0] = a12;
	rotMatrix[2][0] = a13;
	rotMatrix[0][1] = a21;
	rotMatrix[1][1] = a22;
	rotMatrix[2][1] = a23;
	rotMatrix[0][2] = a31;
	rotMatrix[1][2] = a32;
	rotMatrix[2][2] = a33;

	transform(_chain,rotMatrix);
	return;
}

void protein::eulerRotate(const double _phi, const double _theta, const double _psi)
{
	// calculate rotation matrix
	double a11 = cos(_psi)*cos(_phi) - cos(_theta)*sin(_phi)*sin(_psi);
	double a12 = cos(_psi)*sin(_phi) + cos(_theta)*cos(_phi)*sin(_psi);
	double a13 = sin(_psi)*sin(_theta);
	double a21 = -1*sin(_psi)*cos(_phi) - cos(_theta)*sin(_phi)*cos(_psi);
	double a22 = -1*sin(_psi)*sin(_phi) + cos(_theta)*cos(_phi)*cos(_psi);
	double a23 = cos(_psi)*sin(_theta);
	double a31 = sin(_theta)*sin(_phi);
	double a32 = -1*sin(_theta)*cos(_phi);
	double a33 = cos(_theta);

	dblMat rotMatrix(3,3,0.0);
	rotMatrix[0][0] = a11;
	rotMatrix[1][0] = a12;
	rotMatrix[2][0] = a13;
	rotMatrix[0][1] = a21;
	rotMatrix[1][1] = a22;
	rotMatrix[2][1] = a23;
	rotMatrix[0][2] = a31;
	rotMatrix[1][2] = a32;
	rotMatrix[2][2] = a33;

	for (UInt i = 0; i < itsChains.size(); i ++)
	{
		transform(i,rotMatrix);
	}
	return;
}

void protein::undoEulerRotate(UInt _chain, const double _phi, const double _theta, const double _psi)
{
	// calculate inverse rotation matrix (transpose of rotation matrix)
	double a11 = cos(_psi)*cos(_phi) - cos(_theta)*sin(_phi)*sin(_psi);
	double a12 = -1*sin(_psi)*cos(_phi) - cos(_theta)*sin(_phi)*cos(_psi);
	double a13 = sin(_theta)*sin(_phi);
	double a21 = cos(_psi)*sin(_phi) + cos(_theta)*cos(_phi)*sin(_psi);
	double a22 = -1*sin(_psi)*sin(_phi) + cos(_theta)*cos(_phi)*cos(_psi);
	double a23 = -1*sin(_theta)*cos(_phi);
	double a31 = sin(_theta)*sin(_psi);
	double a32 = sin(_theta)*cos(_psi);
	double a33 = cos(_theta);

	dblMat rotMatrix(3,3,0.0);
	rotMatrix[0][0] = a11;
	rotMatrix[1][0] = a12;
	rotMatrix[2][0] = a13;
	rotMatrix[0][1] = a21;
	rotMatrix[1][1] = a22;
	rotMatrix[2][1] = a23;
	rotMatrix[0][2] = a31;
	rotMatrix[1][2] = a32;
	rotMatrix[2][2] = a33;

	transform(_chain,rotMatrix);
	return;
}

void protein::undoEulerRotate(const double _phi, const double _theta, const double _psi)
{
	// calculate inverse rotation matrix (transpose of rotation matrix)
	double a11 = cos(_psi)*cos(_phi) - cos(_theta)*sin(_phi)*sin(_psi);
	double a12 = -1*sin(_psi)*cos(_phi) - cos(_theta)*sin(_phi)*cos(_psi);
	double a13 = sin(_theta)*sin(_phi);
	double a21 = cos(_psi)*sin(_phi) + cos(_theta)*cos(_phi)*sin(_psi);
	double a22 = -1*sin(_psi)*sin(_phi) + cos(_theta)*cos(_phi)*cos(_psi);
	double a23 = -1*sin(_theta)*cos(_phi);
	double a31 = sin(_theta)*sin(_psi);
	double a32 = sin(_theta)*cos(_psi);
	double a33 = cos(_theta);

	dblMat rotMatrix(3,3,0.0);
	rotMatrix[0][0] = a11;
	rotMatrix[1][0] = a12;
	rotMatrix[2][0] = a13;
	rotMatrix[0][1] = a21;
	rotMatrix[1][1] = a22;
	rotMatrix[2][1] = a23;
	rotMatrix[0][2] = a31;
	rotMatrix[1][2] = a32;
	rotMatrix[2][2] = a33;

	for (UInt i = 0; i < itsChains.size(); i ++)
	{
		transform(i,rotMatrix);
	}
	return;
}

void protein::translate(const UInt _index, const double _x,const double _y,const double _z)
{	dblVec _vec;
	_vec.newsize(3);
	_vec[0] = _x;
	_vec[1] = _y;
	_vec[2] = _z;
	translate(_index, _vec);
}

void protein::translate(const double _x, const double _y, const double _z)
{
	dblVec _vec;
	_vec.newsize(3);
	_vec[0] = _x;
	_vec[1] = _y;
	_vec[2] = _z;
	for (UInt i=0; i< itsChains.size(); i++)
	{	translate(i,_vec);
	}
}

void protein::transform(const UInt _index, const dblMat& _dblMat)
{   if(_index < itsChains.size())
    {   itsChains[_index]->transform(_dblMat);
    }
}

void protein::rotate(const UInt _index, const axis _axis, const double _theta)
{	point origin;
	// The default is to set this point to the origin
	//cout << "ROTATING CHAIN " << _index << endl;
	origin.setCoords(0.0,0.0,0.0);
	dblVec vec = dblVec(3);
	for (int i = 0; i<vec.dim(); i++)
	{   vec[i] = 0.0;
	}
	if (_axis == X_axis)
	{   vec[0] = 1.0;
	}
	else if (_axis == Y_axis)
	{   vec[1] = 1.0;
	}
	else if (_axis == Z_axis)
	{   vec[2]  = 1.0;
	}
	rotate(_index, origin , vec, _theta);
}

void protein::rotate(const axis _axis, const double _theta)
{	point origin;
	// The default is to set this point to the origin
	origin.setCoords(0.0,0.0,0.0);
	dblVec vec = dblVec(3);
	for (int i = 0; i<vec.dim(); i++)
	{   vec[i] = 0.0;
	}
	if (_axis == X_axis)
	{   vec[0] = 1.0;
	}
	else if (_axis == Y_axis)
	{   vec[1] = 1.0;
	}
	else if (_axis == Z_axis)
	{   vec[2]  = 1.0;
	}
	for (UInt i=0; i<itsChains.size(); i++)
	{
		rotate(i, origin , vec, _theta);
	}
}

void protein::rotate(const UInt _index, const point& _point,const dblVec& _R_axis, const double _theta)
{   if( _index < itsChains.size())
    {   itsChains[_index]->rotate(_point, _R_axis, _theta);
    }
}


double protein::getSelfEnergy(UInt _chainIndex, UInt _resIndex)
{
	if (calcSelfEnergy == false && messagesActive) cout << "WARNING:  protein::getSelfEnergy(...) invoked even though calcSelfENergy set to false ... " << endl;
	if (_chainIndex >= 0 && _chainIndex < itsChains.size())
	{
		double selfEnergy;
		selfEnergy = itsChains[_chainIndex]->getSelfEnergy(_resIndex);
		// and add self energy of symmetry linked chains
		UInt independentChainIndex = 0;
		for (UInt i = 0; i < itsIndependentChainsMap.size(); i++)
		{
			if (itsIndependentChainsMap[i] == _chainIndex) independentChainIndex = i;
		}
		UInt numSymmetryLinkedChains = itsChainLinkageMap[independentChainIndex].size();
		if (itsChainLinkageMap[independentChainIndex][0] != -1)
		{
			for (UInt i = 0; i < numSymmetryLinkedChains; i ++)
			{
				selfEnergy += itsChains[itsChainLinkageMap[independentChainIndex][i]]->getSelfEnergy(_resIndex);
			}
		}
		return selfEnergy;
	}
	else
	{
		cout << "Invalid Chain Specifier: " << _chainIndex << endl;
		exit(1);
	}
	return -1;
}

void protein::coilcoil(const double _pitch)
{
    if (_pitch == 0.0)
    {
        cout << "ERROR in protein::coilcoil(...)  pitch cannot be zero!!" << endl;
        cout << "\t protein unchanged." << endl;
        return;
    }
    for (UInt i = 0; i < itsChains.size(); i ++)
    {
        itsChains[i]->coilcoil(_pitch);
    }
    return;
}

void protein::coilcoil(UInt _chain, double _pitch)
{
	if (_pitch == 0.0)
	{
		cout << "ERROR in protein::coilcoil(...)  pitch cannot be zero!!" << endl;
		cout << "\t protein unchanged." << endl;
		return;
	}
	if (_chain >=0 && _chain < itsChains.size())
	{
		itsChains[_chain]->coilcoil(_pitch);
	}
	return;
}

// surface area and solvation energy functions.
void protein::initializeSpherePoints()
{
	for (UInt i = 0; i < itsChains.size(); i ++)
	{
		itsChains[i]->initializeSpherePoints();
	}
	return;
}

void protein::initializeSpherePoints(UInt _chain)
{
	if (_chain >=0 && _chain < itsChains.size() )
	{
		itsChains[_chain]->initializeSpherePoints();
	}
	else
	{
		cout << "ERROR in initializeSpherePoints ... chain index out of range" << endl;
	}
	return;
}

void protein::initializeSpherePoints(UInt _chain, UInt _residue)
{
	if (_chain >= 0 && _chain < itsChains.size())
	{
		itsChains[_chain]->initializeSpherePoints(_residue);
	}
	else
	{
		cout << "ERROR in initializeSpherePoints ... chain index out of range" << endl;
	}
	return;
}

double protein::tabulateSurfaceArea()
{
	double surfaceArea = 0.0;
	for (UInt i = 0; i < itsChains.size(); i++)
	{
		surfaceArea += itsChains[i]->tabulateSurfaceArea();
	}

	return surfaceArea;
}

double protein::tabulateSurfaceArea(UInt _chain)
{
	double surfaceArea = 0.0;
	if (_chain >= 0 && _chain < itsChains.size() )
	{
		surfaceArea = itsChains[_chain]->tabulateSurfaceArea();
	}
	else
	{
		cout << "ERROR in tabulateSurfaceArea ... chain index out of range." << endl;
	}

	return surfaceArea;
}

double protein::tabulateSurfaceArea(UInt _chain, UInt _residue)
{
	double surfaceArea = 0.0;
	if (_chain >=0 && _chain < itsChains.size() )
	{
		surfaceArea = itsChains[_chain]->tabulateSurfaceArea(_residue);
	}
	else
	{
		cout << "ERROR in tabulateSurfaceArea ... chain index out of range." << endl;
	}

	return surfaceArea;
}

double protein::tabulateSurfaceArea(UInt _chainIndex, UInt _residueIndex, UInt _atomIndex)
{
	double surfaceArea = 0.0;
	if (_chainIndex >=0 && _chainIndex < itsChains.size() )
	{
		surfaceArea = itsChains[_chainIndex]->tabulateSurfaceArea(_residueIndex, _atomIndex);
	}
	else
	{
		cout << "ERROR in tabulateSurfaceArea ... chain index out of range." << endl;
	}

	return surfaceArea;
}

double protein::getItsSolvationParam()
{
	return itsSolvationParam;
}

void protein::setItsSolvationParam(UInt _param)
{
	itsSolvationParam = _param;
}

void protein::removeSpherePoints()
{
	for (UInt i = 0; i < itsChains.size(); i++ )
	{
		itsChains[i]->removeIntraChainSpherePoints();
		for (UInt j = 0; j < itsChains.size(); j ++)
		{
			if ( i != j ) itsChains[i]->removeInterChainSpherePoints(itsChains[j]);
		}
	}
	return;
}

void protein::removeSpherePoints(UInt _chain)
{
	if  (_chain >= 0 && _chain < itsChains.size() )
	{
		itsChains[_chain]->removeIntraChainSpherePoints();
		for (UInt j = 0; j < itsChains.size(); j ++)
		{
			if (_chain != j) itsChains[_chain]->removeInterChainSpherePoints(itsChains[j]);
		}
	}
	else
	{
		cout << "ERROR in removeSpherePoints ... chain index out of range." << endl;
	}
	return;
}

void protein::removeSpherePoints(UInt _chain, UInt _residue)
{
	if (_chain >=0 && _chain < itsChains.size())
	{
		itsChains[_chain]->removeIntraChainSpherePoints(_residue);
		for (UInt j =0; j < itsChains.size(); j ++)
		{
			if (_chain != j) itsChains[_chain]->removeInterChainSpherePoints(_residue, itsChains[j]);
		}
	}
	else
	{
		cout << "ERROR in removeSpherePoints ... chain index out of range." << endl;
	}
	return;
}

residue* protein::superimposeGLY(const UInt _chain, const UInt _residue)
{
	// get coordinates of a glycine complete with hydrogens superimposed onto the mainchain of this residue
	if (_chain >=0 && _chain < itsChains.size() )
	{
		return itsChains[_chain]->superimposeGLY(_residue);
	}
	else
	{
		cout << "ERROR in protein::superimposeGLY ... chain index out of range." << endl;
		exit(1);
	}
}

double protein::calculateHCA_O_hBondEnergy()
{
	double hbondEnergy = 0.0;
	for (UInt i = 0; i < itsChains.size(); i ++)
	{
		for (UInt j = 0; j < itsChains.size(); j ++)
		{
			if (i!=j)hbondEnergy += itsChains[i]->calculateHCA_O_hBondEnergy(itsChains[j]);
		}
	}
	return hbondEnergy;
}

dblVec protein::getBackBoneCentroid()
{
	dblVec centroid(3);
	centroid[0] = 0.0; centroid[1] = 0.0; centroid[2] = 0.0;
	for (UInt i = 0; i < itsChains.size(); i ++)
	{
		centroid = centroid + itsChains[i]->getBackBoneCentroid();
	}

	centroid = centroid / (double)itsChains.size();
	return centroid;
}

typedef UIntVec::iterator iterUIntVec;

double protein::getResPairEnergy(const UInt _chain1, const UInt _res1, const UInt _chain2, const UInt _res2)
{
	//if (messagesActive) cout << "pair energy " << _chain1 << " " << _res1 << " " << _chain2 << " " << _res2 << ":  " ;
	if (_chain1 >= 0 && _chain1 < itsChains.size() && _chain2 >= 0 && _chain2 < itsChains.size())
	{
		double energy = itsChains[_chain1]->getInterEnergy(_res1, itsChains[_chain2], _res2);
		//cout << _chain1 << " " << _res1 << " " << _chain2 << " " << _res2 << " " <<energy << endl;
		return energy;
	}
	else
	{
		cout << _chain1 << " " << _res1 << " " << _chain2 << " " << _res2 << endl;
		cout << "ERROR:  chain indices out of range in getResPairEnergy" << endl;
		exit(1);
	}
}

vector <dblVec> protein::saveCoords( UInt chainIndex, UInt resIndex)
{
	UInt nAtoms = getNumAtoms(chainIndex, resIndex);
	vector <dblVec> allCoords;
	for (UInt i=0; i<nAtoms; i++)
	{
		dblVec coords = getCoords(chainIndex, resIndex, i);
		allCoords.push_back(coords);
	}
	return allCoords;
}

void protein::setAllCoords( UInt chainIndex, UInt resIndex, vector<dblVec> allCoords)
{
	UInt nAtoms = getNumAtoms(chainIndex, resIndex);
	for (UInt i=0; i<nAtoms; i++)
	{
		setCoords(chainIndex, resIndex, i, allCoords[i]);
	}
}

void protein::protOpt(bool _backboneRelaxation)
{
	protRelax(1000);
	protOpt(_backboneRelaxation);
}

void protein::protOpt(bool _backboneRelaxation, UIntVec _frozenResidues, UIntVec _activeChains)
{
	protRelax(_frozenResidues, _activeChains);
	protOpt(_backboneRelaxation, _frozenResidues, _activeChains);
}

void protein::protMin(bool _backbone)
{
	// Sidechain and backslide optimization with a local dielectric scaling of electrostatics and corresponding Born/Gill implicit solvation energy
	
	//--Initialize variables for loop, calculate starting energy and build energy vectors-----
	saveCurrentState();
	UInt randchain, randres, resnum, backboneOrSidechain = 1;
	UInt clashes, clashesStart, bbClashes, bbClashesStart, chainNum = getNumChains(), plateau = 1000;
	double Energy, pastEnergy = protEnergy(), deltaEnergy, sPhi, sPsi, nobetter = 0.0, KT = KB*Temperature();
	double rotX, rotY, rotZ, transX, transY, transZ;
	vector < vector <double>> currentSidechainConf, newSidechainConf; srand (time(NULL)); vector <double> backboneAngles(2);
	bool sidechainTest, backboneTest, cofactorTest, revert, energyTest, boltzmannAcceptance;
	vector <dblVec> currentCoords;
	//--Run optimizaiton loop to local minima defined by an RT plateau------------------------
	do{
		//--choose random residue and set variables
		randchain = rand() % chainNum, resnum = getNumResidues(randchain), randres = rand() % resnum;
		clashesStart = getNumHardClashes(); nobetter++;
		backboneTest = false, sidechainTest = false, cofactorTest = false, energyTest = false, revert = true;
		if (isCofactor(randchain, randres))
		{
			//--Rock and Roll cofactor in site
			cofactorTest = true;
			currentCoords.clear();
			currentCoords = saveCoords(randchain, randres);
			rotX = rand() % 2, rotY = rand() % 2, rotZ = rand() % 2;
			transX = (rand() % 30)/100, transY = (rand() % 30)/100, transZ = (rand() % 30)/100;
			rotateChainRelative(randchain,X_axis,rotX), rotateChainRelative(randchain,Y_axis,rotY), rotateChainRelative(randchain,Z_axis,rotZ);
			translateChain(randchain, transX, transY, transZ);
			clashes = getNumHardClashes();
			if (clashes <= clashesStart){
					energyTest = true; revert = false;
			}
		}
		else{
			//--Backbone conformation trial--------------------------------------------------------
			if (_backbone) {backboneOrSidechain = rand() % 2;}
			if (randres > 0 && randres < resnum-2 && backboneOrSidechain == 0){
				backboneTest = true; bbClashesStart = getNumHardBackboneClashes();
				sPhi = getPhi(randchain,randres), sPsi = getPsi(randchain,randres);
				backboneAngles = getRandConformationFromBackboneType(sPhi, sPsi);
				setDihedral(randchain,randres,backboneAngles[0],0,0); setDihedral(randchain,randres,backboneAngles[1],1,0);
				bbClashes = getNumHardBackboneClashes();
				if (bbClashes <= bbClashesStart){
					clashes = getNumHardClashes();
					if (clashes <= clashesStart){
						energyTest = true; revert = false;
					}
				}
			}
			//--Sidechain conformation trial--------------------------------------------------------
			else{
				sidechainTest = true;
				currentSidechainConf = getSidechainDihedrals(randchain, randres);
				newSidechainConf = randContinuousSidechainConformation(randchain, randres);
				setSidechainDihedralAngles(randchain, randres, newSidechainConf);
				clashes = getNumHardClashes();
				if (clashes <= clashesStart){
					energyTest = true; revert = false;
				}
			}
		}
		//--Energy-Test-------------------------------------------------------------------------
		if (energyTest){
			Energy = protEnergy();
			deltaEnergy = Energy - pastEnergy;
			boltzmannAcceptance = boltzmannEnergyCriteria(deltaEnergy);
			if (boltzmannAcceptance){
				pastEnergy = Energy;
				if (deltaEnergy < -KT){nobetter = 0;}
			}
			else{revert = true;}
		}
		//--Revert conformation-----------------------------------------------------------------
		if (revert){
			if(cofactorTest)
			{
				setAllCoords(randchain, randres, currentCoords);
			}
			if(backboneTest){
				setDihedral(randchain,randres,sPhi,0,0);
				setDihedral(randchain,randres,sPsi,1,0);
			}
			if(sidechainTest){
				setSidechainDihedralAngles(randchain, randres, currentSidechainConf);
			}
		}
	} while (nobetter < plateau);
	return;
}

void protein::protMin(bool _backbone, UIntVec _frozenResidues, UIntVec _activeChains)
{
	// Sidechain and backslide optimization with a local dielectric scaling of electrostatics and corresponding Born/Gill implicit solvation energy
	
	//--Initialize variables for loop, calculate starting energy and build energy vectors-----
	saveCurrentState();
	UInt randchain, randres, resnum, backboneOrSidechain = 1;
	UInt clashes, clashesStart, bbClashes, bbClashesStart, chainNum = _activeChains.size(), plateau = 1000;
	double Energy, pastEnergy = protEnergy(), deltaEnergy, sPhi, sPsi,nobetter = 0.0, KT = KB*Temperature();
	double rotX, rotY, rotZ, transX, transY, transZ;
	vector < vector <double>> currentSidechainConf, newSidechainConf; srand (time(NULL)); vector <double> backboneAngles(2);
	bool sidechainTest, backboneTest, revert, cofactorTest, energyTest, skip, boltzmannAcceptance;
	vector <dblVec> currentCoords;
	//--Run optimizaiton loop to local minima defined by an RT plateau------------------------
	do{
		//--choose random residue not frozen of active chains
			randchain = _activeChains[rand() % chainNum];
			do{
				skip = false;
				randres = rand() % getNumResidues(randchain);
				for (UInt i = 0; i < _frozenResidues.size(); i++)
				{
					if (randres == _frozenResidues[i]) {skip = true; break;}
				}
			} while (skip);
		clashesStart = getNumHardClashes(); resnum = getNumResidues(randchain); nobetter++;
		backboneTest = false, sidechainTest = false, cofactorTest = false, energyTest = false, revert = true;
		
		if (isCofactor(randchain, randres))
		{
			//--Rock and Roll cofactor in site
			cofactorTest = true;
			currentCoords.clear();
			currentCoords = saveCoords(randchain, randres);
			rotX = rand() % 2, rotY = rand() % 2, rotZ = rand() % 2;
			transX = (rand() % 30)/100, transY = (rand() % 30)/100, transZ = (rand() % 30)/100;
			rotateChainRelative(randchain,X_axis,rotX), rotateChainRelative(randchain,Y_axis,rotY), rotateChainRelative(randchain,Z_axis,rotZ);
			translateChain(randchain, transX, transY, transZ);
			clashes = getNumHardClashes();
			if (clashes <= clashesStart){
					energyTest = true; revert = false;
			}
		}
		else{
			//--Backbone conformation trial--------------------------------------------------------
			if (_backbone) {backboneOrSidechain = rand() % 2;}
			if (randres > 0 && randres < resnum-2 && backboneOrSidechain == 0){
				backboneTest = true; bbClashesStart = getNumHardBackboneClashes();
				sPhi = getPhi(randchain,randres), sPsi = getPsi(randchain,randres);
				backboneAngles = getRandConformationFromBackboneType(sPhi, sPsi);
				setDihedral(randchain,randres,backboneAngles[0],0,0); setDihedral(randchain,randres,backboneAngles[1],1,0);
				bbClashes = getNumHardBackboneClashes();
				if (bbClashes <= bbClashesStart){
					clashes = getNumHardClashes();
					if (clashes <= clashesStart){
						energyTest = true; revert = false;
					}
				}
			}
			//--Sidechain conformation trial--------------------------------------------------------
			else{
				sidechainTest = true;
				currentSidechainConf = getSidechainDihedrals(randchain, randres);
				newSidechainConf = randContinuousSidechainConformation(randchain, randres);
				setSidechainDihedralAngles(randchain, randres, newSidechainConf);
				clashes = getNumHardClashes();
				if (clashes <= clashesStart){
					energyTest = true; revert = false;
				}
			}
		}
		//--Energy-Test-------------------------------------------------------------------------
		if (energyTest){
			Energy = protEnergy();
			deltaEnergy = Energy - pastEnergy;
			boltzmannAcceptance = boltzmannEnergyCriteria(deltaEnergy);
			if (boltzmannAcceptance){
				pastEnergy = Energy;
				if (deltaEnergy < -KT){nobetter = 0;}
			}
			else{revert = true;}
		}
		//--Revert conformation-----------------------------------------------------------------
		if (revert){
			if(cofactorTest)
			{
				setAllCoords(randchain, randres, currentCoords);
			}
			if(backboneTest){
				setDihedral(randchain,randres,sPhi,0,0);
				setDihedral(randchain,randres,sPsi,1,0);
			}
			if(sidechainTest){
				setSidechainDihedralAngles(randchain, randres, currentSidechainConf);
			}
		}
	} while (nobetter < plateau);
	return;
}

void protein::protMin(bool _backbone, UInt chainIndex, UInt resIndex)
{
	// Sidechain and backslide optimization with a local dielectric scaling of electrostatics and corresponding Born/Gill implicit solvation energy
	
	//--Initialize variables for loop, calculate starting energy and build energy vectors-----
	saveCurrentState();
	UInt randchain = chainIndex, randres = resIndex, resnum, backboneOrSidechain = 1;
	UInt clashes, clashesStart, bbClashes, bbClashesStart, plateau = 1000;
	double Energy, pastEnergy = protEnergy(), deltaEnergy, sPhi, sPsi,nobetter = 0.0, KT = KB*Temperature();
	double rotX, rotY, rotZ, transX, transY, transZ;
	vector < vector <double>> currentSidechainConf, newSidechainConf; srand (time(NULL)); vector <double> backboneAngles(2);
	bool sidechainTest, backboneTest, revert, cofactorTest, energyTest, boltzmannAcceptance;
	vector <dblVec> currentCoords;
	//--Run optimizaiton loop to local minima defined by an RT plateau------------------------
	do{
		clashesStart = getNumHardClashes(); resnum = getNumResidues(randchain); nobetter++;
		backboneTest = false, sidechainTest = false,  cofactorTest = false, energyTest = false, revert = true;
		
		if (isCofactor(randchain, randres))
		{
			//--Rock and Roll cofactor in site
			cofactorTest = true;
			currentCoords.clear();
			currentCoords = saveCoords(randchain, randres);
			rotX = rand() % 2, rotY = rand() % 2, rotZ = rand() % 2;
			transX = (rand() % 30)/100, transY = (rand() % 30)/100, transZ = (rand() % 30)/100;
			rotateChainRelative(randchain,X_axis,rotX), rotateChainRelative(randchain,Y_axis,rotY), rotateChainRelative(randchain,Z_axis,rotZ);
			translateChain(randchain, transX, transY, transZ);
			clashes = getNumHardClashes();
			if (clashes <= clashesStart){
					energyTest = true; revert = false;
			}
		}
		else{
			//--Backbone conformation trial--------------------------------------------------------
			if (_backbone) {backboneOrSidechain = rand() % 2;}
			if (randres > 0 && randres < resnum-2 && backboneOrSidechain == 0){
				backboneTest = true; bbClashesStart = getNumHardBackboneClashes();
				sPhi = getPhi(randchain,randres), sPsi = getPsi(randchain,randres);
				backboneAngles = getRandConformationFromBackboneType(sPhi, sPsi);
				setDihedral(randchain,randres,backboneAngles[0],0,0); setDihedral(randchain,randres,backboneAngles[1],1,0);
				bbClashes = getNumHardBackboneClashes();
				if (bbClashes <= bbClashesStart){
					clashes = getNumHardClashes();
					if (clashes <= clashesStart){
						energyTest = true; revert = false;
					}
				}
			}
			//--Sidechain conformation trial--------------------------------------------------------
			else{
				sidechainTest = true;
				currentSidechainConf = getSidechainDihedrals(randchain, randres);
				newSidechainConf = randContinuousSidechainConformation(randchain, randres);
				setSidechainDihedralAngles(randchain, randres, newSidechainConf);
				clashes = getNumHardClashes();
				if (clashes <= clashesStart){
					energyTest = true; revert = false;
				}
			}
		}
		//--Energy-Test-------------------------------------------------------------------------
		if (energyTest){
			Energy = protEnergy();
			deltaEnergy = Energy - pastEnergy;
			boltzmannAcceptance = boltzmannEnergyCriteria(deltaEnergy);
			if (boltzmannAcceptance){
				pastEnergy = Energy;
				if (deltaEnergy < -KT){nobetter = 0;}
			}
			else{revert = true;}
		}
		//--Revert conformation-----------------------------------------------------------------
		if (revert){
			if(cofactorTest)
			{
				setAllCoords(randchain, randres, currentCoords);
			}
			if(backboneTest){
				setDihedral(randchain,randres,sPhi,0,0);
				setDihedral(randchain,randres,sPsi,1,0);
			}
			if(sidechainTest){
				setSidechainDihedralAngles(randchain, randres, currentSidechainConf);
			}
		}
	} while (nobetter < plateau);
	return;
}

void protein::protRelax(UInt _plateau)
{   // Sidechain and backrub optimization with a local dielectric scaling of electrostatics and corresponding Born/Gill implicit solvation energy
	//_plateau: the number of consecutive optimization cycles without an energy decrease (default: 150 for general purpose optimization)
	
	saveCurrentState();
	UInt pastProtClashes = getNumHardClashes();
	if (pastProtClashes > 0)
	{	
		//--Initialize variables for loop, calculate starting energy and build energy vectors---------------
		UInt randchain, randres, randrestype, resnum, randrot, chainNum = getNumChains(), protClashes, resClashes, medResC, nobetter = 0;
		vector < vector <double> > currentRot; vector <UIntVec> allowedRots; srand (time(NULL));

		//--Run optimizaiton loop to relative minima, determined by _plateau----------------------------
		do
		{   //--choose random residue
			randchain = rand() % chainNum;
			resnum = getNumResidues(randchain);
			randres = rand() % resnum;
			randrestype = getTypeFromResNum(randchain, randres);
			nobetter++;
	
			//--Rotamer optimization-----------------------------------------------------------------------
			medResC = getMedianResidueNumHardClashes();
			resClashes = getNumHardClashes(randchain, randres);
			if (resClashes > medResC)
			{
				currentRot = getSidechainDihedrals(randchain, randres);
				allowedRots = getAllowedRotamers(randchain, randres, randrestype);
	
				//--Try a max of one rotamer per branchpoint and keep if an improvement, else revert
				for (UInt b = 0; b < residue::getNumBpt(randrestype); b++)
				{
					if (allowedRots[b].size() > 0)
					{
						randrot = rand() % allowedRots[b].size();
						setRotamerWBC(randchain, randres, b, allowedRots[b][randrot]);
						protClashes = getNumHardClashes();
						if (protClashes < pastProtClashes)
						{
							nobetter = 0, pastProtClashes = protClashes; break;
						}
						else
						{
							setSidechainDihedralAngles(randchain, randres, currentRot);
						}
					}
				}
			}
		} while (nobetter < _plateau);
	}
	return;
}

void protein::protRelax(UIntVec _frozenResidues, UIntVec _activeChains)
{   // Sidechain and backrub optimization with a local dielectric scaling of electrostatics and corresponding Born/Gill implicit solvation energy
	//_plateau: the number of consecutive optimization cycles without an energy decrease (default: 150 for general purpose optimization)
	saveCurrentState();
	UInt pastProtClashes = getNumHardClashes();
	if (pastProtClashes > 0)
	{	
		//--Initialize variables for loop, calculate starting energy and build energy vectors---------------
		UInt randchain, randres, randrestype, randrot, chainNum = _activeChains.size(), protClashes, resClashes, medResC, _plateau = 1000, nobetter = 0;
		vector < vector <double> > currentRot; vector <UIntVec> allowedRots; srand (time(NULL));
		bool skip;
		//--Run optimizaiton loop to relative minima, determined by _plateau----------------------------
		do
		{
			//--choose random residue not frozen of active chains
			randchain = _activeChains[rand() % chainNum];
			do{
				skip = false;
				randres = rand() % getNumResidues(randchain);
				for (UInt i = 0; i < _frozenResidues.size(); i++)
				{
					if (randres == _frozenResidues[i]) {skip = true; break;}
				}
			} while (skip);
			randrestype = getTypeFromResNum(randchain, randres);
			nobetter++;
	
			//--Rotamer optimization-----------------------------------------------------------------------
			medResC = getMedianResidueNumHardClashes();
			resClashes = getNumHardClashes(randchain, randres);
			if (resClashes > medResC)
			{
				currentRot = getSidechainDihedrals(randchain, randres);
				allowedRots = getAllowedRotamers(randchain, randres, randrestype);
	
				//--Try a max of one rotamer per branchpoint and keep if an improvement, else revert
				for (UInt b = 0; b < residue::getNumBpt(randrestype); b++)
				{
					if (allowedRots[b].size() > 0)
					{
						randrot = rand() % allowedRots[b].size();
						setRotamerWBC(randchain, randres, b, allowedRots[b][randrot]);
						protClashes = getNumHardClashes();
						if (protClashes < pastProtClashes)
						{
							nobetter = 0, pastProtClashes = protClashes; break;
						}
						else
						{
							setSidechainDihedralAngles(randchain, randres, currentRot);
						}
					}
				}
			}
		} while (nobetter < _plateau);
	}
	return;
}

void protein::protSampling(UInt iterations)
{
	//--Initialize variables for loop, calculate starting energy and build energy vectors-----
	saveCurrentState();
	UInt randchain, randres, resnum, changes = 0, backboneOrSidechain = 1;
	UInt clashes, clashesStart, bbClashes, bbClashesStart, chainNum = getNumChains();
	double Energy, pastEnergy = protEnergy(), deltaEnergy, sPhi, sPsi;
	vector < vector <double>> currentSidechainConf, newSidechainConf; srand (time(NULL)); vector <double> backboneAngles(2);
	bool sidechainTest, backboneTest, revert, energyTest, boltzmannAcceptance;
	
	//--Run sampling loop to number of iterations------------------------
	do
	{
		//--choose random residue and set variables
		randchain = rand() % chainNum, resnum = getNumResidues(randchain), randres = rand() % resnum;
		clashesStart = getNumHardClashes(); backboneOrSidechain = rand() % 2;
		backboneTest = false, sidechainTest = false, energyTest = false, revert = true;
		
		//--Backbone conformation trial--------------------------------------------------------
		if (randres > 0 && randres < resnum-2 && backboneOrSidechain == 0){
			backboneTest = true; bbClashesStart = getNumHardBackboneClashes();
			sPhi = getPhi(randchain,randres), sPsi = getPsi(randchain,randres);
			backboneAngles = getRandConformationFromBackboneType(sPhi, sPsi);
			setDihedral(randchain,randres,backboneAngles[0],0,0); setDihedral(randchain,randres,backboneAngles[1],1,0);
			bbClashes = getNumHardBackboneClashes();
			if (bbClashes <= bbClashesStart){
				clashes = getNumHardClashes();
				if (clashes <= clashesStart){
					energyTest = true; revert = false;
				}
			}
		}
		//--Sidechain conformation trial--------------------------------------------------------
		else{
			sidechainTest = true;
			currentSidechainConf = getSidechainDihedrals(randchain, randres);
			newSidechainConf = randContinuousSidechainConformation(randchain, randres);
			setSidechainDihedralAngles(randchain, randres, newSidechainConf);
			clashes = getNumHardClashes();
			if (clashes <= clashesStart){
				energyTest = true; revert = false;
			}
		}
		//--Energy-Test-------------------------------------------------------------------------
		if (energyTest){
			Energy = protEnergy();
			deltaEnergy = Energy - pastEnergy;
			boltzmannAcceptance = boltzmannEnergyCriteria(deltaEnergy);
			if (boltzmannAcceptance){
				pastEnergy = Energy; changes++;
				cout << changes << " " << Energy << endl;
			}
			else{revert = true;}
		}
		//--Revert conformation-----------------------------------------------------------------
		if (revert){
			if(backboneTest){
				setDihedral(randchain,randres,sPhi,0,0);
				setDihedral(randchain,randres,sPsi,1,0);
			}
			if(sidechainTest){
				setSidechainDihedralAngles(randchain, randres, currentSidechainConf);
			}
		}
	} while(changes < iterations);
	return;
}

double protein::getResiduesPerTurn(double phi, double psi)
{
	double phipsisum = phi+psi;
	double handedness;
	if ((phipsisum > 0 && phipsisum < 180)|| phipsisum < -180){handedness = -1.0;}
	else{handedness = 1.0;}
	double angleSumHalfRad = ((phipsisum)/2)*PI/180;
	double radAngle = acos(-0.3333333-0.6666666*cos(2*angleSumHalfRad));
	double radAngletoDeg = radAngle*180/PI;
	double residuesPerTurn = (360/radAngletoDeg)*handedness;
	return residuesPerTurn;
}

UInt protein::getBackboneSequenceType(double RPT, double phi)
{
	if (RPT <= -4.48 && phi <= 0)					{return 0;}
	if (RPT > -4.48  && RPT <= -3.86 && phi <= 0)	{return 1;}
	if (RPT > -3.86  && RPT <= -3.24 && phi <= 0)	{return 2;}
	if (RPT > -3.24  && RPT <= -2.62 && phi <= 0)	{return 3;}
	if (RPT > -2.62  && RPT <= -2.00 && phi <= 0)	{return 4;}
	if (RPT >  2.00  && RPT <=  2.62 && phi <= 0)	{return 5;}
	if (RPT >  2.62  && RPT <=  3.24 && phi <= 0)	{return 6;}
	if (RPT >  3.24  && RPT <=  3.86 && phi <= 0)	{return 7;}
	if (RPT >  3.86  && RPT <=  4.48 && phi <= 0)	{return 8;}
	if (RPT >  4.48 && phi <= 0)					{return 9;}
	if (RPT <= -4.48 && phi > 0)					{return 10;}
	if (RPT > -4.48  && RPT <= -3.86 && phi > 0)	{return 11;}
	if (RPT > -3.86  && RPT <= -3.24 && phi > 0)	{return 12;}
	if (RPT > -3.24  && RPT <= -2.62 && phi > 0)	{return 13;}
	if (RPT > -2.62  && RPT <= -2.00 && phi > 0)	{return 14;}
	if (RPT >  2.00  && RPT <=  2.62 && phi > 0)	{return 15;}
	if (RPT >  2.62  && RPT <=  3.24 && phi > 0)	{return 16;}
	if (RPT >  3.24  && RPT <=  3.86 && phi > 0)	{return 17;}
	if (RPT >  3.86  && RPT <=  4.48 && phi > 0)	{return 18;}
	if (RPT >  4.48 && phi > 0)						{return 19;}
	return 0;
}

vector <double> protein::getRandPhiPsifromBackboneSequenceType(UInt _RPTType)
{
	vector <double> angles(2);
	int b, psi, phi;

	if (_RPTType == 0){
		do{
			b = 124 + (rand() % 57);
			phi = ((rand() % 180)+1)-180;
			psi = -1 * phi + b;
		}while(psi > 180 || psi < -180);
		angles[0] = phi; angles[1] = psi;
		return angles;
	}
	if (_RPTType == 1){
		do{
			b = 95 + (rand() % 29);
			phi = ((rand() % 180)+1)-180;
			psi = -1 * phi + b;
		}while(psi > 180 || psi < -180);
		angles[0] = phi; angles[1] = psi;
		return angles;
	}
	if (_RPTType == 2){
		do{
			b = 58 + (rand() % 37);
			phi = ((rand() % 180)+1)-180;
			psi = -1 * phi + b;
		}while(psi > 180 || psi < -180);
		angles[0] = phi; angles[1] = psi;
		return angles;
	}
	if (_RPTType == 3){
		do{
			b = (rand() % 58);
			phi = ((rand() % 180)+1)-180;
			psi = -1 * phi + b;
		}while(psi > 180 || psi < -180);
		angles[0] = phi; angles[1] = psi;
		return angles;
	}
	if (_RPTType == 4){
		do{
			b = -58 + (rand() % 58);
			phi = ((rand() % 180)+1)-180;
			psi = -1 * phi + b;
		}while(psi > 180 || psi < -180);
		angles[0] = phi; angles[1] = psi;
		return angles;
	}
	if (_RPTType == 5){
		do{
			b = -95 + (rand() % 37);
			phi = ((rand() % 180)+1)-180;
			psi = -1 * phi + b;
		}while(psi > 180 || psi < -180);
		angles[0] = phi; angles[1] = psi;
		return angles;
	}
	if (_RPTType == 6){
		do{
			b = -124 + (rand() % 29);
			phi = ((rand() % 180)+1)-180;
			psi = -1 * phi + b;
		}while(psi > 180 || psi < -180);
		angles[0] = phi; angles[1] = psi;
		return angles;
	}
	if (_RPTType == 7){
		do{
			b = -180 + (rand() % 56);
			phi = ((rand() % 180)+1)-180;
			psi = -1 * phi + b;
		}while(psi > 180 || psi < -180);
		angles[0] = phi; angles[1] = psi;
		return angles;
	}
	if (_RPTType == 8){
		do{
			b = 124 + (rand() % 57);
			phi = (rand() % 180)+1;
			psi = -1 * phi + b;
		}while(psi > 180 || psi < -180);
		angles[0] = phi; angles[1] = psi;
		return angles;
	}
	if (_RPTType == 9){
		do{
			b = 95 + (rand() % 29);
			phi = (rand() % 180)+1;
			psi = -1 * phi + b;
		}while(psi > 180 || psi < -180);
		angles[0] = phi; angles[1] = psi;
		return angles;
	}
	if (_RPTType == 10){
		do{
			b = 58 + (rand() % 37);
			phi = (rand() % 180)+1;
			psi = -1 * phi + b;
		}while(psi > 180 || psi < -180);
		angles[0] = phi; angles[1] = psi;
		return angles;
	}
	if (_RPTType == 11){
		do{
			b = (rand() % 58);
			phi = (rand() % 180)+1;
			psi = -1 * phi + b;
		}while(psi > 180 || psi < -180);
		angles[0] = phi; angles[1] = psi;
		return angles;
	}
	if (_RPTType == 12){
		do{
			b = -58 + (rand() % 58);
			phi = (rand() % 180)+1;
			psi = -1 * phi + b;
		}while(psi > 180 || psi < -180);
		angles[0] = phi; angles[1] = psi;
		return angles;
	}
	if (_RPTType == 13){
		do{
			b = -95 + (rand() % 37);
			phi = (rand() % 180)+1;
			psi = -1 * phi + b;
		}while(psi > 180 || psi < -180);
		angles[0] = phi; angles[1] = psi;
		return angles;
	}
	if (_RPTType == 14){
		do{
			b = -124 + (rand() % 29);
			phi = (rand() % 180)+1;
			psi = -1 * phi + b;
		}while(psi > 180 || psi < -180);
		angles[0] = phi; angles[1] = psi;
		return angles;
	}
	if (_RPTType == 15){
		do{
			b = -180 + (rand() % 56);
			phi = (rand() % 180)+1;
			psi = -1 * phi + b;
		}while(psi > 180 || psi < -180);
		angles[0] = phi; angles[1] = psi;
		return angles;
	}
	return angles;
}
vector <double> protein::getRandConformationFromBackboneType(double _phi, double _psi)
{
	double RPT = getResiduesPerTurn(_phi,_psi);
	UInt RPTType = getBackboneSequenceType(RPT,_phi);
	return getRandPhiPsifromBackboneSequenceType(RPTType);
}
void protein::optimizeRotamers()
{
	vector < UIntVec > activePositions;
	vector < UIntVec > rotamerArray;
	activePositions.resize(0);
	rotamerArray.resize(0);
	for (UInt i = 0; i < itsIndependentChainsMap.size(); i ++)
	{
		UInt indChain = itsIndependentChainsMap[i];
		UIntVec activeRes = itsChains[indChain]->getActiveResidues();
		for (UInt j = 0; j < activeRes.size(); j ++)
		{
			UIntVec tempRot = itsChains[indChain]->getAllowedRotamers(activeRes[j], getTypeFromResNum(indChain, activeRes[j]), 0);
			if (tempRot.size() > 1)
			{
				UIntVec position;
				position.resize(0);
				position.push_back(indChain);
				position.push_back(activeRes[j]);
				activePositions.push_back(position);
			}
			else if (tempRot.size() == 1)
			{
				setRotamer(indChain, activeRes[j],0,tempRot[0]);
			}
		}
	}
	if (activePositions.size() > 0)
	{
		rotamerArray = rotamerDEE(activePositions);
		vector < UIntVec > tempActivePositions;
		tempActivePositions.resize(0);
		vector < UIntVec > tempRotamerArray;
		tempRotamerArray.resize(0);
		for (UInt i = 0; i < activePositions.size(); i ++)
		{
			if (rotamerArray[i].size() > 1)
			{
				tempRotamerArray.push_back(rotamerArray[i]);
				tempActivePositions.push_back(activePositions[i]);
			}
			else if (rotamerArray[i].size() == 1)
			{
				setRotamer(activePositions[i][0], activePositions[i][1], 0, rotamerArray[i][0]);
			}
		}
		if (tempRotamerArray.size() != 0)
		{
			optimizeRotamers(tempActivePositions, tempRotamerArray);
		}
		else  // only one solution from DEE - map that onto the protein
		{
			for (UInt i = 0; i < activePositions.size(); i ++)
			{
				setRotamer(activePositions[i][0], activePositions[i][1], 0, rotamerArray[i][0]);
			}
		}
	}
	else
	{
		cout << "ERROR:  no positions to optimize rotamers for." << endl;
	}
	return;
}

void protein::optimizeRotamers(vector <UIntVec> _activePositions)
{
	vector <UIntVec> activePositions;
	activePositions.resize(0);
	for (UInt i = 0; i < _activePositions.size(); i ++)
	{
		UIntVec tempRot = itsChains[_activePositions[i][0]]->getAllowedRotamers(_activePositions[i][1], getTypeFromResNum(_activePositions[i][0], _activePositions[i][1]), 0);
		if (tempRot.size() > 1)
		{
			activePositions.push_back(_activePositions[i]);
		}
		else if (tempRot.size() == 1)
		{
			setRotamer(_activePositions[i][0], _activePositions[i][1],0 , tempRot[0]);
		}
	}

	vector < UIntVec > rotamerArray;
	rotamerArray.resize(0);
	if (activePositions.size() > 0)
	{
		rotamerArray = rotamerDEE(activePositions);
		vector < UIntVec > tempActivePositions;
		tempActivePositions.resize(0);
		vector < UIntVec > tempRotamerArray;
		tempRotamerArray.resize(0);
		for (UInt i = 0; i < activePositions.size(); i ++)
		{
			if (rotamerArray[i].size() > 1)
			{
				tempRotamerArray.push_back(rotamerArray[i]);
				tempActivePositions.push_back(activePositions[i]);
			}
			else if (rotamerArray[i].size() == 1)
			{
				setRotamer(activePositions[i][0], activePositions[i][1], 0, rotamerArray[i][0]);
			}
		}
		if (tempRotamerArray.size() != 0)
		{
			optimizeRotamers(tempActivePositions, tempRotamerArray);
		}
		else  // only one solution from DEE - map that onto the protein
		{
			for (UInt i = 0; i < activePositions.size(); i ++)
			{
				setRotamer(activePositions[i][0], activePositions[i][1], 0, rotamerArray[i][0]);
			}
		}
	}
	else
	{
		cout << "ERROR:  no positions to optimize rotamers for." << endl;
	}
	return;
}

void protein::optimizeRotamers(vector <UIntVec> _activePositions, vector <UIntVec> _rotamerArray)
{
	double lowestEnergy = 1E10;
	UIntVec currentRotamerArray;
	currentRotamerArray.resize(0);
	UIntVec bestRotamerArray;
	// initialize bestRotamerArray
	bestRotamerArray.resize(0);
	for (UInt i = 0; i < _rotamerArray.size(); i ++)
	{
		bestRotamerArray.push_back(_rotamerArray[i][0]);
	}
	// run through rotamer combinations ...
	for (UInt i = 0; i < _rotamerArray[0].size(); i ++)
	{
		UInt index = 0;
		setRotamer(_activePositions[0][0], _activePositions[0][1], 0, _rotamerArray[0][i]);
		if (_activePositions.size() != 1) // if the first is not the only position
		{
			currentRotamerArray.push_back(_rotamerArray[0][i]);
			bestRotamerArray = getEnergySurface(_activePositions, _rotamerArray, currentRotamerArray, bestRotamerArray, index + 1, lowestEnergy);
			currentRotamerArray.pop_back();
		}
	}
	// set protein to best rotamer set
	if (bestRotamerArray.size() == _activePositions.size())
	{
		for (UInt i = 0; i < _activePositions.size(); i ++)
		{
			setRotamer(_activePositions[i][0], _activePositions[i][1], 0, bestRotamerArray[i]);
		}
		// optimizeSmallRotations(_activePositions);
	}
	else
	{
		cout << "ERROR in optimizeRotamers... size of best array doesnt match the number of active positions" << endl;
	}
	return;
}

UIntVec protein::getEnergySurface(vector <UIntVec> _activePositions, vector <UIntVec> _rotamerArray, UIntVec _currentArray, UIntVec _bestArray, UInt _index, double& _lowestEnergy)
{
	if (_index < _activePositions.size())
	{
		for (UInt i = 0; i < _rotamerArray[_index].size(); i ++)
		{
			setRotamer(_activePositions[_index][0], _activePositions[_index][1], 0, _rotamerArray[_index][i]);
			_currentArray.push_back(_rotamerArray[_index][i]);
			_bestArray = getEnergySurface(_activePositions, _rotamerArray, _currentArray, _bestArray, _index + 1, _lowestEnergy);
			_currentArray.pop_back();
		}
		return _bestArray;
	}
	else
	{
		double energy = intraEnergy();
		if (messagesActive)
		{
			for (UInt i = 0; i < _currentArray.size(); i ++ )
			{
				cout << _currentArray[i] << " ";
			}
			cout << "\tLOW:  " << _lowestEnergy << " CUR:  " << energy << endl;
		}
		if (_lowestEnergy > energy)
		{
			_lowestEnergy = energy;
			_bestArray = _currentArray;
		}
		return _bestArray;
	}
}

vector <UIntVec> protein::rotamerDEE()
{
	// create list of active positions
	vector < UIntVec > activePositions;
	vector < UIntVec > rotamerArray;
	for (UInt i = 0; i < itsIndependentChainsMap.size(); i ++)
	{
		UInt indChain = itsIndependentChainsMap[i];
		UIntVec activeRes = itsChains[i]->getActiveResidues();
		for (UInt j = 0; j < activeRes.size(); j ++)
		{
			UIntVec tempRot = itsChains[indChain]->getAllowedRotamers(activeRes[j], getTypeFromResNum(indChain, activeRes[j]), 0);
			if (tempRot.size() > 1)
			{
				UIntVec position;
				position.resize(0);
				position.push_back(indChain);
				position.push_back(activeRes[j]);
				activePositions.push_back(position);
			}
		}
	}
	rotamerArray = rotamerDEE(activePositions);
	return rotamerArray;
}

vector < UIntVec > protein::rotamerDEE(vector <UIntVec> _activePositions)
{
	vector <UIntVec> rotamerArray;
	rotamerArray.resize(0);
	vector < vector < bool > > flagArray;
	flagArray.resize(0);
	for (UInt i = 0; i < _activePositions.size(); i ++)
	{
		rotamerArray.push_back(itsChains[_activePositions[i][0]]->getAllowedRotamers(_activePositions[i][1], getTypeFromResNum(_activePositions[i][0], _activePositions[i][1]), 0));
		vector < bool > tempFlagList;
		for (UInt j = 0; j < rotamerArray[i].size(); j ++)
		{
			tempFlagList.push_back(true);
		}
		flagArray.push_back(tempFlagList);
	}
	UInt count = 1;
	for (UInt i = 0; i < rotamerArray.size(); i ++)
	{
		count *= rotamerArray[i].size();
	}
	if (messagesActive) cout << "Starting number of rotamers: " << count << endl;
	double Eir, Eit, Esum, Emin;
	for (UInt i = 0; i < _activePositions.size(); i ++)
	{
		UInt chain1 = _activePositions[i][0];
		UInt res1 = _activePositions[i][1];
		for (UInt r = 0; r < rotamerArray[i].size(); r ++)
		{
			bool irDE = false;
			setRotamer(chain1, res1, 0, rotamerArray[i][r]);
			Eir = getSelfEnergy(chain1, res1);
			if (Eir > 1E3) irDE=true;
			if (!irDE)
			{
				for (UInt t = 0; t < rotamerArray[i].size(); t ++)
				{
					if (r != t && !irDE)
					{
						setRotamer(chain1,res1, 0, rotamerArray[i][t]);
						Eit = getSelfEnergy(chain1, res1);
						Esum = 0.0;
						for (UInt j = 0; j < _activePositions.size(); j++)
						{
							UInt chain2 = _activePositions[j][0];
							UInt res2 = _activePositions[j][1];

							UInt indChainIndex1 = 0;
							UInt indChainIndex2 = 0;

							for (UInt n = 0; n < itsIndependentChainsMap.size(); n ++)
							{
								if (itsIndependentChainsMap[n] == chain1) indChainIndex1 = n;
								if (itsIndependentChainsMap[n] == chain2) indChainIndex2 = n;
							}
							UInt numSymLinkedChains1 = itsChainLinkageMap[indChainIndex1].size();
							UInt numSymLinkedChains2 = itsChainLinkageMap[indChainIndex2].size();
							if (i != j)
							{
								Emin = 1E20;
								for (UInt s = 0; s < rotamerArray[j].size(); s ++)
								{
									setRotamer(chain2, res2, 0, rotamerArray[j][s]);
									setRotamer(chain1, res1, 0, rotamerArray[i][r]);
									double Etemp = getResPairEnergy(chain1, res1, chain2, res2);
									for (UInt n = 0; n < numSymLinkedChains2; n ++)
									{
										if (itsChainLinkageMap[indChainIndex2][n] != -1)
										{
											UInt chain2sym = itsChainLinkageMap[indChainIndex2][n];
											Etemp += getResPairEnergy(chain1, res1, chain2sym, res2);
										}
									}
									setRotamer(chain1, res1, 0, rotamerArray[i][t]);
									Etemp -= getResPairEnergy(chain1, res1, chain2, res2);
									for (UInt n = 0; n < numSymLinkedChains2; n ++)
									{
										if (itsChainLinkageMap[indChainIndex2][n] != -1)
										{
											UInt chain2sym = itsChainLinkageMap[indChainIndex2][n];
											Etemp -= getResPairEnergy(chain1, res1, chain2sym, res2);
										}
									}
									if (Emin > Etemp) Emin = Etemp;
								}
								Esum += Emin;
							}
							else if (i == j)
							{
								Emin = 1E20;
								setRotamer(chain1, res1, 0, rotamerArray[i][r]);
								double Etemp = 0.0;
								for (UInt n = 0; n < numSymLinkedChains1; n ++)
								{
									if (itsChainLinkageMap[indChainIndex1][n] != -1)
									{
										UInt chain1sym = itsChainLinkageMap[indChainIndex1][n];
										Etemp += getResPairEnergy(chain1, res1, chain1sym, res1);
									}
								}
								setRotamer(chain1, res1, 0, rotamerArray[i][t]);
								for (UInt n = 0; n < numSymLinkedChains1; n ++)
								{
									if (itsChainLinkageMap[indChainIndex1][n] != -1)
									{
										UInt chain1sym = itsChainLinkageMap[indChainIndex1][n];
										Etemp -= getResPairEnergy(chain1, res1, chain1sym, res1);
									}
								}
								if (Emin > Etemp) Emin = Etemp;
								Esum += Emin;
							}
						}
						double total = Eir - Eit + Esum;
						if (total > 0)
						{
							irDE = true; // ir is dead endin
						}
					}
				}
			}
			if  (irDE)
			{
				if (messagesActive) cout << "residue " << chain1 << " " << res1 << " rotamer " << rotamerArray[i][r] << " is deadEnding" <<endl;
				flagArray[i][r] = false;
			}
		}
	}
	// generate final rotamer array
	vector <UIntVec> finalArray;
	finalArray.resize(0);
	for (UInt i = 0; i < rotamerArray.size(); i ++)
	{
		UIntVec tempArray;
		tempArray.resize(0);
		for (UInt j = 0; j < rotamerArray[i].size(); j ++)
		{
			if (flagArray[i][j])
			{
				tempArray.push_back(rotamerArray[i][j]);
			}
		}
		if (tempArray.size() == 0) // DEE failed ...
		{
			for (UInt j = 0; j < rotamerArray[i].size(); j ++)
			{
				tempArray.push_back(rotamerArray[i][j]);
			}
		}
		finalArray.push_back(tempArray);
	}

	count = 1;
	for (UInt i = 0; i < finalArray.size(); i ++)
	{
		if (finalArray[i].size() != 0) count *= finalArray[i].size();
		else
		{
			cout << "ERROR: " << _activePositions[i][0] << " " << _activePositions[i][1] << " has zero rotamers.  Error\n";
			exit(1);
		}
	}
	if (messagesActive) cout << "DEE:  ending number of rotamers = " << count << endl;
	return finalArray;
}

void protein::optimizeSmallRotations(UInt _steps, double _stepSize)
{
	vector < UIntVec > activePositions;
	for (UInt i = 0; i < itsIndependentChainsMap.size(); i ++)
	{
		UInt indChain = itsIndependentChainsMap[i];
		UIntVec activeRes = itsChains[indChain]->getActiveResidues();
		for (UInt j = 0; j < activeRes.size(); j ++)
		{
			UInt numChis = itsChains[indChain]->getNumChis(j,0);
			if (numChis > 0)
			{
				UIntVec position;
				position.resize(0);
				position.push_back(indChain);
				position.push_back(activeRes[j]);
				activePositions.push_back(position);
			}
		}
	}
	optimizeSmallRotations(activePositions, _steps, _stepSize);
	return;
}

void protein::optimizeSmallRotations(vector <UIntVec> _activePositions, UInt _steps, double _stepSize)
{
	vector <UIntVec> activePositions;
	activePositions.resize(0);
	for (UInt i = 0; i < _activePositions.size(); i ++)
	{
		UIntVec tempRot = itsChains[_activePositions[i][0]]->getAllowedRotamers(_activePositions[i][1], getTypeFromResNum(_activePositions[i][0], _activePositions[i][1]), 0);
		if (tempRot.size() > 1)
		{
			activePositions.push_back(_activePositions[i]);
		}
	}

	UInt activePos = 0;
	double lowestEnergy = 1E20;
	vector < vector < double > > bestChiArray;
	bestChiArray.resize(0);
	for (UInt i = 0; i < activePositions.size(); i ++)
	{
		vector <vector <double> > dihedrals = itsChains[activePositions[i][0]]->getDihedrals(activePositions[i][1]);
		bestChiArray.push_back(dihedrals[0]);
	}
	UInt numChis = itsChains[activePositions[activePos][0]]->getNumChis(activePositions[activePos][1],0);
	if (numChis > 0)
	{
		for (UInt chiPos = 0; chiPos < numChis; chiPos ++)
		{
			for (UInt step = 0; step <= _steps; step ++) // to do last pos, need one more step
			{
				if (step == 0) // set dihedral at first position
				{
                    double angle = -0.5 * _steps*(double)_stepSize; //changed div to mult "-1*_steps*(double)_stepSize/2.0;"
					itsChains[activePositions[activePos][0]]->setRelativeChi(activePositions[activePos][1], 0, chiPos, angle);
				}
				else          // increment rotation
				{
					itsChains[activePositions[activePos][0]]->setRelativeChi(activePositions[activePos][1], 0, chiPos, _stepSize);
				}
				if (activePos < activePositions.size())
				{
					bestChiArray = getRotationEnergySurface(activePositions, _steps, _stepSize, activePos + 1, bestChiArray, lowestEnergy);
				}
				else
				{
					double energy = intraEnergy();
					if (lowestEnergy > energy)
					{
						lowestEnergy = energy;
						bestChiArray.resize(0);
						for (UInt i = 0; i < activePositions.size(); i ++)
						{
							vector <vector <double> > dihedrals = itsChains[activePositions[i][0]]->getDihedrals(activePositions[i][1]);
							bestChiArray.push_back(dihedrals[0]);
						}
					}
				}
			}
		}
	}
	for (UInt i = 0; i < activePositions.size(); i ++)
	{
		itsChains[activePositions[i][0]]->setDihedrals(activePositions[i][1], 0, bestChiArray[i]);
	}
	return;
}

vector < vector < double > > protein::getRotationEnergySurface(vector < UIntVec > _active, UInt _steps, double _stepSize, UInt _activePos, vector <vector<double> > _bestChiArray, double &_lowestEnergy)
{
	if (_activePos < _active.size())
	{
		UInt numChis = itsChains[_active[_activePos][0]]->getNumChis(_active[_activePos][1],0);
		if (numChis > 0)
		{
			for (UInt chiPos = 0; chiPos < numChis; chiPos ++)
			{
				for (UInt step = 0; step <= _steps; step ++)
				{
					if (step == 0) // set dihedral at first position
					{
                        double angle = -0.5 * _steps*(double)_stepSize; //changed div to mult "-1*_steps*(double)_stepSize/2.0;"
						itsChains[_active[_activePos][0]]->setRelativeChi(_active[_activePos][1], 0, chiPos, angle);
					}
					else          // increment rotation
					{
						itsChains[_active[_activePos][0]]->setRelativeChi(_active[_activePos][1], 0, chiPos, _stepSize);
					}
					_bestChiArray = getRotationEnergySurface(_active, _steps, _stepSize, _activePos + 1, _bestChiArray, _lowestEnergy);
				}
			}
		}
		else
		{
			_bestChiArray = getRotationEnergySurface(_active, _steps, _stepSize, _activePos + 1, _bestChiArray, _lowestEnergy);
		}
	}
	else
	{
		double energy = intraEnergy();
		if (_lowestEnergy > energy)
		{
			_lowestEnergy = energy;
			_bestChiArray.resize(0);
			if (messagesActive) cout << "bestEnergy ... " << _lowestEnergy << endl;
			for (UInt i = 0; i < _active.size(); i ++)
			{
				vector <vector <double> > dihedrals = itsChains[_active[i][0]]->getDihedrals(_active[i][1]);
				_bestChiArray.push_back(dihedrals[0]);
				//listDihedrals();
			}
		}
	}
	return _bestChiArray;
}

void protein::saveState(string _fileName)
{
	if (messagesActive) cout << " writing file " << _fileName << endl;
	pdbWriter(this, _fileName);
	return;
}

void protein::saveState(string& _fileName)
{
	if (messagesActive) cout << " writing file " << _fileName << endl;
	pdbWriter(this, _fileName);
	return;
}

void protein::optimizeSmallRotations( UIntVec _activePosition, UInt _steps, double _stepSize)
{

	UIntVec tempRot = itsChains[_activePosition[0]]->getAllowedRotamers(_activePosition[1], getTypeFromResNum(_activePosition[0], _activePosition[1]), 0);
	if (tempRot.size() <= 1)
	{
		if (messagesActive) cout << "no small rotations to optimize.  " << endl;
		return;
	}

	double lowestEnergy = 1E20;
	vector < double > bestChiArray;
	bestChiArray.resize(0);
	vector < vector < double > > tempDihedral;
	tempDihedral = itsChains[_activePosition[0]]->getDihedrals(_activePosition[1]);
	bestChiArray = tempDihedral[0];
	UInt numChis = itsChains[_activePosition[0]]->getNumChis(_activePosition[1],0);
	if (numChis > 0)
	{
		UInt chiPos = 0;
		for (UInt step = 0; step <= _steps; step ++) // to do last pos, need one more step
		{
			if (step == 0) // set dihedral at first position
			{
                double angle = -0.5 * _steps*(double)_stepSize; //changed div to mult "-1*_steps*(double)_stepSize/2.0;"
				itsChains[_activePosition[0]]->setRelativeChi(_activePosition[1], 0, chiPos, angle);
			}
			else          // increment rotation
			{
				itsChains[_activePosition[0]]->setRelativeChi(_activePosition[1], 0, chiPos, _stepSize);
			}
			if (chiPos < numChis)
			{
				bestChiArray = getRotationEnergySurface(_activePosition, _steps, _stepSize, chiPos + 1, bestChiArray, lowestEnergy);
			}
			else
			{
				double energy = getPositionEnergy(_activePosition);
				if (lowestEnergy > energy)
				{
					lowestEnergy = energy;
					bestChiArray.resize(0);
					vector <vector <double> > dihedrals = itsChains[_activePosition[0]]->getDihedrals(_activePosition[1]);
					bestChiArray = dihedrals[0];
				}
			}
		}
	}
	itsChains[_activePosition[0]]->setDihedrals(_activePosition[1], 0, bestChiArray);
	return;
}

vector < double >  protein::getRotationEnergySurface(UIntVec  _active, UInt _steps, double _stepSize, UInt _chiPos, vector <double>  _bestChiArray, double &_lowestEnergy)
{
	UInt numChis = itsChains[_active[0]]->getNumChis(_active[1],0);
	if (_chiPos < numChis)
	{
		for (UInt step = 0; step <= _steps; step ++)
		{
			if (step == 0) // set dihedral at first position
			{
                double angle = -0.5 * _steps*(double)_stepSize; //changed div to mult "-1*_steps*(double)_stepSize/2.0;"
				itsChains[_active[0]]->setRelativeChi(_active[1], 0, _chiPos, angle);
			}
			else          // increment rotation
			{
				itsChains[_active[0]]->setRelativeChi(_active[1], 0, _chiPos, _stepSize);
			}
			_bestChiArray = getRotationEnergySurface(_active, _steps, _stepSize, _chiPos + 1, _bestChiArray, _lowestEnergy);
		}
	}
	else
	{
		double energy = getPositionEnergy(_active);
		vector <vector <double> > dihedrals = itsChains[_active[0]]->getDihedrals(_active[1]);
		if (_lowestEnergy > energy)
		{
			_lowestEnergy = energy;
			_bestChiArray.resize(0);
			if (messagesActive) cout << _active[0] << " " << _active[1] << " bestEnergy ... " << _lowestEnergy;
			_bestChiArray = dihedrals[0];
			if (messagesActive)
			{
				for (UInt i = 0; i < _bestChiArray.size(); i ++)
				{
					cout << " " << _bestChiArray[i];
				}
				cout << endl;
			}
		}
	}
	return _bestChiArray;
}

double protein::getHBondEnergy(const UInt _chain1, const UInt _res1, const UInt _chain2, const UInt _res2)
{
    enum aminoAcid {A, R, N, D, C, Q, E, G, H, I, L, K, M, F, P, S, T, W, Y, V};

    UInt type1 = getTypeFromResNum(_chain1, _res1);
    UInt type2 = getTypeFromResNum(_chain2, _res2);

    vector <dblVec> donorList(0);
    vector <dblVec> donorBaseList(0);
    vector <dblVec> acceptorList(0);
    vector <dblVec> acceptorBaseList(0);
    vector <double> donorAngleList(0);
    vector <double> acceptorAngleList(0);


    for (UInt i = 1; i <= 2; i ++)
    {
        UInt chain, res, type;
        if (i == 1)
        {
            chain = _chain1;
            res = _res1;
            type = type1;
        }
        if (i == 2)
        {
            chain = _chain2;
            res = _res2;
            type = type2;
        }

        if (type == H)
        {
            dblVec donor = getCoords(chain,res,"NE2");
            dblVec temp1 = getCoords(chain,res,"CE1");
            dblVec temp2 = getCoords(chain,res,"CD2");
            dblVec donorBase = (temp1 + temp2) * 0.5; // converted div to multi "donorBase = (temp1 + temp2)/2.0;"
            dblVec acceptor = getCoords(chain,res,"ND1");
            temp2 = getCoords(chain,res,"CG");
            dblVec acceptorBase = (temp1 + temp2) * 0.5; // converted div to mult "acceptorBase = (temp1 + temp2) / 2.0;"
            donorList.push_back(donor);
            donorBaseList.push_back(donorBase);
            acceptorList.push_back(acceptor);
            acceptorBaseList.push_back(acceptorBase);
            donorAngleList.push_back(180.0);
            acceptorAngleList.push_back(180.0);
        }
        if (type == T)
        {
            dblVec donor = getCoords(chain,res,"OG1");
            dblVec donorBase = getCoords(chain,res,"CB");
            dblVec acceptor = getCoords(chain,res,"OG1");
            dblVec acceptorBase = getCoords(chain,res,"CB");
            donorList.push_back(donor);
            donorBaseList.push_back(donorBase);
            acceptorList.push_back(acceptor);
            acceptorBaseList.push_back(acceptorBase);
            donorAngleList.push_back(109.0);
            acceptorAngleList.push_back(109.0);
        }
        if (type == S)
        {
            dblVec donor = getCoords(chain,res,"OG");
            dblVec donorBase = getCoords(chain,res,"CB");
            dblVec acceptor = getCoords(chain,res,"OG");
            dblVec acceptorBase = getCoords(chain,res,"CB");
            donorList.push_back(donor);
            donorBaseList.push_back(donorBase);
            acceptorList.push_back(acceptor);
            acceptorBaseList.push_back(acceptorBase);
            donorAngleList.push_back(109.0);
            acceptorAngleList.push_back(109.0);
        }
        if (type == Y)
        {
            dblVec donor = getCoords(chain,res,"OH");
            dblVec donorBase = getCoords(chain,res,"CZ");
            dblVec acceptor = getCoords(chain,res,"OH");
            dblVec acceptorBase = getCoords(chain,res,"CZ");
            donorList.push_back(donor);
            donorBaseList.push_back(donorBase);
            acceptorList.push_back(acceptor);
            acceptorBaseList.push_back(acceptorBase);
            donorAngleList.push_back(109.0);
            acceptorAngleList.push_back(109.0);
        }
      if (type == N)
        {
            dblVec donor = getCoords(chain,res,"ND2");
            dblVec donorBase = getCoords(chain,res,"CG");
            dblVec acceptor = getCoords(chain,res,"OD1");
            dblVec acceptorBase = getCoords(chain,res,"CG");
            donorList.push_back(donor);
            donorBaseList.push_back(donorBase);
            acceptorList.push_back(acceptor);
            acceptorBaseList.push_back(acceptorBase);
            donorAngleList.push_back(120.0);
            acceptorAngleList.push_back(180.0);
        }
        if (type == Q)
        {
            dblVec donor = getCoords(chain,res,"NE2");
            dblVec donorBase = getCoords(chain,res,"CD");
            dblVec acceptor = getCoords(chain,res,"OOE1");
            dblVec acceptorBase = getCoords(chain,res,"CD");
            donorList.push_back(donor);
            donorBaseList.push_back(donorBase);
            acceptorList.push_back(acceptor);
            acceptorBaseList.push_back(acceptorBase);
            donorAngleList.push_back(120.0);
            acceptorAngleList.push_back(180.0);
        }
        if (type == D)
        {
            dblVec acceptor1 = getCoords(chain,res,"OD1");
            dblVec acceptor2 = getCoords(chain,res,"OD2");
            dblVec acceptorBase = getCoords(chain,res,"CG");
            acceptorList.push_back(acceptor1);
            acceptorList.push_back(acceptor2);
            acceptorBaseList.push_back(acceptorBase);
            acceptorBaseList.push_back(acceptorBase);
            acceptorAngleList.push_back(180.0);
            acceptorAngleList.push_back(180.0);
        }
        if (type == E)
        {
            dblVec acceptor1 = getCoords(chain,res,"OE1");
            dblVec acceptor2 = getCoords(chain,res,"OE2");
            dblVec acceptorBase = getCoords(chain,res,"CD");
            acceptorList.push_back(acceptor1);
            acceptorList.push_back(acceptor2);
            acceptorBaseList.push_back(acceptorBase);
            acceptorBaseList.push_back(acceptorBase);
            acceptorAngleList.push_back(180.0);
            acceptorAngleList.push_back(180.0);
        }
        if (type == W)
        {
            dblVec donor = getCoords(chain,res,"NE1");
            dblVec temp1 = getCoords(chain,res,"CD1");
            dblVec temp2 = getCoords(chain,res,"CE2");
            dblVec donorBase = (temp1 + temp2) * 0.5; // converted div to multi "donorBase = (temp1 + temp2)/2.0;"
            donorList.push_back(donor);
            donorBaseList.push_back(donorBase);
            donorAngleList.push_back(180.0);
        }
        dblVec carbonyl = getCoords(chain,res,"O");
        dblVec carbonylC = getCoords(chain,res,"C");
        acceptorList.push_back(carbonyl);
        acceptorBaseList.push_back(carbonylC);
        acceptorAngleList.push_back(120);
    }
    double energy = 0.0;
    if (donorList.size() == 0) return 0.0;

   for (UInt i = 0; i < donorList.size(); i ++)
    {
        for (UInt j = 0; j < acceptorList.size(); j ++)
        {
            dblVec DDB = donorList[i] - donorBaseList[i];
            dblVec DA = donorList[i] - acceptorList[j];
            dblVec AAB = acceptorList[j] - acceptorBaseList[j];

            double magDDB = sqrt(CMath::dotProduct(DDB,DDB));
            double magDA = sqrt(CMath::dotProduct(DA,DA));
            double magAAB = sqrt(CMath::dotProduct(AAB,AAB));

            double thisE;
            if (magDA < 1e-50) thisE = 0.0;
            else
            {

                double donorAngle = acos( CMath::dotProduct(DA,DDB)/(magDDB*magDA) );
                double acceptorAngle = acos( CMath::dotProduct(DA,AAB)/(magDA*magAAB) );

                double distRatio = 2.8 / magDA;
                thisE = 2.0 *(5.0*pow(distRatio,12.0)-6.0*pow(distRatio,10.0));
                double angleFactor1 = cos(donorAngleList[i]*PI * 0.00555555555-donorAngle); //converted div to mult "cos(donorAngleList[i]*PI/180.0-donorAngle)"
                double angleFactor2 = cos(acceptorAngleList[j]*PI * 0.00555555555-acceptorAngle);//converted div to mult "cos(acceptorAngleList[j]*PI/180.0-acceptorAngle)"
                thisE = thisE * pow(angleFactor1, 2.0) * pow (angleFactor2, 2.0);
            }
            energy += thisE;
        }
    }
    return energy;
}
//this function will calculate the hammingdistance between two sequences obtained from two pdb files
double protein::getHammingDistance(vector<string>seq1,vector<string>seq2) 
{ 
	double count = 0.0;double percent=0.0;double countsim=0.0;
	for (UInt i=0;i<seq1.size();i++){
        	if(seq1[i]==seq2[i]){
                	if (seq1[i]!="-" || seq2[i]!="X"){
                        	 count++;
                }
         }
        }
	countsim=seq1.size()-count;
	percent=((countsim/double(seq1.size()))*100);
   	return percent;
}  


