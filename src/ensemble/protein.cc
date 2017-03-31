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
    energies.clear();
    energies.resize(0);
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
    energies.clear();
    energies.resize(0);
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

//******************testing junk****************
void protein::accessChainZeroResZero()
{
	itsChains[0]->accessResZero();
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

double protein::getPositionSoluteEnergy(vector <int> _position)
{
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
       double tempE = itsChains[_position[1]]->getPositionIntraSoluteEnergy(_position);
       intraEnergy += tempE;
       for (UInt i = 0; i < itsChains.size(); i++)
       {
           if ((int)i != _position[1])
           {
               tempE = itsChains[_position[1]]->getPositionInterSoluteEnergy(_position, itsChains[i]);
               intraEnergy += tempE;
           }
       }
   }
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

int protein::setAngleLocal(const UInt _chain, const UInt _res, double _angle, double deltaTheta, UInt angleType, int distance, int direction)
{
	if (_chain < itsChains.size())
	{
		return itsChains[_chain]->setAngleLocal(_res, _angle, deltaTheta, angleType, distance, direction);
	}
	else
	{
		cout << "chain index out of range" << endl;
		return -1;
	}
}

int protein::setDihedralLocal(const UInt _chainIndex, const UInt _resIndex, double _deltaTheta, UInt _angleType)
{
	if (_chainIndex < itsChains.size())
	{
		return itsChains[_chainIndex]->setDihedralLocal(_resIndex, _deltaTheta, _angleType);
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

double protein::getPositionSoluteEnergy(UInt _chain, UInt _residue, bool _updateDielectrics)
{
	if (_updateDielectrics)
	{
        this->updatePositionDielectrics(_chain, _residue);
	}
	vector <int> position;
	position.resize(0);
	position.push_back(-1);
	position.push_back((int)_chain);
	position.push_back((int)_residue);
	double positionEnergy=getPositionSoluteEnergy(position);
	return positionEnergy;
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

double protein::intraSoluteEnergy(bool _updateDielectrics)
{
	double intraEnergy = 0.0;
	if (_updateDielectrics)
	{
        this->updateDielectrics();
	}
	for(UInt i=0; i<itsChains.size(); i++)
	{
		intraEnergy += itsChains[i]->intraSoluteEnergy();
		for(UInt j=i+1; j<itsChains.size(); j++)
		{
			intraEnergy += itsChains[i]->interSoluteEnergy(itsChains[j]);

		}
	}
	return intraEnergy;
}

double protein::intraSoluteEnergy(bool _updateDielectrics, UInt _activeChain)
{
    double intraEnergy = 0.0;
    if (_updateDielectrics)
    {
        this->updateDielectrics();
    }
    intraEnergy += itsChains[_activeChain]->intraSoluteEnergy();
    for(UInt j=0; j<itsChains.size(); j++)
    {
        if (j != _activeChain)
        {
            intraEnergy += itsChains[_activeChain]->interSoluteEnergy(itsChains[j]);
        }
    }
    return intraEnergy;
}

double protein::interSoluteEnergy(bool _updateDielectrics, UInt _chain1, UInt _chain2)
{
	if (_updateDielectrics)
	{
        this->updateDielectrics();
	}
	double interEnergy = itsChains[_chain1]->interSoluteEnergy(itsChains[_chain2]);
	return interEnergy;
}

vector <double> protein::calculateDielectric(UInt _chainIndex, UInt _residueIndex, UInt _atomIndex)
{
    vector <double> polarization(2);
    vector <double> _polarization(2);
    vector <double> dielectric(2);
    double waterPol = residueTemplate::getPolarizability(52);
    double waterVol = residueTemplate::getVolume(52);
    double totalVol = residue::cutoffCubeVolume;
    polarization[0] = 0.0;
    polarization[1] = 0.0;

    // get volume and polarizabilities through protein around atom
	for(UInt i=0; i<itsChains.size(); i++)
	{
        _polarization = itsChains[_chainIndex]->calculateDielectric(itsChains[i], _residueIndex, _atomIndex);
        polarization[0] += _polarization[0];
        polarization[1] += _polarization[1];
	}

    // calculate local dielectric for atom
	double totalWaterVol = totalVol-(polarization[0]/1.89);
    int waters = totalWaterVol/waterVol;
    double totalWaterPol = waters*waterPol;
	double die = 1+4*3.14*((waters)/totalVol)*(totalWaterPol+polarization[1]);
    if (die < 2) { die = 2.0;}
    if (die > 78) { die = 78.0;}
    dielectric[0] = die;
    dielectric[1] = waters;
    return dielectric;
}

vector <double> protein::calculateDielectric(chain* _chain, residue* _residue, atom* _atom)
{
    vector <double> polarization(2);
    vector <double> _polarization(2);
    vector <double> dielectric(2);
    double waterPol = residueTemplate::getPolarizability(52);
    double waterVol = residueTemplate::getVolume(52);
    double totalVol = residue::cutoffCubeVolume;
    polarization[0] = 0.0;
    polarization[1] = 0.0;

    // get volume and polarizabilities through protein around atom
    for(UInt i=0; i<itsChains.size(); i++)
    {
        _polarization = _chain->calculateDielectric(itsChains[i], _residue, _atom);
        polarization[0] += _polarization[0];
        polarization[1] += _polarization[1];
    }

    // calculate local dielectric for atom
	double totalWaterVol = totalVol-(polarization[0]/1.89);
    int waters = totalWaterVol/waterVol;
    double totalWaterPol = waters*waterPol;
	double die = 1+4*3.14*((waters)/totalVol)*(totalWaterPol+polarization[1]);
    if (die < 2) { die = 2.0;}
    if (die > 78) { die = 78.0;}
    dielectric[0] = die;
    dielectric[1] = waters;
    return dielectric;
}

vector <double> protein::calculateChainIndependentDielectric(chain* _chain, residue* _residue, atom* _atom)
{
    vector <double> polarization(2);
    vector <double> _polarization(2);
    vector <double> dielectric(2);
    double waterPol = residueTemplate::getPolarizability(52);
    double waterVol = residueTemplate::getVolume(52);
    double totalVol = residue::cutoffCubeVolume;
    polarization[0] = 0.0;
    polarization[1] = 0.0;

    // get volume and polarizabilities through protein around atom
     _polarization = _chain->calculateDielectric(_chain, _residue, _atom);
     polarization[0] += _polarization[0];
     polarization[1] += _polarization[1];

    // calculate local dielectric for atom
	double totalWaterVol = totalVol-(polarization[0]/1.89);
    int waters = totalWaterVol/waterVol;
    double totalWaterPol = waters*waterPol;
	double die = 1+4*3.14*((waters)/totalVol)*(totalWaterPol+polarization[1]);
    if (die < 2) { die = 2.0;}
    if (die > 78) { die = 78.0;}
    dielectric[0] = die;
    dielectric[1] = waters;
    return dielectric;
}

void protein::updateChainIndependentDielectrics(UInt _chainIndex)
{
	vector <double> dielectric(2);
	for(UInt i=0; i<itsChains[_chainIndex]->itsResidues.size(); i++)
	{
		for(UInt j=0; j<itsChains[_chainIndex]->itsResidues[i]->itsAtoms.size(); j++)
		{
			dielectric = this->calculateChainIndependentDielectric(itsChains[_chainIndex], itsChains[_chainIndex]->itsResidues[i], itsChains[_chainIndex]->itsResidues[i]->itsAtoms[j]);
			itsChains[_chainIndex]->itsResidues[i]->itsAtoms[j]->setDielectric(dielectric[0]);
			itsChains[_chainIndex]->itsResidues[i]->itsAtoms[j]->setNumberofWaters(dielectric[1]);

		}
	}
}

vector <double> protein::calculateResidueIndependentDielectric(residue* _residue, atom* _atom)
{
	vector <double> polarization(2);
	vector <double> _polarization(2);
	vector <double> dielectric(2);
	double waterPol = residueTemplate::getPolarizability(52);
	double waterVol = residueTemplate::getVolume(52);
	double totalVol = residue::cutoffCubeVolume;
	polarization[0] = 0.0;
	polarization[1] = 0.0;

	// get volume and polarizabilities through protein around atom
	 _polarization = _residue->calculateDielectric(_atom);
	 polarization[0] += _polarization[0];
	 polarization[1] += _polarization[1];

	// calculate local dielectric for atom
	double totalWaterVol = totalVol-(polarization[0]/1.89);
	int waters = totalWaterVol/waterVol;
	double totalWaterPol = waters*waterPol;
	double die = 1+4*3.14*((waters)/totalVol)*(totalWaterPol+polarization[1]);
	if (die < 2) { die = 2.0;}
	if (die > 78) { die = 78.0;}
	dielectric[0] = die;
	dielectric[1] = waters;
	return dielectric;
}

void protein::updateResidueIndependentDielectrics(UInt _chainIndex, UInt _resIndex)
{
	vector <double> dielectric(2);
	for(UInt j=0; j<itsChains[_chainIndex]->itsResidues[_resIndex]->itsAtoms.size(); j++)
	{
		dielectric = calculateResidueIndependentDielectric(itsChains[_chainIndex]->itsResidues[_resIndex], itsChains[_chainIndex]->itsResidues[_resIndex]->itsAtoms[j]);
		itsChains[_chainIndex]->itsResidues[_resIndex]->itsAtoms[j]->setDielectric(dielectric[0]);
		itsChains[_chainIndex]->itsResidues[_resIndex]->itsAtoms[j]->setNumberofWaters(dielectric[1]);
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

//Functions used in fast (non-redundant) energy calculation (protEnergy) //////////////////
void protein::updateDielectrics()
{
	vector <double> dielectric(2);
	for(UInt i=0; i<itsChains.size(); i++)
	{
		for(UInt j=0; j<itsChains[i]->itsResidues.size(); j++)
		{
			for(UInt k=0; k<itsChains[i]->itsResidues[j]->itsAtoms.size(); k++)
			{
				dielectric = this->calculateDielectric(itsChains[i], itsChains[i]->itsResidues[j], itsChains[i]->itsResidues[j]->itsAtoms[k]);
				itsChains[i]->itsResidues[j]->itsAtoms[k]->setDielectric(dielectric[0]);
				itsChains[i]->itsResidues[j]->itsAtoms[k]->setNumberofWaters(dielectric[1]);
			}
		}
	}
}

void protein::updatePositionDielectrics(UInt _chainIndex, UInt _residueIndex)
{
	vector <double> dielectric(2);
	for(UInt i=0; i<itsChains[_chainIndex]->itsResidues[_residueIndex]->itsAtoms.size(); i++)
	{
		dielectric = this->calculateDielectric(itsChains[_chainIndex], itsChains[_chainIndex]->itsResidues[_residueIndex], itsChains[_chainIndex]->itsResidues[_residueIndex]->itsAtoms[i]);
		itsChains[_chainIndex]->itsResidues[_residueIndex]->itsAtoms[i]->setDielectric(dielectric[0]);
		itsChains[_chainIndex]->itsResidues[_residueIndex]->itsAtoms[i]->setNumberofWaters(dielectric[1]);
	}
}

void protein::updateEnergyDatabase(vector < vector < vector <double> > > &_energies)
{
	if (_energies.empty()) // build energy database
	{
		buildResidueEnergyPairs(_energies);
	}
	else // update energy database
	{
		updateProtEnergy(_energies);
	}
}

double protein::protEnergy()
{
	updateEnergyDatabase(energies);

	double protEnergy = 0;
	for (UInt i = 0; i < energies.size(); i++) // total energy database
	{
		for (UInt j = 0; j < energies[i].size(); j++)
		{
			for (UInt k = 0; k < energies[i][j].size(); k++)
			{
				protEnergy += energies[i][j][k];
			}
		}
	}
	return protEnergy;
}

double protein::resEnergy(UInt chainIndex, UInt resIndex)
{
	if (itsChains[chainIndex]->itsResidues[resIndex]->getMoved() != 0)
	{
		updateEnergyDatabase(energies);
	}
	double resEnergy = 0;
	UInt chaini, chainj, resi, resj, k;
	for (chaini = 0; chaini < itsChains.size(); chaini++)
	{
		for (resi = 0; resi < itsChains[chaini]->itsResidues.size(); resi++)
		{
			k = 0;
			if (chaini == chainIndex && resi == resIndex)
			{
				resEnergy += energies[chaini][resi][k];
			}
			for (chainj = 0; chainj < chaini+1; chainj++)
			{
				for (resj = 0; resj < getNumResidues(chainj); resj++)
				{
					if (chaini == chainj && resi == resj)
					{
						break;
					}
					else
					{
						k++;
						if ((chaini == chainIndex && resi == resIndex) || (chainj == chainIndex && resj == resIndex))
						{
							resEnergy += energies[chaini][resi][k];
						}
					}
				}
			}
		}
	}
	return resEnergy;
}

double protein::getMedianResEnergy()
{
	double median, resE;
	vector <double> resEnergies;
	updateEnergyDatabase(energies);
	for (UInt i = 0; i < itsChains.size(); i++)
	{
		for (UInt j = 0; j < itsChains[i]->itsResidues.size(); j++)
		{
			resE = resEnergy(i,j);
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

double protein::getMedianResEnergy(UIntVec _activeChains)
{
	double median, resE;
	vector <double> resEnergies;
	updateEnergyDatabase(energies);
	for (UInt i = 0; i < _activeChains.size(); i++)
	{
		for (UInt j = 0; j < itsChains[_activeChains[i]]->itsResidues.size(); j++)
		{
			resE = resEnergy(_activeChains[i],j);
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

double protein::getMedianResEnergy(UIntVec _activeChains, UIntVec _activeResidues)
{
	double median, resE;
	vector <double> resEnergies;
	updateEnergyDatabase(energies);
	for (UInt i = 0; i < _activeChains.size(); i++)
	{
		for (UInt j = 0; j < _activeResidues.size(); j++)
		{
			resE = resEnergy(_activeChains[i], _activeResidues[j]);
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

double protein::getMedianDeltaH()
{
	double median, resE;
	vector <double> resEnergies;
	for (UInt i = 0; i < itsChains.size(); i++)
	{
		for (UInt j = 0; j < itsChains[i]->itsResidues.size(); j++)
		{
			resE = deltaH(i,j);
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


void protein::buildResidueEnergyPairs(vector < vector < vector <double> > > &_energies)
{
	_energies.clear();
	double Energy;
	vector < vector < vector <double> > > chainE;
	vector < vector <double> > resE;
	vector <double> E;

	//populate energy vector with starting energies
	if (residueTemplate::itsAmberElec.getScaleFactor() != 0.0 || residue::getHydroSolvationScaleFactor() != 0.0 || residue::getElectroSolvationScaleFactor() != 0.0)
	{
		updateDielectrics();
	}
	UInt chaini, chainj, resi, resj;
	for (chaini = 0; chaini < itsChains.size(); chaini++)
	{
		resE.clear();
		for (resi = 0; resi < getNumResidues(chaini); resi++)
		{
			E.clear();
			for (chainj = 0; chainj < chaini+1; chainj++)
			{
				for (resj = 0; resj < getNumResidues(chainj); resj++)
				{
					//cout << chaini << " " << resi << " " << chainj << " " << resj<< " | ";
					if (chaini == chainj && resi == resj)
					{
						Energy = itsChains[chaini]->itsResidues[resi]->intraSoluteEnergy();
						E.push_back(Energy);
						break;
					}
					else
					{
						Energy = itsChains[chaini]->itsResidues[resi]->interSoluteEnergy(itsChains[chainj]->itsResidues[resj]);
						E.push_back(Energy);
					}
				}
			}
			resE.push_back(E);
			//cout << endl;
		}
		chainE.push_back(resE);
	}
	for (chaini = 0; chaini < itsChains.size(); chaini++)
	{
		for (resi = 0; resi < itsChains[chaini]->itsResidues.size(); resi++)
		{
			itsChains[chaini]->itsResidues[resi]->setMoved(0);
		}
	}
	_energies = chainE;
	return;
}

void protein::updateProtEnergy(vector < vector < vector <double> > > &_energies)
{
	double residueEnergy;

	//update energies of residues transformed since last calculation
	UInt chaini, chainj, resi, resj, k;
	for (chaini = 0; chaini < itsChains.size(); chaini++)
	{
		for (resi = 0; resi < itsChains[chaini]->itsResidues.size(); resi++)
		{
			k = 0;
			for (chainj = 0; chainj < chaini+1; chainj++)
			{
				for (resj = 0; resj < getNumResidues(chainj); resj++)
				{
					if (chaini == chainj && resi == resj)
					{
						if (itsChains[chaini]->itsResidues[resi]->getMoved() != 0)
						{
							if (residueTemplate::itsAmberElec.getScaleFactor() != 0.0 || residue::getHydroSolvationScaleFactor() != 0.0 || residue::getElectroSolvationScaleFactor() != 0.0)
							{
								updatePositionDielectrics(chaini, resi);
							}
							residueEnergy = itsChains[chaini]->itsResidues[resi]->intraSoluteEnergy();
							_energies[chaini][resi][k] = residueEnergy;
						}
						k++;
						break;
					}
					else
					{
						if (itsChains[chaini]->itsResidues[resi]->getMoved() != 0 || itsChains[chainj]->itsResidues[resj]->getMoved() != 0)
						{
							if (itsChains[chainj]->itsResidues[resj]->getMoved() != 0 && (residueTemplate::itsAmberElec.getScaleFactor() != 0.0 || residue::getHydroSolvationScaleFactor() != 0.0 || residue::getElectroSolvationScaleFactor() != 0.0))
							{
								updatePositionDielectrics(chainj, resj);
							}
							residueEnergy = itsChains[chaini]->itsResidues[resi]->interSoluteEnergy(itsChains[chainj]->itsResidues[resj]);
							_energies[chaini][resi][k] = residueEnergy;
						}
					}
					k++;
				}
			}
		}
	}
	for (chaini = 0; chaini < itsChains.size(); chaini++)
	{
		for (resi = 0; resi < itsChains[chaini]->itsResidues.size(); resi++)
		{
			itsChains[chaini]->itsResidues[resi]->setMoved(0);
		}
	}
	return;
}

double protein::getFreeAminoAcidEnergy()
{
	double refEnergy = 0;
	UInt chaini, resi;
	for (chaini = 0; chaini < itsChains.size(); chaini++)
	{
		for (resi = 0; resi < itsChains[chaini]->itsResidues.size(); resi++)
		{
			refEnergy += getFreeAminoAcidEnergy(chaini,resi);
		}
	}
	return refEnergy;
}

double protein::getFreeAminoAcidEnergy(UInt _chainIndex, UInt _resIndex)   // protEnergy of free amino acids in current conformation
{
	updateResidueIndependentDielectrics(_chainIndex, _resIndex);
	double refEnergy = itsChains[_chainIndex]->itsResidues[_resIndex]->intraSoluteEnergy();
	itsChains[_chainIndex]->itsResidues[_resIndex]->setMoved(1);
	return refEnergy;
}

double protein::deltaH()
{
    double deltaH;
	deltaH = protEnergy()-getFreeAminoAcidEnergy();
    return deltaH;
}

double protein::deltaH(UInt chainIndex, UInt resIndex)
{
	double deltaH = resEnergy(chainIndex, resIndex)-getFreeAminoAcidEnergy(chainIndex, resIndex);
    return deltaH;
}

//-/////////////////////////////////////////////////////////////////////////////

double protein::BBEnergy()
{
	double energy = 0.0;
	for (UInt i = 0; i < itsChains.size(); i ++)
	{
		energy += itsChains[i]->BBEnergy();
	}
	return energy;
}

vector <double> protein::chainBindingEnergy()
{
    double bindingEnergy, complexEnergy, intraChainEnergy = 0.0;
    vector <double> Energy;
    complexEnergy = intraSoluteEnergy(true);
    Energy.push_back(complexEnergy);
    UInt numChains = this->getNumChains();
    for (UInt j = 0; j < numChains; j++)
    {
        this->updateChainIndependentDielectrics(j);
        intraChainEnergy += itsChains[j]->intraSoluteEnergy();
    }
    bindingEnergy = complexEnergy - intraChainEnergy;
    Energy.push_back(bindingEnergy);
    return Energy;
}

double protein::bindingPositionSoluteEnergy(UInt _chain, UInt _residue, UInt _otherChain)
{
	double bindingPositionE, intraChainPositionE, complexPositionE = 0.0;
	vector <int> position;
	position.resize(0);
	position.push_back(-1);
	position.push_back((int)_chain);
	position.push_back((int)_residue);
    this->updateChainIndependentDielectrics(_chain);
	intraChainPositionE = itsChains[_chain]->getPositionIntraSoluteEnergy(position);
    this->updateDielectrics();
	complexPositionE += itsChains[_chain]->getPositionIntraSoluteEnergy(position);
	complexPositionE += itsChains[_chain]->getPositionInterSoluteEnergy(position, itsChains[_otherChain]);
	bindingPositionE = complexPositionE - intraChainPositionE;
	return bindingPositionE;
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

void protein::protOpt(bool _backbone)
{   // Sidechain and backrub optimization with a local dielectric scaling of electrostatics and corresponding Born/Gill implicit solvation energy
	//_plateau: the number of consecutive optimization cycles without an energy decrease (default: 150 for general purpose optimization)

	//--Initialize variables for loop, calculate starting energy and build energy vectors---------------
	UInt randchain, randres, randrestype, randrot, chainNum = getNumChains(), keep, nobetter = 0, _plateau = 200;
	double deltaTheta = 0, Energy, resE, medResE, pastEnergy = protEnergy(), currentBetaChi, energyBuffer = 0.05;
	vector < vector <double> > currentRot; vector <UIntVec> allowedRots; srand (time(NULL));

	//--Run optimizaiton loop to relative minima, determined by _plateau----------------------------
	do
	{   //--choose random residue
		randchain = rand() % chainNum;
		randres = rand() % getNumResidues(randchain);
		randrestype = getTypeFromResNum(randchain, randres);
		nobetter++;

		//--Backrub optimization-----------------------------------------------------------------------
		if (nobetter > _plateau && _backbone)
		{
			resE = resEnergy(randchain, randres), medResE = getMedianResEnergy();
			if (resE > medResE)
			{
				//--randomly choose degree change (-1 or +1) of backrub dihedral (Ca-Cb angle)
				do
				{ deltaTheta = ((rand() % 3) -1);
				} while (deltaTheta == 0);

				//--transform angle while energy improves, until energy degrades, then revert one step
				do
				{
					keep = 0;
					currentBetaChi = getBetaChi(randchain, randres);
					setBetaChi(randchain, randres, deltaTheta+currentBetaChi);
					Energy = protEnergy();
					if (Energy < (pastEnergy-energyBuffer))
					{
						nobetter = 0, keep = 1, pastEnergy = Energy;
					}
				} while (keep == 1);
				setBetaChi(randchain, randres, currentBetaChi);
			}
		}

		//--Rotamer optimization-----------------------------------------------------------------------
		resE = resEnergy(randchain, randres), medResE = getMedianResEnergy();
		if (resE > medResE)
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
					Energy = protEnergy();
					if (Energy < (pastEnergy-energyBuffer))
					{
						nobetter = 0, pastEnergy = Energy; break;
					}
					else
					{
						setSidechainDihedralAngles(randchain, randres, currentRot);
					}
				}
			}
		}
	} while (nobetter < _plateau * 1.2);
	return;
}

void protein::protOpt(bool _backbone, UIntVec _frozenResidues, UIntVec _activeChains) //_activeChain only optimized and energy calculated for it
{   // Sidechain and backrub optimization with a local dielectric scaling of electrostatics and corresponding Born/Gill implicit solvation energy
	//_plateau: the number of consecutive optimization cycles without an energy decrease (default: 150 for general purpose optimization)

	//--Initialize variables for loop, calculate starting energy and build energy vectors---------------
	UInt randchain, randres, randrestype, randrot, chainNum = _activeChains.size(), keep, nobetter = 0, _plateau = 150;
	double deltaTheta = 0, Energy, resE, medResE, pastEnergy = protEnergy(), currentBetaChi, energyBuffer = 0.05;
	vector < vector <double> > currentRot; vector <UIntVec> allowedRots; srand (time(NULL));
	bool skip;

	//--Run optimizaiton loop to relative minima, determined by _plateau----------------------------
	do
	{   //--choose random residue not frozen of active chains
		randchain = _activeChains[rand() % chainNum];
		do
		{	skip = false;
			randres = rand() % getNumResidues(randchain);
			for (UInt i = 0; i < _frozenResidues.size(); i++)
			{
				if (randres == _frozenResidues[i]) {skip = true;}
			}
		} while (skip);
		randrestype = getTypeFromResNum(randchain, randres);
		nobetter++;

		//--Backrub optimization-----------------------------------------------------------------------
		if (nobetter > _plateau && _backbone)
		{
			resE = resEnergy(randchain, randres), medResE = getMedianResEnergy(_activeChains);
			if (resE > medResE)
			{
				//--randomly choose degree change (-1 or +1) of backrub dihedral (Ca-Cb angle)
				do
				{ deltaTheta = ((rand() % 3) -1);
				} while (deltaTheta == 0);

				//--transform angle while energy improves, until energy degrades, then revert one step
				do
				{
					keep = 0;
					currentBetaChi = getBetaChi(randchain, randres);
					setBetaChi(randchain, randres, deltaTheta+currentBetaChi);
					Energy = protEnergy();
					if (Energy < (pastEnergy-energyBuffer))
					{
						nobetter = 0, keep = 1, pastEnergy = Energy;
					}
				} while (keep == 1);
				setBetaChi(randchain, randres, currentBetaChi);
			}
		}

		//--Rotamer optimization-----------------------------------------------------------------------
		resE = resEnergy(randchain, randres), medResE = getMedianResEnergy(_activeChains);
		if (resE > medResE)
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
					Energy = protEnergy();
					if (Energy < (pastEnergy-energyBuffer))
					{
						nobetter = 0, pastEnergy = Energy; break;
					}
					else
					{
						setSidechainDihedralAngles(randchain, randres, currentRot);
					}
				}
			}
		}
	} while (nobetter < _plateau * 1.2);
	return;
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

void protein::optimizeRotamersPN()
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
        rotamerArray = rotamerDEE();
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

                cout << "Optimized!!!" << endl;
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
    double PI = 3.1415;

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

UInt protein::getNumHardClashes()
{
	UInt numClashes = 0;
	for (UInt i = 0; i < itsChains.size(); i ++)
	{
		numClashes += itsChains[i]->getNumHardClashes();
		for (UInt j = i +1; j < itsChains.size(); j ++)
		{
			numClashes += itsChains[i]->getNumHardClashes(itsChains[j]);
		}
	}
	return numClashes;
}
