#include "chain.h"

typedef vector<chainPosition*>::iterator iterCP;
typedef vector<UInt>::iterator iterUINT;
typedef vector<int>::iterator iterINT;
typedef vector<vector<int> >::iterator iterINTVEC;
UInt chain::howMany = 0;

void chain::initialize()
{	itsResidues.resize(0);
	itsChainPositions.resize(0);
	itsSecondaryStructures.resize(0);
	itsRepackActivePositionMap.resize(0);
	itsIndependentPositions.resize(0);
	itsResidueLinkageMap.resize(0);
	itsSpaceLink.newsize(3);
	itsSpinLink.newsize(3);
	itsLastTargetResidue = -1;
	for (UInt i = 0; i < 3; i++)
	{
		itsSpaceLink[i] = 1.0; // set default space link of 1,1,1
		itsSpinLink[i] = 1.0;  // set default spin link of 1,1,1
	}
	chainModBuffer undo;
	chainModBuffer redo;
	itsBuffers.push_back(undo);
	itsBuffers.push_back(redo);
	resetUndoBuffer();
	resetRedoBuffer();
	stateBuffers.resize(0);
	howMany++;
}

chain::chain()
{
	initialize();
	itsChainID = '0';
}

chain::chain(char _id)
{
	initialize();
	itsChainID = _id;
}

// deep copy constructor
chain::chain(const chain& _rhs)
{
	//cout << "chain deep copy constructor called." << endl;
	initialize();
	for (UInt i=0;i<_rhs.itsResidues.size();i++)
	{	addResidue(new residue( *(_rhs.itsResidues[i])));
	}
	//cout << "residues added." << endl;
	for (UInt i=0;i<_rhs.itsChainPositions.size();i++)
	{	if (_rhs.itsChainPositions[i])
		{	addChainPosition(new chainPosition( *(_rhs.itsChainPositions[i])));
		}
		else
		{	addChainPosition(0);
		}
	}
	//cout << "chainpositions added." << endl;
	itsRepackActivePositionMap = _rhs.itsRepackActivePositionMap;
    itsIndependentPositions = _rhs.itsIndependentPositions;
    itsResidueLinkageMap = _rhs.itsResidueLinkageMap;
	itsSpaceLink = _rhs.itsSpaceLink;

	for (UInt i=0;i<_rhs.itsSecondaryStructures.size();i++)
	{	if (_rhs.itsSecondaryStructures[i])
		{
			addSecondaryStructure(new secondaryStructure( *(_rhs.itsSecondaryStructures[i])));
		}
		else
		{	addSecondaryStructure(0);
		}
	}
	//cout << "super secondary structures added." << endl;
	itsChainID = _rhs.itsChainID;
	//buffers
	itsBuffers[0] = _rhs.itsBuffers[0];
	itsBuffers[1] = _rhs.itsBuffers[1];
	stateBuffers = _rhs.itsBuffers;
}

chain::~chain()
{
	for (UInt i=0; i<itsResidues.size(); i++)
	{	delete itsResidues[i];
	}
	for (UInt i=0; i<itsChainPositions.size(); i++)
	{	if (itsChainPositions[i])
			delete itsChainPositions[i];
	}
	for (UInt i=0; i<itsSecondaryStructures.size(); i++)
	{	if (itsSecondaryStructures[i])
			delete itsSecondaryStructures[i];
	}
    howMany--;
}

void chain::removeResidue(UInt _resNum)
{
    delete itsResidues[_resNum];
    itsResidues.resize(itsResidues.size()-1);
    delete itsChainPositions[_resNum];
    itsChainPositions.resize(itsChainPositions.size()-1);
    delete itsSecondaryStructures[_resNum];
    itsSecondaryStructures.resize(itsSecondaryStructures.size()-1);
}

void chain::add(residue* _pResidue)
{
	addResidue(_pResidue);
	// A Null chainPosition is a flag for an residue which has not
	// yet been "activated" for mutation/modification
	addChainPosition(0);
	secondaryStructure* pSecStruc = new secondaryStructure(itsResidues.size()-1);
	addSecondaryStructure(pSecStruc);
	UInt index = itsResidues.size() -1;
    itsIndependentPositions.push_back(index);

    vector <int> tempIntVec;
    tempIntVec.resize(0);
    tempIntVec.push_back(-1);
    itsResidueLinkageMap.push_back(tempIntVec);
}

void chain::activateChainPosition(UInt _indexInChain)
{	if (_indexInChain < itsResidues.size())
	{
		chainPosition* newCP = new chainPosition(_indexInChain,itsResidues[_indexInChain]->getTypeIndex());
		//this sets the backbone-dependent rotamer libraries
		newCP->setRotamerLibIndex(itsSecondaryStructures[_indexInChain]->getSecondaryStructure());
		//cout << "the secondary structure index for " << _indexInChain << " is " << itsSecondaryStructures[_indexInChain]->getSecondaryStructure() << endl;
		itsChainPositions[_indexInChain] = newCP;
		//now let's make sure that only alpha amino acids are allowed for an alpha
		//only beta for a beta, etc.
        //UInt resTypeIndex = itsResidues[_indexInChain]->getTypeIndex();
        //UInt numBpt = residue::getNumBpt(resTypeIndex);
		//cout << "Num branchpoints in actual residue ==" << numBpt << endl;
        //UInt initialNumAllowedRes = newCP->getNumAllowedRes();
		vector<UInt> notAllowedSet;
		notAllowedSet.resize(0);
        //for (UInt i=0; i< initialNumAllowedRes; i++) // this prevents using branchpoints as dihedral pivot points for complex residues, so commented out.  doug p. 2015
        //{
            //UInt identity = newCP->itsAllowedResidues[i].getIdentity();
            //UInt numBptInPossibleRes = residue::getNumBpt(identity);
			//cout << "numBptInPossibleRes = " << numBptInPossibleRes << endl;
            //if ( numBptInPossibleRes != numBpt )
            //{
            //	notAllowedSet.push_back(identity);
            //}

        //}
		/*for (UInt i=0; i<notAllowedSet.size(); i++)
		{
			cout << "noAllowedSet[" << i << "] = " << notAllowedSet[i] << endl;
		}
		*/
		for (UInt i=0; i<notAllowedSet.size(); i++)
		{
			setResNotAllowed(_indexInChain,notAllowedSet[i]);
		}
	}
}

void chain::activateAllChainPositions()
{	for(UInt i=0; i<itsResidues.size(); i++)
	{	activateChainPosition(i);
	}
}

void chain::makeAtomSilent(const UInt _resIndex, const UInt _atomIndex)
{
	if (_resIndex >=0 && _resIndex < itsResidues.size())
	{
		itsResidues[_resIndex]->makeAtomSilent(_atomIndex);
	}
	else
		cout << "ERROR in chain::makeAtomSilent ...\n\t" << _resIndex << " is not a valid residue number." << endl;
	return;
}

void chain::makeResidueSilent(const UInt _resIndex)
{
	if (_resIndex >=0 && _resIndex < itsResidues.size())
	{
		itsResidues[_resIndex]->makeResidueSilent();
	}
	else
		cout << "ERROR in chain::makeAtomSilent ...\n\t" << _resIndex << " is not a valid residue number." << endl;
	return;
}

bool chain::activateForRepacking(UInt _indexInChain)
{
	if (_indexInChain < itsResidues.size())
	{	
		for (UInt i=0; i<itsRepackActivePositionMap.size(); i++)
		{	if ( itsRepackActivePositionMap[i] == _indexInChain)
			{
				return true;
			}
		}

	}
	return false;
}
bool chain::activateForRepacking(UInt _indexStart, UInt _indexEnd)
{	for(UInt i=_indexStart; i<_indexEnd; i++)
	{	if(!activateForRepacking(i))
		{	return false;
		}
	}
	return true;
}

bool chain::activateAllForRepacking()
{	for(UInt i=0; i<itsResidues.size(); i++)
	{	if(!activateForRepacking(i))
		{	return false;
		}
	}
	return true;
}

void chain::randomizeSystem(ran& _ran)
{	for (UInt i=0; i<itsRepackActivePositionMap.size(); i++)
	{	UInt pos = itsRepackActivePositionMap[i];
		int id = chooseTargetIdentity(pos,_ran);
		if (id >=0)
		{	mutate(pos,id);
			int bpt = chooseTargetBranchpoint(pos,_ran);
			if (bpt >=0)
			{	int rot = chooseTargetRotamer(pos,bpt,_ran);
				if (rot >=0)
				{	setRotamerWithoutBuffering(pos,bpt,rot);
				}
			}
		}
		resetUndoBuffer();
	}
}

void chain::makeAllAlanine()
{	for (UInt i=0; i<itsRepackActivePositionMap.size(); i++)
	{	UInt pos = itsRepackActivePositionMap[i];
		mutateWithoutBuffering(pos,0);
	}
	return;
}

void chain::finishChainBuild()
{
	assignSecondaryStructure();
	activateAllChainPositions();
}

void chain::assignSecondaryStructure()
{
	for (unsigned int i=0; i< itsResidues.size(); i++)
	{
		itsSecondaryStructures[i]->assign(this);
	}
}

void chain::listSecondaryStructure() const
{
	for (unsigned int i=0; i<itsSecondaryStructures.size(); i++)
	{
		if (itsSecondaryStructures[i])
		{
			cout << itsResidues[i]->getResNum() << "  ";
			cout << itsSecondaryStructures[i]->getSecondaryStructure();
			cout << endl;
		}
		else
		{
			cout << "ERROR: secondary structure for residue ";
			cout << itsResidues[i]->getResNum() << " undefined";
			cout << endl;
		}
	}
}

void chain::listDihedrals()
{
	vector< vector< double> > dihedrals;
	for (UInt i=0; i<itsResidues.size(); i++)
	{	cout <<"dihedrals: res " << i << "  ";
		dihedrals = getSidechainDihedralAngles(i);
		for (UInt j=0; j<dihedrals.size(); j++)
		{
			for (UInt k=0; k<dihedrals[j].size();k++)
			{
				cout << dihedrals[j][k] << " ";
			}
			cout << endl;
		}
		cout << endl;
	}
}

void chain::listChiDefinitions() const
{
}

void chain::mutate(const UInt _indexInChain, const UInt _aaType)
{   // Check if the residue is even in our chain
    if (_indexInChain >= 0 && itsChainPositions[_indexInChain])
    {   // Now check if this position is allowed to be mutated
        // to this residue type
        if (itsChainPositions[_indexInChain]->residueIsAllowed(_aaType))
        {
            residue* pOldRes = itsResidues[_indexInChain];
            //buffer everything
            bufferResidueIntoUndoBuffer(_indexInChain);
            // mutate to new residue and place new pointer in chain
            itsResidues[_indexInChain] = pOldRes->mutate( _aaType );
            delete pOldRes;
            itsChainPositions[_indexInChain]->setCurrentResIndex(_aaType);

            int numBpt = residue::getNumBpt(_aaType);
            vector<UInt> rotState;
            rotState.resize(0);
            if (numBpt >= 0)
            {
                for (UInt i=0; i < (UInt)numBpt; i++)
                {
                    UInt numChis = getNumChis(_indexInChain,i);
                    if (numChis > 0)
                    {
                        UIntVec allowedRotamers = getAllowedRotamers(_indexInChain, _aaType, i);
                        UInt rotIndex = allowedRotamers[0];
                        if (rotIndex >=0)
                        {
                            //cout << rotIndex << " ";
                            setRotamerWithoutBuffering(_indexInChain,i,rotIndex);
                        }
                    }
                }
            }
            return;
        }
        else
        {   cout << "Mutation not allowed..." << endl;
            cout << endl << "ABORTED" << endl;
            return;
        }
    }
    else
    {   cout << "Residue not in active set : " << _indexInChain  << endl;
        cout << endl << "ABORTED" << endl;
        return;
    }
}


void chain::mutateWithoutBuffering(const UInt _indexInChain, const UInt _aaType)
{	// Check if the residue is even in our chain
	if (_indexInChain >= 0 && itsChainPositions[_indexInChain])
	{	// Now check if this position is allowed to be mutated
		// to this residue type
		if (itsChainPositions[_indexInChain]->residueIsAllowed(_aaType))
		{	residue* pOldRes = itsResidues[_indexInChain];
			// mutate to new residue and place new pointer in chain
			itsResidues[_indexInChain] = pOldRes->mutate( _aaType );
			delete pOldRes;
			itsChainPositions[_indexInChain]->setCurrentResIndex(_aaType);
			int numBpt = residue::getNumBpt(_aaType);
			vector<UInt> rotState;
			rotState.resize(0);
			if (numBpt >= 0)
			{
				for (UInt i=0; i < (UInt)numBpt; i++)
				{
					UInt numChis = getNumChis(_indexInChain,i);
					if (numChis > 0)
					{
						UIntVec allowedRotamers = getAllowedRotamers(_indexInChain, _aaType, i);
						UInt rotIndex = allowedRotamers[0];
						if (rotIndex >=0)
						{
							setRotamerWithoutBuffering(_indexInChain,i,rotIndex);
						}
					}
				}
			}
			return;
		}
		else
		{	cout << "Mutation not allowed..." << endl;
			cout << endl << "ABORTED" << endl;
			return;
		}
	}
	else
	{	cout << "Residue not in active set : " << _indexInChain  << endl;
		cout << endl << "ABORTED" << endl;
		return;
	}
}

void chain::fixBrokenResidue(const UInt _indexInChain)
{	if (_indexInChain < itsResidues.size())
	{	// Now check if this position is allowed to be mutated
		// to this residue type

		residue* pOldRes = itsResidues[_indexInChain];
        vector < vector <double> > currentRot = getSidechainDihedralAngles(_indexInChain);
		itsResidues[_indexInChain] = pOldRes->mutate( itsResidues[_indexInChain]->getTypeIndex() );
        setSidechainDihedralAngles(_indexInChain, currentRot);
		itsResidues[_indexInChain]->isArtificiallyBuilt = true;
		delete pOldRes;
	}
}

void chain::rebuildResidue(const UInt _indexInChain)
{	if (_indexInChain < itsResidues.size())
	{
		vector < vector <double> > currentRot;
		residue* pOldRes = itsResidues[_indexInChain];
		UInt itsResType = itsResidues[_indexInChain]->getTypeIndex();
		currentRot = getSidechainDihedralAngles(_indexInChain);
		bool fliptoD = isDAminoAcid(pOldRes);
		if (fliptoD){itsResidues[_indexInChain] = pOldRes->mutateNew(itsResType+Daa);}
		else{itsResidues[_indexInChain] = pOldRes->mutateNew(itsResType);}
		setSidechainDihedralAngles(_indexInChain, currentRot);
		itsResidues[_indexInChain]->isArtificiallyBuilt = true;
		delete pOldRes;
	}
}

void chain::rebuildResiduesInChain()
{
	for (UInt i=0; i< itsResidues.size(); i++)
	{
		rebuildResidue(i);
	}
}

bool chain::isDAminoAcid(residue* currentRes)
{
	bool flip = false;
	UInt resType = currentRes->getTypeIndex();
	if (currentRes->isL(resType) && currentRes->getNumAtoms() > 4 ){
			
		dblVec Ncoords = currentRes->getCoords(0);
		dblVec Ccoords = currentRes->getCoords(2);
		dblVec CAcoords = currentRes->getCoords(1);
		dblVec CBcoords; if (resType == 19){CBcoords = currentRes->getCoords(6);} //proline
		else{CBcoords = currentRes->getCoords(4);}
		
		dblVec N_CA = CAcoords - Ncoords;
		dblVec N_C = Ccoords - Ncoords;
		dblVec CA_CB = CBcoords - CAcoords;

		// Calculate Cross Products
		dblVec N_CA_X_N_C(3);
		N_CA_X_N_C = CMath::cross(N_CA,N_C);

		double chiralityAngle = CMath::dotProduct(N_CA_X_N_C, CA_CB);
		if (chiralityAngle > 0){flip = true;}
	}
	return flip;
}

vector<chainModBuffer> chain::performRandomMutation(ran& _ran)
{
	if (itsRepackActivePositionMap.size() == 0)
	{	// if this fails, there are no active residues!
		resetUndoBuffer();
		resetRedoBuffer();
		cout << endl << "ABORTED(level 1rm)" << endl;
		return itsBuffers;
	}
	int residueToModify = chooseTargetResidue(_ran);
	if (residueToModify == -1)
	{	// if this test fails, there are no active residues!
		resetUndoBuffer();
		resetRedoBuffer();
		cout << endl << "ABORTED(level 2rm)" << endl;
		return itsBuffers;
	}
	cout << "Res #" << residueToModify << " ";
	int identity = chooseTargetIdentity(residueToModify, _ran);
	if (identity == -1)
	{	// if this test fails, you may still be able to try another residue
		resetUndoBuffer();
		resetRedoBuffer();
		cout << endl << "ABORTED(level 3rm)" << endl;
		return itsBuffers;
	}
	// mutate function will buffer current rotamer state
	mutate(residueToModify,identity);
	cout << "To " << identity << " ";
	// now we need to see if this needs to have its rotamer state set
	int numBpt = residue::getNumBpt(identity);
	vector<UInt> rotState;
	rotState.resize(0);
	if (numBpt >= 0)
	{
		for (UInt i=0; i < (UInt)numBpt; i++)
		{
			int rotIndex = chooseTargetRotamer(residueToModify,i,_ran);
			if (rotIndex >=0)
			{
				cout << rotIndex << " ";
				setRotamerWithoutBuffering(residueToModify,i,rotIndex);
				rotState.push_back(rotIndex);
			}
		}
	}
	bufferResidueIntoRedoBuffer(residueToModify);
/*
	itsBuffers[1].setIndexInChainBuffer(residueToModify);
	itsBuffers[1].setResidueIdentityBuffer(identity);
	itsBuffers[1].setRotamerIndexBuffer(rotState);
	vector<vector<double> > tempAngles = itsResidues[residueToModify]->getSidechainDihedralAngles();
	itsBuffers[1].setSidechainDihedralAngleBuffer(tempAngles);
	*/
	return itsBuffers;
}

vector<chainModBuffer> chain::performRandomMutation(ran& _ran, vector <int> _position)
{
    if (itsRepackActivePositionMap.size() == 0)
    {   // if this fails, there are no active residues!
        resetUndoBuffer();
        resetRedoBuffer();
        cout << endl << "ABORTED(level 1rm)" << endl;
        return itsBuffers;
    }
    int residueToModify = _position[2];
	itsLastTargetResidue = residueToModify;
    if (residueToModify == -1)
    {   // if this test fails, there are no active residues!
        resetUndoBuffer();
        resetRedoBuffer();
        cout << endl << "ABORTED(level 2rm)" << endl;
        return itsBuffers;
    }
    cout << "Res #" << residueToModify << " ";
    int identity = chooseTargetIdentity(residueToModify, _ran);
    if (identity == -1)
    {   // if this test fails, you may still be able to try another residue
        resetUndoBuffer();
        resetRedoBuffer();
        cout << endl << "ABORTED(level 3rm)" << endl;
        return itsBuffers;
    }
    // mutate function will buffer current rotamer state
    mutate(residueToModify,identity);
    cout << "To " << identity << " ";
    // now we need to see if this needs to have its rotamer state set
    int numBpt = residue::getNumBpt(identity);
    vector<UInt> rotState;
    rotState.resize(0);
    if (numBpt >= 0)
    {
        for (UInt i=0; i < (UInt)numBpt; i++)
        {
            int rotIndex = chooseTargetRotamer(residueToModify,i,_ran);
            if (rotIndex >=0)
            {
                cout << rotIndex << " ";
                setRotamerWithoutBuffering(residueToModify,i,rotIndex);
                rotState.push_back(rotIndex);
            }
        }
    }
    bufferResidueIntoRedoBuffer(residueToModify);
/*
    itsBuffers[1].setIndexInChainBuffer(residueToModify);
    itsBuffers[1].setResidueIdentityBuffer(identity);
    itsBuffers[1].setRotamerIndexBuffer(rotState);
    vector<vector<double> > tempAngles = itsResidues[residueToModify]->getSidechainDihedralAngles();
    itsBuffers[1].setSidechainDihedralAngleBuffer(tempAngles);
    */
    return itsBuffers;
}




void chain::commitLastMutation()
{
	resetUndoBuffer();
	resetRedoBuffer();
#ifdef _CHAIN_DEBUG
	listDihedrals();
#endif
}

void chain::repeatModification(const chainModBuffer& _externalRedoBuffer)
{
	//cout << "chain::repeatModification(chainModBuffer _redoBuffer) called" << endl;
	//cout << "this function not yet implemented" << endl;
	resetUndoBuffer();
	resetRedoBuffer();

	int index = _externalRedoBuffer.getIndexInChainBuffer();
	//cout << "externalRedoBufferIndex = " << index << endl;

	if (index != -1)
	{
		// buffer the residue that's there into the undo buffer
		bufferResidueIntoUndoBuffer(index);
		UInt bufferResType = _externalRedoBuffer.getResidueIdentityBuffer();
		// compare this to what is there already
		if (bufferResType != itsResidues[index]->getTypeIndex())
		{	// we need to perform a mutation
			//cout << "repeating MUT" << endl;
			mutateWithoutBuffering(index,bufferResType);
		}
		vector<UInt> bufferRotamerIndex = _externalRedoBuffer.getRotamerIndexBuffer();
		if (bufferRotamerIndex != itsChainPositions[index]->getCurrentRotamerIndex())
		{
			// we need to set the rotameric state
			//cout << "repeating ROTC" << endl;
			setRotamerWithoutBuffering(index,bufferRotamerIndex);
		}

		vector<vector<double> > bufferDihedrals = _externalRedoBuffer.getSidechainDihedralAngleBuffer();
		if (bufferDihedrals != itsResidues[index]->getSidechainDihedralAngles())
		{
			//cout << "repeating DIHEDRAL CHANGE" << endl;
			for (UInt i=0; i<bufferDihedrals.size(); i++)
			{	for (UInt j=0; j<bufferDihedrals.size(); j++)
				{
					setDihedralWithoutBuffering(index,i,j,bufferDihedrals[i][j]);
				}
			}
		}
	}
}

void chain::undoLastMutation()
{	int index = itsBuffers[0].getIndexInChainBuffer();
	if (index < 0)
	{	cout << "Cannot undo mutation: indexInChainBuffer lost" << endl;
		itsBuffers[0].printAll();
		resetUndoBuffer();
		resetRedoBuffer();
		return;
	}

	int resID = itsBuffers[0].getResidueIdentityBuffer();
	if (resID <0)
	{	cout << "Cannot undo: residueIdentityBuffer lost" << endl;
		itsBuffers[0].printAll();
		resetUndoBuffer();
		resetRedoBuffer();
		return;
	}

	mutate(index,resID);

	itsChainPositions[index]->setCurrentResIndex(itsResidues[index]->getTypeIndex());
	vector< UInt> rotIndex = itsBuffers[0].getRotamerIndexBuffer();
	setRotamerWithoutBuffering(index,rotIndex);

	vector<vector<double> > angles = itsBuffers[0].getSidechainDihedralAngleBuffer();
	UInt bpt = angles.size();
	if (bpt != 0)
	{	for (UInt i=0; i<bpt; i++)
		{	UInt numAngles = angles[i].size();
			{	for (UInt j=0;j<numAngles;j++)
				{	double angle = angles[i][j];
					itsResidues[index]->setChi(i,j,angle);
				}
			}
		}
	}
	resetUndoBuffer();
	resetRedoBuffer();

#ifdef _CHAIN_DEBUG
	listDihedrals();
#endif
	return;
}

void chain::setRotamerWithoutBuffering(const UInt _indexInChain, vector<UInt> _rotamerIndex)
{	if (_indexInChain >=0 && _indexInChain < itsChainPositions.size())
	{	if (itsChainPositions[_indexInChain])
		{	UInt lib = itsChainPositions[_indexInChain]->getRotamerLibIndex();
			for (UInt i=0;i<_rotamerIndex.size(); i++)
			{	itsResidues[_indexInChain]->setRotamerWithCheck(lib,i, _rotamerIndex[i]);
				itsChainPositions[_indexInChain]->setCurrentRotamerIndex(_rotamerIndex[i]);
			}
		}
	}
}

void chain::setRotamerWithoutBuffering(const UInt _indexInChain, vector<UInt> _rotamerIndex, vector< vector< double> > dihedralAngles)
{	if (_indexInChain >=0 && _indexInChain < itsChainPositions.size())
	{	if (itsChainPositions[_indexInChain])
		{	UInt lib = itsChainPositions[_indexInChain]->getRotamerLibIndex();
			for (UInt i=0;i<_rotamerIndex.size(); i++)
			{	itsResidues[_indexInChain]->setRotamerWithCheck(lib,i, _rotamerIndex[i]);
				itsChainPositions[_indexInChain]->setCurrentRotamerIndex(_rotamerIndex[i]);
			}
		}
	}
}

void chain::setRotamerWithoutBuffering(const UInt _indexInChain, const UInt _bpt, const UInt _rotamerIndex)
{	if (_indexInChain >=0 && _indexInChain < itsChainPositions.size())
	{	if (itsChainPositions[_indexInChain])
		{	UInt lib = itsChainPositions[_indexInChain]->getRotamerLibIndex();
			itsResidues[_indexInChain]->setRotamer(lib,_bpt, _rotamerIndex);
			itsChainPositions[_indexInChain]->setCurrentRotamerIndex(_rotamerIndex);
		}
	}
}

void chain::setRotamerWithBuffering(const UInt _indexInChain, const UInt _bpt, const UInt _rotamerIndex)
{	if (_indexInChain >=0 && _indexInChain < itsChainPositions.size())
	{	if (itsChainPositions[_indexInChain])
		{	UInt lib = itsChainPositions[_indexInChain]->getRotamerLibIndex();
			bufferResidueIntoUndoBuffer(_indexInChain);
			itsResidues[_indexInChain]->setRotamerWithCheck(lib,_bpt, _rotamerIndex);
			itsChainPositions[_indexInChain]->setCurrentRotamerIndex(_rotamerIndex);
			bufferResidueIntoRedoBuffer(_indexInChain);
		}
	}
}

void chain::setRotamerWithBuffering(const UInt _indexInChain, vector<UInt> _rotamer)
{	if (_indexInChain >=0 && _indexInChain < itsChainPositions.size())
	{	if (itsChainPositions[_indexInChain])
		{	UInt lib = itsChainPositions[_indexInChain]->getRotamerLibIndex();
			bufferResidueIntoUndoBuffer(_indexInChain);
			for (UInt i=0; i< _rotamer.size(); i++)
			{
				itsResidues[_indexInChain]->setRotamerWithCheck(lib,i, _rotamer[i]);
			}
			itsChainPositions[_indexInChain]->setCurrentRotamerIndex(_rotamer);
			bufferResidueIntoRedoBuffer(_indexInChain);
		}
	}
}

void chain::setPolarHRotamerWithoutBuffering(const UInt _indexInChain, const UInt _rotamerIndex)
{	if (_indexInChain >=0 && _indexInChain < itsChainPositions.size())
	{	if (itsChainPositions[_indexInChain])
		{
			if(itsResidues[_indexInChain]->getHasPolarHRotamers())
			{
				itsResidues[_indexInChain]->setPolarHRotamerWithCheck(_rotamerIndex);
				itsChainPositions[_indexInChain]->setCurrentPolarHRotamerIndex(_rotamerIndex);
			}
		}
	}
}

void chain::setAllHydrogensOn(const bool _hydrogensOn)
{
	for(UInt i=0;i<itsResidues.size();i++)
	{
		itsResidues[i]->setHydrogensOn(_hydrogensOn);
	}
}

void chain::setAllPolarHydrogensOn(const bool _polarHydrogensOn)
{
	for(UInt i=0;i<itsResidues.size();i++)
	{
		itsResidues[i]->setPolarHydrogensOn(_polarHydrogensOn);
	}
}

double chain::netCharge()
{
    double nCharge = 0.0;
    for(UInt i=0;i<itsResidues.size();i++)
    {
        nCharge += itsResidues[i]->netCharge();
    }
    return nCharge;
}

vector<chainModBuffer> chain::performRandomRotamerChange(ran& _ran)
{
	if (itsRepackActivePositionMap.size() == 0)
	{
		 resetUndoBuffer();
		 resetRedoBuffer();
		 cout << endl << "ABORTED chain class level 1" << endl;
		 return itsBuffers;
	}
	int residueToModify = chooseTargetResidue(_ran);
	cout << residueToModify << " ";
	if (residueToModify == -1)
	{
		 resetUndoBuffer();
		 resetRedoBuffer();
		 cout << endl << "ABORTED chain class level 2" << endl;
		 return itsBuffers;
	}
	cout << itsResidues[residueToModify]->getType() << " ";
	//cout << itsChainPositions[residueToModify]->getCurrentResIndex() << "!" << endl;
	int branchpoint = chooseTargetBranchpoint(residueToModify,_ran);
	cout << branchpoint << " ";
	if (branchpoint < 0)
	{
		 resetUndoBuffer();
		 resetRedoBuffer();
		 cout << endl << "ABORTED chain class level 3" << endl;
		 return itsBuffers;
	}
	int rotIndex = chooseTargetRotamer(residueToModify,branchpoint,_ran);
	cout << rotIndex << " ";
	if (rotIndex < 0)
	{
		 resetUndoBuffer();
		 resetRedoBuffer();
		 cout << endl << "ABORTED" << endl;
		 return itsBuffers;
	}
	// setRotamerWithBuffering takes care of buffering residue
	setRotamerWithBuffering(residueToModify,branchpoint,rotIndex);
	return itsBuffers;
}

void chain::commitLastRotamerChange()
{	resetUndoBuffer();
	resetRedoBuffer();
#ifdef _CHAIN_DEBUG
	listDihedrals();
#endif
}

vector<chainModBuffer> chain::performRandomRotamerChange(ran& _ran, vector <int> _position)
{
    if (itsRepackActivePositionMap.size() == 0)
    {
         resetUndoBuffer();
         resetRedoBuffer();
         cout << endl << "ABORTED in chain class level 1" << endl;
         return itsBuffers;
    }
    int residueToModify = _position[2];
	itsLastTargetResidue = residueToModify;
    cout << residueToModify << " ";
    if (residueToModify == -1)
    {
         resetUndoBuffer();
         resetRedoBuffer();
         cout << endl << "ABORTED in chain class level 2" << endl;
         return itsBuffers;
    }
    cout << itsResidues[residueToModify]->getType() << " ";
    //cout << itsChainPositions[residueToModify]->getCurrentResIndex() << "!" << endl;
    int branchpoint = chooseTargetBranchpoint(residueToModify,_ran);
    cout << branchpoint << " ";
    if (branchpoint < 0)
    {
         resetUndoBuffer();
         resetRedoBuffer();
         cout << endl << "ABORTED in chain class level 3" << endl;
         return itsBuffers;
    }
    int rotIndex = chooseTargetRotamer(residueToModify,branchpoint,_ran);
    cout << rotIndex << " ";
    if (rotIndex < 0)
    {
         resetUndoBuffer();
         resetRedoBuffer();
         cout << endl << "ABORTED in chain class level 4" << endl;
         return itsBuffers;
    }
    // setRotamerWithBuffering takes care of buffering residue
    setRotamerWithBuffering(residueToModify,branchpoint,rotIndex);
    return itsBuffers;
}


void chain::undoLastRotamerChange()
{	int index = itsBuffers[0].getIndexInChainBuffer();
	if (index < 0)
	{	cout<< "Cannot Undo Rotamer Change: indexInChainBuffer lost" << endl;
		resetUndoBuffer();
		resetRedoBuffer();
		return;
	}
	vector<UInt> rotamerIndex = itsBuffers[0].getRotamerIndexBuffer();
	vector<vector<double> > sidechainDihedralAngles = itsBuffers[0].getSidechainDihedralAngleBuffer();
	setRotamerWithoutBuffering(index,rotamerIndex,sidechainDihedralAngles);

// here we need to recoup the exact last state of the residue once
// small rotations are allowed. this will be done through the
// sidechainDihedralAngles vector

	if (sidechainDihedralAngles.size() == 0)
	{	cout << "Cannot undo: sideChainDihedralAngleBuffer lost" << endl;
		resetUndoBuffer();
		resetRedoBuffer();
		return;
	}
	for (UInt i=0; i<sidechainDihedralAngles.size(); i++)
	{	for (UInt j=0; j<sidechainDihedralAngles[i].size(); j++)
		{	itsResidues[index]->setChi(i,j,sidechainDihedralAngles[i][j]);
		}
	}
	resetUndoBuffer();
	resetRedoBuffer();
#ifdef _CHAIN_DEBUG
	listDihedrals();
#endif
}

vector<chainModBuffer> chain::performRandomRotamerRotation(ran& _ran)
{
	if (itsRepackActivePositionMap.size() == 0)
	{
		resetUndoBuffer();
		resetRedoBuffer();
		cout << endl << "ABORTED" << endl;
		return itsBuffers;
	}
	int residueToModify = chooseTargetResidue(_ran);
	cout << residueToModify << " ";
	if (residueToModify == -1)
	{
		resetUndoBuffer();
		resetRedoBuffer();
		cout << endl << "ABORTED" << endl;
		return itsBuffers;
	}
	cout << itsResidues[residueToModify]->getType() << " ";
	int branchpoint = chooseTargetBranchpoint(residueToModify,_ran);
	cout << branchpoint << " ";
	if (branchpoint < 0)
	{
		resetUndoBuffer();
		resetRedoBuffer();
		cout << endl << "ABORTED" << endl;
		return itsBuffers;
	}
	int dihedralIndex = chooseTargetDihedralIndex(residueToModify,
		branchpoint,_ran);
	cout << dihedralIndex << " ";
	if (dihedralIndex < 0)
	{
		resetUndoBuffer();
		resetRedoBuffer();
		cout << endl << "ABORTED" << endl;
		return itsBuffers;
	}
	double newAngle = chooseTargetDihedralAngle(residueToModify,
		branchpoint, dihedralIndex,_ran);
	cout << newAngle << " ";
	setDihedralWithBuffering(residueToModify, branchpoint,
		dihedralIndex, newAngle);
	return itsBuffers;
}


vector<chainModBuffer> chain::performRandomRotamerRotation(ran& _ran, vector <int> _position)
{
    if (itsRepackActivePositionMap.size() == 0)
    {
        resetUndoBuffer();
        resetRedoBuffer();
        cout << endl << "ABORTED" << endl;
        return itsBuffers;
    }
    int residueToModify = _position[2];
	itsLastTargetResidue = residueToModify;
    cout << residueToModify << " ";
    if (residueToModify == -1)
    {
        resetUndoBuffer();
        resetRedoBuffer();
        cout << endl << "ABORTED" << endl;
        return itsBuffers;
    }
    cout << itsResidues[residueToModify]->getType() << " ";
    int branchpoint = chooseTargetBranchpoint(residueToModify,_ran);
    cout << branchpoint << " ";
    if (branchpoint < 0)
    {
        resetUndoBuffer();
        resetRedoBuffer();
        cout << endl << "ABORTED" << endl;
        return itsBuffers;
    }
    int dihedralIndex = chooseTargetDihedralIndex(residueToModify,
        branchpoint,_ran);
    cout << dihedralIndex << " ";
    if (dihedralIndex < 0)
    {
        resetUndoBuffer();
        resetRedoBuffer();
        cout << endl << "ABORTED" << endl;
        return itsBuffers;
    }
    double newAngle = chooseTargetDihedralAngle(residueToModify,
        branchpoint, dihedralIndex,_ran);
    cout << newAngle << " ";
    setDihedralWithBuffering(residueToModify, branchpoint,
        dihedralIndex, newAngle);
    return itsBuffers;
}

void chain::commitLastRotamerRotation()
{   resetUndoBuffer();
    resetRedoBuffer();
#ifdef _CHAIN_DEBUG
    listDihedrals();
#endif
}



void chain::undoLastRotamerRotation()
{	int index = itsBuffers[0].getIndexInChainBuffer();
	if (index < 0)
	{	cout << "Cannot undo rotation: indexInChainBuffer lost" << endl;
		resetUndoBuffer();
		resetRedoBuffer();
		return;
	}
	vector<vector<double> > sidechainDihedralAngles = itsBuffers[0].getSidechainDihedralAngleBuffer();
	if (sidechainDihedralAngles.size() == 0)
	{	cout << "Cannot undo: sideChainDihedralAngleBuffer lost" << endl;
		resetUndoBuffer();
		resetRedoBuffer();
		return;
	}
	for (UInt i=0; i<sidechainDihedralAngles.size(); i++)
	{	for (UInt j=0; j<sidechainDihedralAngles[i].size(); j++)
		{	itsResidues[index]->setChi(i,j,sidechainDihedralAngles[i][j]);
		}
	}
	resetUndoBuffer();
	resetRedoBuffer();
#ifdef _CHAIN_DEBUG
	listDihedrals();
#endif
}

vector <dblVec> chain::saveCoords(UInt resIndex)
{
	UInt nAtoms = getNumAtoms(resIndex);
	vector <dblVec> allCoords;
	for (UInt i=0; i<nAtoms; i++)
	{
		dblVec coords = getCoords(resIndex, i);
		allCoords.push_back(coords);
	}
	return allCoords;
}

void chain::setAllCoords(UInt resIndex, vector<dblVec> allCoords)
{
	UInt nAtoms = getNumAtoms(resIndex);
	for (UInt i=0; i<nAtoms; i++)
	{
		setCoords(resIndex, i, allCoords[i]);
	}
}

vector <chainModBuffer> chain::saveCurrentState()
{
	stateBuffers.resize(0);
	for (UInt i = 0; i < itsResidues.size(); i ++)
	{
		if (!isCofactor(i)){
			chainModBuffer temp;
			temp.setIndexInChainBuffer(i);
			temp.setResidueIdentityBuffer(itsResidues[i]->getTypeIndex());
			vector<UInt> tempRotamer;
			tempRotamer = itsChainPositions[i]->getCurrentRotamerIndex();
			temp.setRotamerIndexBuffer(tempRotamer);
			vector< vector< double> > tempAngles;
			tempAngles = itsResidues[i]->getSidechainDihedralAngles();
			temp.setSidechainDihedralAngleBuffer(tempAngles);
			vector< double> tempBBAngles(3);
			tempBBAngles.push_back(itsResidues[i]->getPhi());
			tempBBAngles.push_back(itsResidues[i]->getPsi());
			tempBBAngles.push_back(itsResidues[i]->getBetaChi());
			temp.setBackboneDihedralAngleBuffer(tempBBAngles);
			stateBuffers.push_back(temp);
		}
	}	
	return stateBuffers;
}

vector <chainModBuffer> chain::saveCurrentState(vector <int> _indices)
{
	stateBuffers.resize(0);
	for (UInt x = 0; x < _indices.size(); x ++)
	{
		int i = _indices[x];
		if (i >= 0 && (UInt)i < itsResidues.size())
		{
			if (!isCofactor(i)){
				chainModBuffer temp;
				temp.setIndexInChainBuffer(i);
				temp.setResidueIdentityBuffer(itsResidues[i]->getTypeIndex());
				vector<UInt> tempRotamer;
				tempRotamer = itsChainPositions[i]->getCurrentRotamerIndex();
				temp.setRotamerIndexBuffer(tempRotamer);
				vector< vector< double> > tempAngles;
				tempAngles = itsResidues[i]->getSidechainDihedralAngles();
				temp.setSidechainDihedralAngleBuffer(tempAngles);
				vector< double> tempBBAngles(3);
				tempBBAngles.push_back(itsResidues[i]->getPhi());
				tempBBAngles.push_back(itsResidues[i]->getPsi());
				tempBBAngles.push_back(itsResidues[i]->getBetaChi());
				temp.setBackboneDihedralAngleBuffer(tempBBAngles);
				stateBuffers.push_back(temp);
			}
		}
	}
	return stateBuffers;
}

void chain::commitState()
{
	stateBuffers.resize(0);
	return;
}

void chain::undoState()
{
	for (UInt x = 0; x < stateBuffers.size(); x ++)
	{
		int index = stateBuffers[x].getIndexInChainBuffer();
		if (index < 0)
		{
			cout << "Cannot undo mutation: indexInChainBuffer lost" << endl;
			stateBuffers[x].printAll();
			stateBuffers.resize(0);
			return;
		}
		int resID = stateBuffers[x].getResidueIdentityBuffer();
		if (resID <0)
		{
			cout << "Cannot undo: residueIdentityBuffer lost" << endl;
			stateBuffers[x].printAll();
			stateBuffers.resize(0);
			return;
		}
		residue* pOldRes = itsResidues[index];
		itsResidues[index] = pOldRes->mutate(resID);
		delete pOldRes;

		itsChainPositions[index]->setCurrentResIndex(itsResidues[index]->getTypeIndex());
		
		vector<double> bbAngles = stateBuffers[x].getBackboneDihedralAngleBuffer();
		itsResidues[index]->setPhi(bbAngles[0]);
		itsResidues[index]->setPsi(bbAngles[1]);
		itsResidues[index]->setBetaChi(bbAngles[2]);
		
		vector< UInt> rotIndex = stateBuffers[x].getRotamerIndexBuffer();
		setRotamerWithoutBuffering(index,rotIndex);

		vector<vector<double> > angles = stateBuffers[x].getSidechainDihedralAngleBuffer();

		UInt bpt = angles.size();
		if (bpt != 0)
		{
			for (UInt i=0; i < bpt; i++)
			{
				UInt numAngles = angles[i].size();
				if (numAngles !=0)
				{
					for (UInt j=0;j<numAngles;j++)
					{
						double angle = angles[i][j];
						itsResidues[index]->setChi(i,j,angle);
					}
				}
			}
		}
		stateBuffers[x].resetAllBuffers();
	}
#ifdef _CHAIN_DEBUG
	listDihedrals();
#endif
	stateBuffers.resize(0);
	return;
}

void chain::copyState(vector <chainModBuffer> _externalRedoBuffer)
{
	//cout << "chain::repeatModification(chainModBuffer _redoBuffer) called" << endl;
	//cout << "this function not yet implemented" << end
	stateBuffers.resize(0);

	for (UInt x = 0; x < _externalRedoBuffer.size(); x ++)
	{
		int index = _externalRedoBuffer[x].getIndexInChainBuffer();
		//cout << "externalRedoBufferIndex = " << index << endl;

		if (index != -1)
		{
			// buffer the residue that's there into the undo buffer
			bufferResidueIntoUndoBuffer(index);
			UInt bufferResType = _externalRedoBuffer[x].getResidueIdentityBuffer();
			// compare this to what is there already
			if (bufferResType != itsResidues[index]->getTypeIndex())
			{	// we need to perform a mutation
				//cout << "repeating MUT" << endl;
				mutateWithoutBuffering(index,bufferResType);
			}
			vector<double> bufferBBDihedrals = _externalRedoBuffer[x].getBackboneDihedralAngleBuffer();
			if (bufferBBDihedrals != itsResidues[index]->getBackboneAngles())
			{
				itsResidues[index]->setPhi(bufferBBDihedrals[0]);
				itsResidues[index]->setPsi(bufferBBDihedrals[1]);
			}
			vector<UInt> bufferRotamerIndex = _externalRedoBuffer[x].getRotamerIndexBuffer();
			if (bufferRotamerIndex != itsChainPositions[index]->getCurrentRotamerIndex())
			{
				// we need to set the rotameric state
				//cout << "repeating ROTC" << endl;
				setRotamerWithoutBuffering(index,bufferRotamerIndex);
			}

			vector<vector<double> > bufferDihedrals = _externalRedoBuffer[x].getSidechainDihedralAngleBuffer();
			if (bufferDihedrals != itsResidues[index]->getSidechainDihedralAngles())
			{
				//cout << "repeating DIHEDRAL CHANGE" << endl;
				for (UInt i=0; i<bufferDihedrals.size(); i++)
				{	for (UInt j=0; j<bufferDihedrals.size(); j++)
					{
						setDihedralWithoutBuffering(index,i,j,bufferDihedrals[i][j]);
					}
				}
			}
		}
	}
}


void chain::listAllowedRotamers(UInt _indexInChain) const
{	if (_indexInChain < itsChainPositions.size())
	{	if (itsChainPositions[_indexInChain])
		{	itsChainPositions[_indexInChain]->listAllowedRotamers();
		}
	}
}

void chain::setResNotAllowed(const UInt _indexInChain, const UInt _aaType)
{	if( _indexInChain >= 0 && _indexInChain < itsChainPositions.size() )
	{	if (itsChainPositions[_indexInChain])
		{	if (_aaType < residue::getDataBaseSize())
			{	itsChainPositions[_indexInChain]->setResNotAllowed(_aaType);
        	}
		}
		else
		{	//cout << "Residue not in active set..." << endl;
		}
	}
}

void chain::setRotamerNotAllowed(const UInt _indexInChain, const UInt _aaType, const UInt _bpt, const UInt _rotamer)
{
	if ( _indexInChain >= 0 && _indexInChain < itsChainPositions.size() )
	{
		itsChainPositions[_indexInChain]->setRotamerNotAllowed( _aaType, _bpt, _rotamer );
	}
	else
	{
		cout << "Error in chain::setRotamerNotAllowed." << endl;
		cout << " position " << _indexInChain << " is out of range on chain." << endl;
	}
	return;
}

UIntVec chain::getAllowedRotamers(const UInt _indexInChain, const UInt _aaType, const UInt _bpt)
{
	UIntVec allowedRotamers;
	allowedRotamers.resize(0);
	if ( _indexInChain >=0 && _indexInChain < itsChainPositions.size() )
	{
		allowedRotamers = itsChainPositions[_indexInChain]->getAllowedRotamers(_aaType, _bpt);
	}
	else
	{
        cout << "Error in chain::getAllowedRotamers." << endl;
        cout << " position " << _indexInChain << " is out of range on chain." << endl;
	}

	return allowedRotamers;
}

vector <UIntVec> chain::getAllowedRotamers(const UInt _indexInChain, const UInt _aaType)
{
    vector <UIntVec> allowedRotamers;
    allowedRotamers.resize(0);
    if ( _indexInChain >=0 && _indexInChain < itsChainPositions.size() )
    {
        allowedRotamers = itsChainPositions[_indexInChain]->getAllowedRotamers(_aaType);
    }
    else
    {
        cout << "Error in chain::getAllowedRotamers." << endl;
        cout << " position " << _indexInChain << " is out of range on chain." << endl;
    }

    return allowedRotamers;
}

UIntVec chain::getResAllowed(const UInt _indexInChain)
{
	UIntVec allowedResidues;
	if (_indexInChain >= 0 && _indexInChain < itsChainPositions.size())
	{
		if (itsChainPositions[_indexInChain]) // is residue active?
		{
			allowedResidues = itsChainPositions[_indexInChain]->getResAllowed();
			return allowedResidues;
		}
		else
		{
			cout << "Residue " << _indexInChain << " is not in the active set " << endl;
			allowedResidues.resize(0);
			return allowedResidues; // return null
		}
	}
	else
	{
		cout << "Residue " << _indexInChain << " is not in the chain " << endl;
		allowedResidues.resize(0);
		return allowedResidues; // return null
	}
}

void chain::setResAllowed(const UInt _indexInChain, const UInt _aaType)
{	if( _indexInChain >= 0 && _indexInChain < itsChainPositions.size() )
	{	if (itsChainPositions[_indexInChain])
		{	if (_aaType < residue::getDataBaseSize())
			{	itsChainPositions[_indexInChain]->setResAllowed(_aaType);
        	}
		}
		else
		{	//cout << "Residue not in active set..." << endl;
		}
	}
}

void chain::setOnlyNativeIdentity(const UInt _indexInChain)
{	if ( _indexInChain >=0 && _indexInChain < itsChainPositions.size() )
	{	if (itsChainPositions[_indexInChain])
		{	UInt nativeType = itsResidues[_indexInChain]->getTypeIndex();
			cout << "Native type: " << nativeType << endl;
			itsChainPositions[_indexInChain]->setOnlyAllowedIdentity(nativeType);
		}
	}
}

// selects a position to modify without actually changing any local values of the object
int chain::chooseNextTargetPosition(ran& _ran)
{
	UInt activePositionSize = itsRepackActivePositionMap.size();
	if (activePositionSize != 0)
    {
		UInt choice = UInt(_ran.getNext()*activePositionSize);
        return itsRepackActivePositionMap[choice];
    }
    itsLastTargetResidue = -1;
    return -1;
}

int chain::chooseTargetResidue(ran& _ran)
{	UInt activePositionSize = itsRepackActivePositionMap.size();
	// cout << "apmsize = " << activePositionSize << endl;
	if (activePositionSize != 0)
	{	UInt choice = UInt(_ran.getNext()*activePositionSize);
		itsLastTargetResidue = itsRepackActivePositionMap[choice];
		return itsRepackActivePositionMap[choice];
	}
	itsLastTargetResidue = -1;
	return -1;
}

int chain::chooseTargetIdentity(const UInt _indexInChain, ran& _ran)
{	if (_indexInChain >= 0 && _indexInChain < itsChainPositions.size())
	{	if (itsChainPositions[_indexInChain])
		{	return (itsChainPositions[_indexInChain])->chooseResidueFromLibrary(_ran);
		}
		else
		{	cout << "Error from chain::chooseTargetIdentity" << endl;
			cout << "Residue not in active set : " << _indexInChain << endl;
		}
	}
	return -1;
}

int chain::chooseTargetBranchpoint(UInt _indexInChain, ran& _ran)
{
	if (itsChainPositions[_indexInChain])
	{
		int bpt =  (itsChainPositions[_indexInChain])->chooseBranchpoint(_ran);
		return bpt;
	}
	else
	{	cout << "Error from chain::chooseTargetBranchpoint" << endl;
		cout << "Residue not in active set : " << itsBuffers[0].getIndexInChainBuffer() << endl;
	}
	return -1;
}

int chain::chooseTargetRotamer(const UInt _indexInChain, const UInt _bpt, ran& _ran)
{
	/*
	cout << "sizeofChainPositions: " << itsChainPositions.size() << endl;
	cout << "_indexInChain: " << _indexInChain << endl;
	for (UInt i=0; i<itsChainPositions.size(); i++)
	{	cout << i << " " << itsChainPositions[i] << endl;
	}
	*/

	if (_indexInChain >=0 && _indexInChain < itsChainPositions.size())
	{	if (itsChainPositions[_indexInChain])
		{	return itsChainPositions[_indexInChain]->chooseRotamerFromLibrary(_bpt,_ran);
		}
		else
		{	cout << "Error from chain::chooseTargetRotamer" << endl;
			cout << "Residue not in active set : " << _indexInChain << endl;
		}
	}
	return -1;
}

void chain::setRelativeChi(const UInt _resIndex, const UInt _bpt, const UInt _chi, const double _angle)
{
	itsResidues[_resIndex]->setChiByDelta(_bpt, _chi, _angle);
	return;
}

void chain::setDihedralWithBuffering(const UInt _index, const UInt _bpt, const UInt _dihedralIndex, const double _angle)
{
	bufferResidueIntoUndoBuffer(_index);
	itsResidues[_index]->setChi(_bpt, _dihedralIndex, _angle);
	bufferResidueIntoRedoBuffer(_index);
	return;
}

void chain::setDihedralWithoutBuffering(const UInt _index, const UInt _bpt,
	const UInt _dihedralIndex, const double _angle)
{
	itsResidues[_index]->setChi(_bpt, _dihedralIndex, _angle);
	return;
}

vector< vector< double> > chain::getSidechainDihedralAngles(UInt _indexInChain)
{
	vector< vector< double> > theAngles;
	theAngles = itsResidues[_indexInChain]->getSidechainDihedralAngles();
	return theAngles;
}

void chain::setSidechainDihedralAngles(UInt _indexInChain, vector <vector <double> > Angles)
{
	for (UInt i=0; i<Angles.size(); i++)
	{	for (UInt j=0; j<Angles[i].size(); j++)
		{	itsResidues[_indexInChain]->setChi(i,j,Angles[i][j]);
		}
	}
	itsResidues[_indexInChain]->setMoved();
	return;
}

int chain::chooseTargetDihedralIndex(const UInt _index, const UInt _bpt, ran& _ran)
{
	UInt resType = itsResidues[_index]->getTypeIndex();
	int numDihedrals = residue::dataBase[resType].getNumberOfChis(_bpt);
	if( numDihedrals > 0 )
	{	return _ran.getNext(0,numDihedrals-1);
	}
	return -1;
}

double chain::chooseTargetDihedralAngle(const UInt _index, const UInt _bpt,
		  const UInt _angleIndex, ran& _ran)
{
	double theIdealAngle = 0.0;
	vector<vector< double> > theActualAngles;
	theActualAngles = itsResidues[_index]->getSidechainDihedralAngles();
	double angleToChange = theActualAngles[_bpt][_angleIndex];
	if (angleToChange >= 0.0 && angleToChange < 120.0 )
	{	theIdealAngle = 60.0;
	}
	else if (angleToChange >= -120.0 && angleToChange < 0.0 )
	{	theIdealAngle = -60.0;
	}
	else if ( (angleToChange >= 120.0 && angleToChange < 180.0) ||
			  (angleToChange >= -180.0 && angleToChange < -120.0) )
	{	theIdealAngle = 180;
	}
	// Uniform deviate - i.e. any angle has an equal probability of
	// being chosen
	double theTargetAngle = _ran.getNext(theIdealAngle - 30,theIdealAngle + 30.0);
	//cout << "\t angle chosen at chain::chooseTargetDihedralAngle is: " << theTargetAngle << "\t";
	return theTargetAngle;
}

void chain::translate(const dblVec& _dblVec)
{	for (UInt i=0; i<itsResidues.size(); i++)
	{	itsResidues[i]->translate(_dblVec);
	}
}

void chain::translate(const double _x,const double _y,const double _z)
{	dblVec vec(3);
	vec[0] = _x;
	vec[1] = _y;
	vec[2] = _z;
	translate(vec);
}

void chain::translateR(const double _x,const double _y,const double _z)
{	dblVec vec(3);
	vec[0] = _z;
	vec[1] = _y;
	vec[2] = _x;
	translate(vec);
}

void chain::transform(const dblMat& _dblMat)
{	for (UInt i=0; i<itsResidues.size(); i++)
	{	itsResidues[i]->transform(_dblMat);
	}
}

void chain::rotate(const axis _axis,const double _theta)
{	// The default is to set this point to the origin
	point origin;
	origin.setCoords(0.0,0.0,0.0);

	dblVec vec(3);
#ifdef USE_SVMT
	for (UInt i=0; i<vec.extent(); i++)
	{	vec[i] = 0.0;
	}
#else
	for (int i=0; i<vec.dim(); i++)
	{	vec[i] = 0.0;
	}
#endif
	if (_axis == X_axis)
	{	vec[0] = 1.0;
	}
	else if (_axis == Y_axis)
	{	vec[1] = 1.0;
	}
	else if (_axis == Z_axis)
	{	vec[2]	= 1.0;
	}
	rotate(origin, vec, _theta);
}

void chain::rotateRelative(const axis _axis,const double _theta)
{	
	dblVec COMcoords(3);
	COMcoords[0] = 0.0, COMcoords[1] = 0.0, COMcoords[2] = 0.0;
	UInt totalAtoms = 0;
	for (UInt i=0; i<itsResidues.size(); i++)
	{	
		UInt nAtoms = itsResidues[i]->getNumAtoms();
		dblVec coords(3);
		for (UInt j=0; j<nAtoms; j++)
		{
			dblVec coords = itsResidues[i]->getCoords(j);
			COMcoords[0] = COMcoords[0]+coords[0], COMcoords[1] = COMcoords[1]+coords[1], COMcoords[2] = COMcoords[2]+coords[2];
			totalAtoms++;
		}
	}
	COMcoords[0] = COMcoords[0]/totalAtoms, COMcoords[1] = COMcoords[1]/totalAtoms, COMcoords[2] = COMcoords[2]/totalAtoms;
	translate(COMcoords[0]*-1,COMcoords[1]*-1,COMcoords[2]*-1);
	point origin;
	origin.setCoords(0,0,0);

	dblVec vec(3);
#ifdef USE_SVMT
	for (UInt i=0; i<vec.extent(); i++)
	{	vec[i] = 0.0;
	}
#else
	for (int i=0; i<vec.dim(); i++)
	{	vec[i] = 0.0;
	}
#endif
	if (_axis == X_axis)
	{	vec[0] = 1.0;
	}
	else if (_axis == Y_axis)
	{	vec[1] = 1.0;
	}
	else if (_axis == Z_axis)
	{	vec[2]	= 1.0;
	}
	rotate(origin, vec, _theta);
	translate(COMcoords[0],COMcoords[1],COMcoords[2]);
}
	
void chain::rotate(const point& _point,const dblVec& _R_axis, const double _theta)
{	for (UInt i=0; i<itsResidues.size(); i++)
	{	itsResidues[i]->rotate(_point, _R_axis, _theta);
	}
}

void chain::updateResiduesPerTurnType()
{
	if (itsResidues.size() > 1){
		for(UInt i=0; i<itsResidues.size(); i++)
		{	
			UInt RPT = getBackboneSequenceType(i);
			itsResidues[i]->setResiduesPerTurnType(RPT);
		}
	}
}

double chain::getResiduesPerTurn(const UInt _resIndex)
{
	double phi, psi, rpt = 0.0;
	if (!isCofactor(_resIndex))
	{
		if (_resIndex == 0){
			psi = itsResidues[_resIndex]->getPsi();
			phi = itsResidues[_resIndex+1]->getPhi();
		}
		else if (_resIndex == itsResidues.size()-1){
			psi = itsResidues[_resIndex-1]->getPsi();
			phi = itsResidues[_resIndex]->getPhi();
		}
		else{
			phi = itsResidues[_resIndex]->getPhi();
			psi = itsResidues[_resIndex]->getPsi();
		}
		return getResiduesPerTurn(phi, psi);
	}
	else{return 0;}
}

double chain::getResiduesPerTurn(double phi, double psi)
{
	double rpt;
	double radSum = ((phi+psi)/2)*PI/180;
	double radDiff = ((psi-phi)/2)*PI/180;
	double radTheta = 2*acos(-0.8235*sin(radSum)-0.0222*sin(radDiff));
	double theta = radTheta*180/PI;
	if (theta > 180) {rpt = 360/(theta-360);} else {rpt = 360/theta;}
	double rise = 1/sin(radTheta/2)*(2.999*cos(radSum)-0.657*cos(radDiff));
	if (rise < 0) {rpt = rpt*-1;}
	return rpt;
}

UInt chain::getBackboneSequenceType(const UInt _resIndex)
{
	double phi;
	if (_resIndex == 0){
		phi = itsResidues[_resIndex+1]->getPhi();
	}
	else{
		phi = itsResidues[_resIndex]->getPhi();
	}
	double RPT = getResiduesPerTurn(_resIndex);
	return getBackboneSequenceType(RPT, phi);
}

UInt chain::getBackboneSequenceType(double RPT, double phi)
{
	if (RPT <= -4.83 && phi <= 0)					{return 0;}
	if (RPT > -4.83  && RPT <= -4.13 && phi <= 0)	{return 1;}
	if (RPT > -4.13  && RPT <= -3.43 && phi <= 0)	{return 2;}
	if (RPT > -3.43  && RPT <= -2.73 && phi <= 0)	{return 3;}
	if (RPT > -2.73  && RPT < -2.00 && phi <= 0)	{return 4;}
	if (RPT >=  2.00  && RPT <=  2.73 && phi <= 0)	{return 5;}
	if (RPT >  2.73  && RPT <=  3.43 && phi <= 0)	{return 6;}
	if (RPT >  3.43  && RPT <=  4.13 && phi <= 0)	{return 7;}
	if (RPT >  4.13  && RPT <=  4.83 && phi <= 0)	{return 8;}
	if (RPT >  4.83 && phi <= 0)					{return 9;}
	if (RPT <= -4.83 && phi > 0)					{return 10;}
	if (RPT > -4.83  && RPT <= -4.13 && phi > 0)	{return 11;}
	if (RPT > -4.13  && RPT <= -3.43 && phi > 0)	{return 12;}
	if (RPT > -3.43  && RPT <= -2.73 && phi > 0)	{return 13;}
	if (RPT > -2.73  && RPT < -2.00 && phi > 0)		{return 14;}
	if (RPT >=  2.00  && RPT <=  2.73 && phi > 0)	{return 15;}
	if (RPT >  2.73  && RPT <=  3.43 && phi > 0)	{return 16;}
	if (RPT >  3.43  && RPT <=  4.13 && phi > 0)	{return 17;}
	if (RPT >  4.13  && RPT <=  4.83 && phi > 0)	{return 18;}
	if (RPT >  4.83 && phi > 0)						{return 19;}
	return 0;
}

double chain::getPhi(const UInt _indexInChain)
{
	double tempdouble;
	if (_indexInChain < itsResidues.size())
	{
		tempdouble = (itsResidues[_indexInChain])->getPhi();
	}
	else
	{
		cout << "Error from chain::getPhi" << endl;
		cout << "Residue indexInChain out of range: " << _indexInChain << endl;
		tempdouble = 1000.0;
	}
	return tempdouble;
}

double chain::getPsi(const UInt _indexInChain)
{
	double tempdouble;
	if (_indexInChain < itsResidues.size())
	{
		tempdouble = (itsResidues[_indexInChain])->getPsi();
	}
	else
	{
		cout << "Error from chain::getPsi" << endl;
		cout << "Residue indexInChain out of range: " << _indexInChain << endl;
		tempdouble = 1000.0;
	}
	return tempdouble;
}

double chain::getAngle(const UInt _indexInChain, UInt angleType)
{
	double tempdouble;
	if (_indexInChain < itsResidues.size())
	{
		tempdouble = (itsResidues[_indexInChain])->getAngle(angleType);
	}
	else
	{
		cout << "Error from chain::getAngle" << endl;
		cout << "Residue indexInChain out of range: " << _indexInChain << endl;
		tempdouble = 1000.0;
	}
	return tempdouble;
}

int chain::setPsi(const UInt _index, double _psi)
{
	if (_index < itsResidues.size())
	{
		return itsResidues[_index]->setPsi(_psi);
	}
	else 
	{
		cout << "Residue index out of range: " << _index << endl;
		return -1;
	}
}

int chain::setPhi(const UInt _index, double _phi)
{
	if (_index < itsResidues.size())
	{
		return itsResidues[_index]->setPhi(_phi);
	}
	else
	{
		cout << "Residue index out of range: " << _index << endl;
		return -1;
	}
}

int chain::setDihedral(const UInt _resIndex, double _dihedral, UInt _angleType, UInt _direction)
{
	if (_resIndex < itsResidues.size())
	{
		return itsResidues[_resIndex]->setDihedral(_dihedral, _angleType, _direction);
	}
	else
	{
		cout << "Residue index out of range: " << _resIndex << endl;
		return -1;
	}
}

double chain::getAmide(const UInt _indexInChain)
{
	double tempdouble;
	if (_indexInChain < itsResidues.size())
	{
		tempdouble = (itsResidues[_indexInChain])->getAmide();
	}
	else
	{
		cout << "Error from chain::getAmide" << endl;
		cout << "Residue indexInChain out of range: " << _indexInChain << endl;
		tempdouble = 1000.0;
	}
	return tempdouble;
}


double chain::getVolume(UInt _method)
{
	double itsVolume = 0.0;
	for (UInt i = 0; i < itsResidues.size(); i ++)
	{
		itsVolume += itsResidues[i]->getVolume(_method);
	}
	return itsVolume;
}

double chain::getInterEnergy(const UInt _residue1, const UInt _atom1, chain* _other, const UInt _residue2, const UInt _atom2)
{
	if (_residue1 >=0 && _residue1 < itsResidues.size() && _residue2 >=0 && _residue2 < _other->itsResidues.size())
	{
		return itsResidues[_residue1]->getIntraEnergy(_atom1, _other->itsResidues[_residue2], _atom2);
	}
	else
	{
		cout << "ERROR in chain::getIntraEnergy(....) *residue index out of range." << endl;
		exit(1);
	}
}

double chain::getInterEnergy(const UInt _residue1, chain* _other, const UInt _residue2)
{
	if (_residue1 >= 0 && _residue1 < itsResidues.size() && _residue2 >= 0 && _residue2 < _other->itsResidues.size())
	{
		return itsResidues[_residue1]->interEnergy(_other->itsResidues[_residue2]);
	}
	else
	{
		cout << "ERROR in chain:getIntraEnergy(...)  residue index out of range." << endl;
		exit(1);
	}
}

double chain::intraEnergy()
{	double intraEnergy = 0.0;
	//cout << "chain::intraEnergy";
	for(UInt i=0; i<itsResidues.size(); i++)
	{	// choose not to include residue internal energy
		//cout << endl << "res " << i << "  ";
		double tempE = itsResidues[i]->intraEnergy();
		intraEnergy += tempE;
		//cout << tempE << " ";
		//double rotE = rotamerEnergy(i);
		//intraEnergy += rotE;
		//cout << rotE << " ";
		double interE = 0.0;
		for(UInt j=i+1; j<itsResidues.size(); j++)
		{	interE += itsResidues[i]->interEnergy(itsResidues[j]);
		}
		intraEnergy += interE;
		//cout << interE << " ";
	}
	//cout << " done!" << endl << "intraEnergy being returned is : ";
	// cout << intraEnergy << endl;
	return intraEnergy;
}

void chain::updateEnergy()
{	
	bool resI, resJ;
	double resEnergy;
	for(UInt i=0; i<itsResidues.size(); i++)
	{	
		resI = itsResidues[i]->getMoved(0);
		for(UInt j=i+1; j<itsResidues.size(); j++)
		{
			resJ =  itsResidues[j]->getMoved(0);
			if (resI || resJ){
				resEnergy = itsResidues[i]->interSoluteEnergy(itsResidues[j]), resEnergy /= 2;
				if(resI){itsResidues[i]->sumEnergy(resEnergy);}
				if(resJ){itsResidues[j]->sumEnergy(resEnergy);}
			}
		}
		if(resI){
			resEnergy = itsResidues[i]->intraSoluteEnergy();
			itsResidues[i]->sumEnergy(resEnergy);
		}
	}
}

void chain::updateEnergy(chain* _other)
{
	bool resI, resJ;
	double resEnergy;
	for(UInt i=0; i<itsResidues.size(); i++)
	{
		resI = itsResidues[i]->getMoved(0);
		for(UInt j=0; j<_other->itsResidues.size(); j++)
		{
			resJ = _other->itsResidues[j]->getMoved(0);
			if (resI || resJ){
				resEnergy = itsResidues[i]->interSoluteEnergy(_other->itsResidues[j]), resEnergy /= 2;
				if(resI){itsResidues[i]->sumEnergy(resEnergy);}
				if(resJ){_other->itsResidues[j]->sumEnergy(resEnergy);}
			}
		}
	}
}

void chain::polarizability()
{
	bool resI, resJ;
	for(UInt i=0; i<itsResidues.size(); i++)
	{	
		resI = itsResidues[i]->getMoved(0);
		for(UInt j=i+1; j<itsResidues.size(); j++)
		{	
			resJ = itsResidues[j]->getMoved(0);
			if (resI || resJ){
				itsResidues[i]->polarizability(itsResidues[j]);
			}
		}
		if (resI){
			itsResidues[i]->polarizability();
		}
	}
}

void chain::polarizability(chain* _other)
{
	bool resI, resJ;
	for(UInt i=0; i<itsResidues.size(); i++)
	{
		resI = itsResidues[i]->getMoved(0);
		for(UInt j=0; j<_other->itsResidues.size(); j++)
		{
			resJ = _other->itsResidues[j]->getMoved(0);
			if (resI || resJ){
				itsResidues[i]->polarizability(_other->itsResidues[j]);
			}
		}
	}
}

void chain::calculateDielectrics()
{
	for(UInt i=0; i<itsResidues.size(); i++)
	{	
		if (itsResidues[i]->getMoved(0)){
			itsResidues[i]->calculateDielectrics();
		}
	}
}

void chain::updateMovedDependence(UInt _EorC)
{
	for(UInt i=0; i<itsResidues.size(); i++)
	{	
		if (itsResidues[i]->getCheckMovedDependence(_EorC)){
			for(UInt j=i+1; j<itsResidues.size(); j++)
			{	
				itsResidues[i]->updateMovedDependence(itsResidues[j], _EorC);
			}
		}
	}
}

void chain::updateMovedDependence(chain* _other, UInt _EorC)
{
	for(UInt i=0; i<itsResidues.size(); i++)
	{
		if (itsResidues[i]->getCheckMovedDependence(_EorC)){
			for(UInt j=0; j<_other->itsResidues.size(); j++)
			{
				itsResidues[i]->updateMovedDependence(_other->itsResidues[j], _EorC);
			}
		}
	}
}


void chain::setMoved(bool _moved, UInt _EorC)
{
	for(UInt i=0; i<itsResidues.size(); i++)
	{	
		itsResidues[i]->setMoved(_moved, _EorC);
	}
}

double chain::getEnergy()
{
	double Energy = 0.0;
	for(UInt i=0; i<itsResidues.size(); i++)
	{	
		Energy += itsResidues[i]->getEnergy();
	}
	return Energy;
}


void chain::updateClashes()
{
	bool resI, resJ;
	UInt clashes;
	for(UInt i=0; i<itsResidues.size(); i++)
	{	
		resI = itsResidues[i]->getMoved(1);
		for(UInt j=i+1; j<itsResidues.size(); j++)
		{
			resJ =  itsResidues[j]->getMoved(1);
			if (resI || resJ){
				clashes = itsResidues[i]->getNumHardClashes(itsResidues[j]), clashes /= 2;
				if(resI){itsResidues[i]->sumClashes(clashes);}
				if(resJ){itsResidues[j]->sumClashes(clashes);}
			}
		}
		if(resI){
			clashes = itsResidues[i]->getNumHardClashes();
			itsResidues[i]->sumClashes(clashes);
		}
	}
}

void chain::updateClashes(chain* _other)
{
	bool resI, resJ;
	UInt clashes;
	for(UInt i=0; i<itsResidues.size(); i++)
	{
		resI = itsResidues[i]->getMoved(1);
		for(UInt j=0; j<_other->itsResidues.size(); j++)
		{
			resJ = _other->itsResidues[j]->getMoved(1);
			if (resI || resJ){
				clashes = itsResidues[i]->getNumHardClashes(_other->itsResidues[j]), clashes /= 2;
				if(resI){itsResidues[i]->sumClashes(clashes);}
				if(resJ){_other->itsResidues[j]->sumClashes(clashes);}
			}
		}
	}
}

UInt chain::getClashes()
{
	UInt clashes = 0;
	for(UInt i=0; i<itsResidues.size(); i++)
	{	
		clashes += itsResidues[i]->getClashes();
	}
	return clashes;
}


void chain::updateBackboneClashes()
{
	bool resI, resJ;
	UInt clashes;
	for(UInt i=0; i<itsResidues.size(); i++)
	{	
		resI = itsResidues[i]->getMoved(2);
		for(UInt j=i+1; j<itsResidues.size(); j++)
		{
			resJ =  itsResidues[j]->getMoved(2);
			if (resI || resJ){
				clashes = itsResidues[i]->getNumHardBackboneClashes(itsResidues[j]), clashes /= 2;
				if(resI){itsResidues[i]->sumBackboneClashes(clashes);}
				if(resJ){itsResidues[j]->sumBackboneClashes(clashes);}
			}
		}
	}
}

void chain::updateBackboneClashes(chain* _other)
{
	bool resI, resJ;
	UInt clashes;
	for(UInt i=0; i<itsResidues.size(); i++)
	{
		resI = itsResidues[i]->getMoved(2);
		for(UInt j=0; j<_other->itsResidues.size(); j++)
		{
			resJ = _other->itsResidues[j]->getMoved(2);
			if (resI || resJ){
				clashes = itsResidues[i]->getNumHardBackboneClashes(_other->itsResidues[j]), clashes /= 2;
				if(resI){itsResidues[i]->sumBackboneClashes(clashes);}
				if(resJ){_other->itsResidues[j]->sumBackboneClashes(clashes);}
			}
		}
	}
}

UInt chain::getBackboneClashes()
{
	UInt clashes = 0;
	for(UInt i=0; i<itsResidues.size(); i++)
	{	
		clashes += itsResidues[i]->getBackboneClashes();
	}
	return clashes;
}

double chain::getPositionIntraEnergy(vector <int> _position)
{
	double intraEnergy = 0.0;
	intraEnergy += itsResidues[_position[2]]->intraEnergy();
	//intraEnergy += rotamerEnergy(_position[2]);
	//cout << " ... at residue level:  calculating intraenergy at position " << _position[0] << " " << _position[1] << " " << _position[2];
	for (UInt i = 0; i < itsResidues.size(); i++)
	{
		if ( (int)i != _position[2]) intraEnergy += itsResidues[_position[2]]->interEnergy(itsResidues[i]);
	}
	//cout << " " << intraEnergy << " i " << itsResidues[_position[2]]->intraEnergy() << " r " << rotamerEnergy(_position[2]) << endl;
	return intraEnergy;
}

double chain::interEnergy(chain* _other)
{	double interEnergy = 0.0;
	//cout << "chain::interEnergy";
	for(UInt i=0; i<itsResidues.size(); i++)
	{
	//	cout << "res " << i << "  ";
		for(UInt j=0; j<_other->itsResidues.size(); j++)
		{	interEnergy += itsResidues[i]->interEnergy(_other->itsResidues[j]);
		}
	}
	//cout << " done!" << endl;
	return interEnergy;
}

double chain::getPositionInterEnergy(vector <int> _position, chain* _other)
{
	double interEnergy = 0.0;
	for (UInt i = 0; i < _other->itsResidues.size(); i++)
	{
		interEnergy+= itsResidues[_position[2]]->interEnergy(_other->itsResidues[i]);
		//cout << interEnergy << endl;
	}

	return interEnergy;
}

double chain::rotamerEnergy()
{	double rotEnergy = 0;
	for(UInt i=0; i<itsResidues.size(); i++)
	{	rotEnergy += rotamerEnergy(i);
	}
	return rotEnergy;
}

double chain::rotamerEnergy(const UInt _indexInChain)
{	double theEnergy = 0.0;
	if (rotamer::getScaleFactor() != 0.0)
	{
		UInt itsType = itsResidues[_indexInChain]->getTypeIndex();
		if(	residue::dataBase[itsType].chiDefinitions.size() > 0 && 
			residue::dataBase[itsType].chiDefinitions[0].size() >= 4)
		{	if(itsChainPositions[_indexInChain] == 0)
			{	return 0.0;
			}
			vector<UInt> rotamerIndex =
					itsChainPositions[_indexInChain]->getCurrentRotamerIndex();
			/*
			for (UInt i=0; i<rotamerIndex.size(); i++)
			{
				cout << "rotIndex " << i << " = " << rotamerIndex[i] << endl;
			}
			*/
			UInt rotamerLibIndex =
					itsChainPositions[_indexInChain]->getRotamerLibIndex();
			for (UInt i=0; i< rotamerIndex.size(); i++)
			{	
				theEnergy = rotamer::getScaleFactor() * residue::dataBase[itsType].itsRotamerLibs[rotamerLibIndex]->getEnergy(i,rotamerIndex[i]);
			}
			/*
			cout << "itsType = " << itsType << " ";
			cout << "libIndex = " << rotamerLibIndex << " ";
			cout << "E = " << theEnergy << " ";
			residue::dataBase[itsType].itsRotamerLibs[rotamerLibIndex]->print();
			*/
			ASSERT(theEnergy < 1e50);
			ASSERT(theEnergy > -1e50);
		}
	}
	return theEnergy;
}
	

void chain::addResidue(residue* _pResidue)
{	if (itsResidues.size() == 0)
	{	_pResidue->setPrevRes(0);
		_pResidue->setNextRes(0);
	}
	else
	{	residue* oldLast = itsResidues[itsResidues.size()-1];
		_pResidue->setPrevRes(oldLast);
		_pResidue->setNextRes(0);
		(_pResidue->getPrevRes())->setNextRes(_pResidue);
	}

	itsResidues.push_back(_pResidue);
}

void chain::addChainPosition(chainPosition* _pChainPosition)
{	itsChainPositions.push_back(_pChainPosition);
}

void chain::addSecondaryStructure(secondaryStructure* _pSecondaryStructure)
{	itsSecondaryStructures.push_back(_pSecondaryStructure);
}

int chain::mapResNumToChainPosition(const int _resNum)
{
	for (UInt i=0; i<itsResidues.size(); i++)
	{	int tempresnum = itsResidues[i]->getResNum();
		if (_resNum == tempresnum)
		{	return i;
		}
	}
	cout << "Error - Residue # " << _resNum << " cannot be mapped" << endl;
	return 0;
}

//BUFFERING
void chain::resetUndoBuffer()
{
	//cout << "resetting undo buffer" << endl;
	itsBuffers[0].resetAllBuffers();
}

void chain::bufferResidueIntoUndoBuffer(UInt _index)
{
	//cout << "Buffering into undo buffer" << endl;
	itsBuffers[0].setIndexInChainBuffer(_index);
	itsBuffers[0].setResidueIdentityBuffer(itsResidues[_index]->getTypeIndex());
	vector<UInt> tempRotamer;
	tempRotamer = itsChainPositions[_index]->getCurrentRotamerIndex();
	itsBuffers[0].setRotamerIndexBuffer(tempRotamer);
	vector< vector< double> > tempAngles;
	tempAngles = itsResidues[_index]->getSidechainDihedralAngles();
	itsBuffers[0].setSidechainDihedralAngleBuffer(tempAngles);
	//cout << "contents of undo buffer" << endl;
	//itsBuffers[0].printAll();
}

void chain::resetRedoBuffer()
{
	//cout << "resetting redo buffer" << endl;
	itsBuffers[1].resetAllBuffers();
}

void chain::bufferResidueIntoRedoBuffer(UInt _index)
{
	//cout << "Buffering into redo buffer" << endl;
	itsBuffers[1].setIndexInChainBuffer(_index);
	itsBuffers[1].setResidueIdentityBuffer(itsResidues[_index]->getTypeIndex());
	vector<UInt> tempRotamer;
	tempRotamer = itsChainPositions[_index]->getCurrentRotamerIndex();
	itsBuffers[1].setRotamerIndexBuffer(tempRotamer);
	vector< vector< double> > tempAngles;
	tempAngles = itsResidues[_index]->getSidechainDihedralAngles();
	itsBuffers[1].setSidechainDihedralAngleBuffer(tempAngles);
	//cout << "contents of redo buffer" << endl;
	//itsBuffers[1].printAll();
}

void chain::setDihedrals(UInt _resIndex, UInt _bpt, vector <double> _dihedrals )
{
	for (UInt i = 0; i < _dihedrals.size(); i++)
	{
		setDihedralWithoutBuffering(_resIndex, _bpt, i, _dihedrals[i]);
	}
}

UInt chain::getNumChis(const UInt _resIndex, const UInt _bpt)
{
    if (_resIndex >=0 && _resIndex < itsResidues.size())
    {
        UInt resType = itsResidues[_resIndex]->getTypeIndex();
        return (UInt)residue::dataBase[resType].getNumberOfChis(_bpt);
    }
    else
    {
        cout << "ERROR in chain::getNumChis(...) - residue index out of range!" << endl;
    }
    return 0;
}

void chain::coilcoil(const double _pitch)
{
    for (UInt i = 0; i < itsResidues.size(); i++)
    {
        itsResidues[i]->coilcoil(_pitch);
    }
    return;
}

void chain::symmetryLinkResidueAtoB(UInt _masterResIndex, UInt _slaveResIndex)
{
    UInt master = 0;
    UInt slave = 0;
    bool masterFoundInIndependentList = false;
    bool slaveFoundInIndependentList = false;
    //remove slave residue from independent residue map
    for (UInt i = 0; i < itsIndependentPositions.size(); i++)
    {
        if (_slaveResIndex == itsIndependentPositions[i])
        {
            slave = i;
            slaveFoundInIndependentList = true;
        }
    }
    if (slaveFoundInIndependentList)
    {
        iterUINT firstUINT;
        firstUINT = itsIndependentPositions.begin();
        if ( slave < itsIndependentPositions.size())
        {
            itsIndependentPositions.erase(firstUINT + slave);
        }
        // also remove slave from resLinkageMap
        vector<int> subLinkages = itsResidueLinkageMap[slave];
        iterINTVEC firstINTVEC;
        firstINTVEC = itsResidueLinkageMap.begin();
        if (slave < itsResidueLinkageMap.size())
        {
            itsResidueLinkageMap.erase(firstINTVEC + slave);
        }
        if (subLinkages.size() == 1 && subLinkages[0] == -1)
        {
            // we're fine - no sublinkages
        }
        else
        {
            for (UInt i = 0; i < subLinkages.size(); i++)
            {
                symmetryLinkResidueAtoB(subLinkages[i], _masterResIndex);  // recursive sublinkages
            }
        }
    }
    // now find master residue and map
    for (UInt i = 0; i < itsIndependentPositions.size(); i++)
    {
        if (_masterResIndex == itsIndependentPositions[i])
        {
            master = i;
            masterFoundInIndependentList = true;
        }
        if (masterFoundInIndependentList)
        {
            //check to see if this is the first one
            if (itsResidueLinkageMap[master].size()==1 && itsResidueLinkageMap[master][0] == -1)
            {
                itsResidueLinkageMap[master][0] = _slaveResIndex;
            }
            else
            {
                itsResidueLinkageMap[master].push_back(_slaveResIndex);
            }
        }
        else
        {
            cout << "Error in chain::symmetryLinkResidueAtoB ..." << endl;
            cout << "\tIndependent residue not found in list." << endl;
        }
    }
}

double chain::getSelfEnergy(UInt _residueIndex)
{
    if (_residueIndex >=0 && _residueIndex < itsResidues.size())
    {
        double selfEnergy = 0.0;
        for (UInt i = 0; i < itsResidues.size(); i++)
        {
            //cout << "Residue " << _residueIndex << " with " << i << endl;
            selfEnergy += itsResidues[_residueIndex]->getSelfEnergy(itsResidues[i]);
        }
        return selfEnergy;
    }
    else
    {
        cout << "ERROR in chain::getSelfEnergy(...) residue index out of range" << endl;
        exit(1);
    }
    return -1;

}

// surface area operations...

void chain::initializeSpherePoints()
{
	for (UInt i = 0; i < itsResidues.size(); i ++)
	{
		itsResidues[i]->initializeSpherePoints();
	}

	return;
}

void chain::initializeSpherePoints(UInt _residue)
{
	if (_residue >= 0 && _residue < itsResidues.size() )
	{
		itsResidues[_residue]->initializeSpherePoints();
	}
	else
	{
		cout << "ERROR in chain::intializeSpherePoints ... residue index out of range." << endl;
	}

	return;
}

double chain::tabulateSurfaceArea()
{
	double surfaceArea = 0.0;
	for (UInt i = 0; i < itsResidues.size(); i++)
	{
		surfaceArea += itsResidues[i]->tabulateSurfaceArea();
	}

	return surfaceArea;
}

double chain::tabulateSurfaceArea(UInt _residue)
{
	double surfaceArea = 0.0;
	if (_residue >= 0 && _residue < itsResidues.size())
	{
		surfaceArea = itsResidues[_residue]->tabulateSurfaceArea();
	}
	else
	{
		cout << "ERROR in chain::tabulateSurfaceArea ... residue index out of range." << endl;
	}

	return surfaceArea;
}

double chain::tabulateSurfaceArea(UInt _residueIndex, UInt _atomIndex)
{
	double surfaceArea = 0.0;
	if (_residueIndex >= 0 && _residueIndex < itsResidues.size())
	{
		surfaceArea = itsResidues[_residueIndex]->tabulateSurfaceArea(_atomIndex);
	}
	else
	{
		cout << "ERROR in chain::tabulateSurfaceArea ... residue index out of range." << endl;
	}
	return surfaceArea;
}
	
void chain::removeIntraChainSpherePoints()
{
	for (UInt i = 0; i < itsResidues.size(); i ++)
	{
		itsResidues[i]->removeIntraResidueSpherePoints();
		for (UInt j = 0; j < itsResidues.size(); j ++)
		{
			if (i != j) itsResidues[i]->removeInterResidueSpherePoints(itsResidues[j]);
		}
	}

	return;
}

void chain::removeIntraChainSpherePoints(UInt _residue)
{
	if (_residue >=0 && _residue < itsResidues.size() )
	{
		itsResidues[_residue]->removeIntraResidueSpherePoints();
		for (UInt i = 0; i < itsResidues.size(); i ++)
		{
			if (i != _residue) itsResidues[_residue]->removeInterResidueSpherePoints(itsResidues[i]);
		}
	}
	else
	{
		cout << "ERROR in chain::removeIntraChainSpherePoints ... residue out of range" << endl;
	}

	return;
}

void chain::removeInterChainSpherePoints(chain* _other)
{
	for (UInt i = 0; i < itsResidues.size(); i ++)
	{
		for (UInt j = 0; j < _other->itsResidues.size(); j ++)
		{
			itsResidues[i]->removeInterResidueSpherePoints(_other->itsResidues[j]);
		}
	}

	return;
}

void chain::removeInterChainSpherePoints(UInt _residue, chain* _other)
{
	if (_residue >= 0 && _residue < itsResidues.size() )
	{
		for (UInt i = 0; i < _other->itsResidues.size(); i++ )
		{
			itsResidues[_residue]->removeInterResidueSpherePoints(_other->itsResidues[i]);
		}
	}
	else
	{
		cout << "ERROR in chain::removeInterChainSpherePoints ... residue index out of range." <<endl;
	}

	return;
}

residue* chain::superimposeGLY(const UInt _residue)
{
	if (_residue >= 0 && _residue < itsResidues.size() )
	{
		return itsResidues[_residue]->superimposeGLY();
	}
	else
	{
		cout << "ERROR in chain::superimposeGLY ... residue index out of range." << endl;
		exit(1);
	}
}

double chain::calculateHCA_O_hBondEnergy(chain* _other)
{
	double hBondEnergy = 0.0;
	for (UInt donor = 0; donor < itsResidues.size(); donor ++)
	{
		if (itsResidues[donor]->getType() == "GLY")
		{
			residue* tempGLY = itsResidues[donor]->superimposeGLY();
			for (UInt acceptor = 0; acceptor < _other->itsResidues.size(); acceptor ++)
			{
				if( itsResidues[donor] != _other->itsResidues[acceptor])
				{
					double temphBondEnergy = tempGLY->calculateHCA_O_hBondEnergy(_other->itsResidues[acceptor]);
					//if (temphBondEnergy < 0) cout << "found good hbond " << donor << " " << acceptor << " " << temphBondEnergy<< endl;
					hBondEnergy += temphBondEnergy;
				}
			}
			delete tempGLY;
		}
	}
	return hBondEnergy;
}

dblVec chain::getBackBoneCentroid()
{
	dblVec centroid(3);
	centroid[0] = 0.0; centroid[1] = 0.0; centroid[2] = 0.0;
	for (UInt i = 0; i < itsResidues.size(); i++)
	{
		centroid = centroid + itsResidues[i]->getBackBoneCentroid();
	}
	double size = (double)itsResidues.size();
	centroid = centroid/size;
	return centroid;
}
