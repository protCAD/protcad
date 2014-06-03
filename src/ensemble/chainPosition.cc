#include "chainPosition.h"

UInt chainPosition::howMany = 0;
typedef vector<allowedResidue>::iterator iterAR;
typedef vector<UIntVec>::iterator iterVEC;


chainPosition::chainPosition(UInt _resNum, UInt _aaType)
{	//cout << "ChainPosition constructor called:"
	//    << "chainPosition::chainPosition()" << endl;
        // Default behavior: Every amino acid in the type base
	// is allowed.  Limiting the allowed amino acids is left
	// initially to the user
	for (UInt i=0; i<residue::getDataBaseSize();i++)
	{	addResidue(i);

	}
	itsResNum = _resNum;
	setCurrentResIndex(_aaType);
	howMany++;
	itsSecondaryStructureIndex = 5; 
	itsRotamerLibIndex = itsSecondaryStructureIndex;
	// Proline - causes segfault when allowed c.s. 8/22/00
	// setResNotAllowed(14);
}

chainPosition::chainPosition(const chainPosition& _rhs)
{	for (UInt i=0; i<itsAllowedResidues.size();i++)
	{	allowedResidue temp(_rhs.itsAllowedResidues[i]);
		itsAllowedResidues.push_back(temp); 
	}
	itsResNum = _rhs.itsResNum;
	itsSecondaryStructureIndex = _rhs.itsSecondaryStructureIndex;
	itsRotamerLibIndex = _rhs.itsRotamerLibIndex;
	howMany++;
}

chainPosition::~chainPosition()
{	//cout<< "ChainPosition destructor called " << endl;
	howMany--;
}

void chainPosition::setResNotAllowed(UInt _aaType)
{	for (UInt i=0; i<itsAllowedResidues.size();i++)
	{	if (_aaType == itsAllowedResidues[i].getIdentity())
		{	UInt currentIdentity = getCurrentAllowedResIdentity();
			removeResidue(i);
			setCurrentResIndex(currentIdentity);
			break;
		}
	}
}

void chainPosition::setResAllowed(UInt _aaType)
{	for (UInt i=0; i<itsAllowedResidues.size();i++)
	{	if (_aaType == itsAllowedResidues[i].getIdentity())
		{
			return;
		}
	}
	UInt currentIdentity = getCurrentAllowedResIdentity();
	addResidue(_aaType);
	setCurrentResIndex(currentIdentity);
}

void chainPosition::setOnlyAllowedIdentity(UInt _aaType)
{	UInt counter = 0;	
	while (itsAllowedResidues.size() > 1)
	{	if (itsAllowedResidues[counter].getIdentity() != _aaType)
		{	
			UInt currentIdentity = getCurrentAllowedResIdentity();
			removeResidue(counter);
			setCurrentResIndex(currentIdentity);
		}
		else
		{	counter++;
		}
	}
	/*
	cout << "-----------------------" << endl;
	cout << "AllowedRes size = " << itsAllowedResidues.size() << endl;
	cout << "Identity " << itsAllowedResidues[0].getIdentity() << endl;
	cout << "Target Identity " << _aaType << endl;
	itsAllowedResidues[0].listAllowedRotamers();
	cout << "-----------------------" << endl;
	*/
}

UInt chainPosition::getCurrentAllowedResIdentity() const
{	return itsAllowedResidues[itsCurrentAllowedResIndex].getIdentity();
}

UIntVec chainPosition::getAllowedRotamers(UInt _aaType, UInt _bpt)
{
	vector <UIntVec> allowedRotamers;
	allowedRotamers.resize(0);
	for (UInt i=0; i<itsAllowedResidues.size(); i++)
    {   
		if (_aaType == itsAllowedResidues[i].getIdentity())
        {  
			allowedRotamers = itsAllowedResidues[i].getAllowedRotamers();
        }
    }
	if (_bpt < allowedRotamers.size()) return allowedRotamers[_bpt];
	else allowedRotamers.resize(0);
	return allowedRotamers[0];
}
	
		
	


void chainPosition::setRotamerNotAllowed(UInt _aaType, UInt _bpt, UInt _rot)
{	for (UInt i=0; i<itsAllowedResidues.size(); i++)
	{	if (_aaType == itsAllowedResidues[i].getIdentity())
		{	itsAllowedResidues[i].setRotamerNotAllowed(_bpt,_rot);
			break;
		}
	}
}

void chainPosition::setRotamerNotAllowed(UInt _aaType, UInt _rot)
{	for (UInt i=0; i<itsAllowedResidues.size(); i++)
	{	if (_aaType == itsAllowedResidues[i].getIdentity())
		{	itsAllowedResidues[i].setRotamerNotAllowed(0,_rot);
			break;
		}
	}
}

bool chainPosition::residueIsAllowed(UInt _aaType)
{	for (UInt i=0; i<itsAllowedResidues.size(); i++)
	{	UInt allowedIdentity = itsAllowedResidues[i].getIdentity();	
		if (_aaType == allowedIdentity)
		{	return true;
		}
	}
	return false;
}

int chainPosition::chooseResidueFromLibrary(ran& _ran)
{	// NOTE: random number input should be a uniform deviate
	// ALSO: the way this function is currently written, residues will be
	// chosen with equal likelyhood from itsAllowedResidues
	UInt numAllowed = itsAllowedResidues.size();
	if(numAllowed > 0)
	{	UInt random = UInt(_ran.getNext() * numAllowed);
		return itsAllowedResidues[random].itsIdentity;
	}
	return -1;
}

int chainPosition::chooseRotamerFromLibrary(UInt _bpt, ran& _ran)
{	
	return itsAllowedResidues[itsCurrentAllowedResIndex].chooseRotamerFromLibrary(_bpt,_ran);
}

int chainPosition::chooseBranchpoint(ran& _ran)
{	return itsAllowedResidues[itsCurrentAllowedResIndex].chooseBranchpoint(_ran);
}

vector<UInt> chainPosition::getCurrentRotamerIndex() const
{	return itsAllowedResidues[itsCurrentAllowedResIndex].getCurrentRotamerIndex();
}

void chainPosition::setCurrentRotamerIndex(const UInt _rotamer)
{	itsAllowedResidues[itsCurrentAllowedResIndex].setCurrentRotamerIndex(_rotamer);
}

void chainPosition::setCurrentRotamerIndex(vector<UInt> _rotamer)
{	itsAllowedResidues[itsCurrentAllowedResIndex].setCurrentRotamerIndex(_rotamer);
}

void chainPosition::setCurrentPolarHRotamerIndex(const UInt _rotamer)
{	itsAllowedResidues[itsCurrentAllowedResIndex].setCurrentPolarHRotamerIndex(_rotamer);
}

void chainPosition::setCurrentResIndex(const UInt _aaType)
{	for (UInt i=0; i<itsAllowedResidues.size(); i++)
	{	if ( _aaType == itsAllowedResidues[i].getIdentity())
		{	itsCurrentAllowedResIndex = i;
			return;
		}
	}
	cout << "Error in chainPosition::setCurrentResIndex()" << endl;
}

void chainPosition::removeResidue(UInt _index)
{	iterAR firstAR;
	firstAR = itsAllowedResidues.begin();
	if ( _index < itsAllowedResidues.size())
	{	itsAllowedResidues.erase(firstAR + _index);
	}
}

void chainPosition::addResidue(UInt _index)
{	
	allowedResidue tempAllowedRes(_index);
	itsAllowedResidues.push_back(tempAllowedRes);
}


void chainPosition::listAllowedRotamers()
{	cout << "Allowed rotamers for residue #" << itsResNum << endl << endl;
	for (UInt i=0; i<itsAllowedResidues.size(); i++)
	{	cout << "AATYPE " << itsAllowedResidues[i].getIdentity() << endl;
		itsAllowedResidues[i].listAllowedRotamers();
	}
	cout << "End of allowed Rotamers" << endl;
}

void chainPosition::setRotamerLibIndex(const UInt _rotamerLibIndex)
{
	if (_rotamerLibIndex == itsRotamerLibIndex)
	{	return;
	}

	for (UInt i=0; i<itsAllowedResidues.size(); i++)
	{
		itsAllowedResidues[i].buildAllowedRotamers(_rotamerLibIndex);
	}
	return;
}

UIntVec chainPosition::getResAllowed()
{
	UIntVec allowedResidues;
	allowedResidues.resize(0);
	
	for (UInt i = 0; i < itsAllowedResidues.size(); i++)
	{
		allowedResidues.push_back(itsAllowedResidues[i].getIdentity());
	}
	
	return allowedResidues;
}


vector< vector < UIntVec > > chainPosition::getAllowedDB() const
{
	vector < vector < UIntVec > > database;
	database.resize(0);
	vector < UIntVec > allowedRot;
	UInt numallowed = itsAllowedResidues.size();
	UInt maxNumAA = residueTemplate::getHowMany();
	for (UInt j=0; j<maxNumAA; j++)
	{	// see if we can find "j" residue type in the itsAllowedResidues vector
		UInt tempindex = 999;
		for (UInt i=0; i<numallowed; i++)
		{	if (itsAllowedResidues[i].getIdentity() == j)
			{	tempindex = i;
				break;
			}
		}
		if (tempindex != 999)
		{	allowedRot.resize(0);			
			allowedRot = itsAllowedResidues[tempindex].getAllowedRotamers();
			database.push_back(allowedRot);
		}
		else
		{	allowedRot.resize(0);
			database.push_back(allowedRot);
		}
	}
	return database;
}
/*   CONFORMATION CODE NOT USED
void chainPosition::listAllowedConformations(UInt _aaType)
{
	for (UInt i = 0; i < itsAllowedResidues.size(); i++)
	{
		if (itsAllowedResidues[i].getIdentity() == _aaType)
			itsAllowedResidues[i].listAllowedConformations();
	}
	return;
}


void chainPosition::setConformationNotAllowed(UInt _aaType, UIntVec _conformation)
{
	for (UInt i = 0; i < itsAllowedResidues.size(); i++)
	{
		if (itsAllowedResidues[i].getIdentity() == _aaType)
			itsAllowedResidues[i].setConformationNotAllowed(_conformation);
	}
}
*/
