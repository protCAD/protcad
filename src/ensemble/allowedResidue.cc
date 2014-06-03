#include "allowedResidue.h"

typedef vector<UInt>::iterator iterINT;
typedef vector<UIntVec>::iterator iterVEC;

allowedResidue::allowedResidue()
{	//cout << "allowedResidue default constructor called:"
	//    << "allowedResidue::allowedResidue()" << endl;
	// Default behavior: Every rotamer is initially allowed
	//buildAllowedConformations(3);
}

allowedResidue::allowedResidue(const UInt _aaType)
{	itsIdentity= _aaType;
	itsCurrentRotamerIndex.push_back(0);
	// default rotamer lib is the backbone independent one
	buildAllowedRotamers(5);
	//buildAllowedConformations(3);
}

allowedResidue::allowedResidue(const UInt _aaType, const UInt _rotamerLibIndex)
{	itsIdentity = _aaType;
	itsCurrentRotamerIndex.push_back(0);
	buildAllowedRotamers(_rotamerLibIndex);
	//buildAllowedConformations(3);
}

//deep copy
allowedResidue::allowedResidue(const allowedResidue& _rhs)
{
	itsAllowedRotamers = _rhs.itsAllowedRotamers;
	itsCurrentRotamerIndex = _rhs.itsCurrentRotamerIndex;
	itsIdentity = _rhs.itsIdentity;
	//itsAllowedConformations = _rhs.itsAllowedConformations;
}

allowedResidue::~allowedResidue()
{	//cout<< "ChainPosition destructor called " << endl;
}

void allowedResidue::buildAllowedRotamers(const UInt _rotamerLibIndex)
{
	itsAllowedRotamers.resize(0);
	UInt bpt = residue::getNumBpt(itsIdentity);
	itsAllowedRotamers.resize(bpt);
	for (UInt i=0;i<bpt;i++)
	{	if ( residue::dataBase[itsIdentity].chiDefinitions.size() > 0 &&
			residue::dataBase[itsIdentity].chiDefinitions[i].size() >= 4 )
		{	UInt numRotamers= (residue::dataBase[itsIdentity].itsRotamerLibs[_rotamerLibIndex])->itsRotamers[i].size();
			UIntVec rocco;
			for (UInt j=0;j<numRotamers;j++)
			{	rocco.push_back(j);
			}
			itsAllowedRotamers[i] = rocco;
		}
	}
//	for (UInt i = 0; i < itsAllowedRotamers.size(); i++)
//	{
//		for (UInt j = 0; j < itsAllowedRotamers[i].size(); j++) cout << " " << itsAllowedRotamers[i][j];
//		cout << endl;
//	}
}

void allowedResidue::setRotamerNotAllowed(const UInt _bpt,const UInt _rot)
{	UInt numbpt = itsAllowedRotamers.size();
	if (_bpt < numbpt)
	{	for (UInt i=0; i<itsAllowedRotamers[_bpt].size(); i++)
		{	if (_rot == itsAllowedRotamers[_bpt][i])
			{	deleteRotamer(_bpt,i);
			}
		}
	}
}

UInt allowedResidue::getRotamerLibSize(const UInt _bpt) const
{	if (_bpt < itsAllowedRotamers.size())
	{	return itsAllowedRotamers[_bpt].size();
	}
	return 0;
}

int allowedResidue::chooseRotamerFromLibrary(const UInt _bpt, ran& _ran) const
{	if (itsAllowedRotamers.size() != 0 )
	{	
#ifdef _ALLOWED_RESIDUE_DEBUG
		cout << "allowedResidue::itsIdentity = " << itsIdentity << endl;
		cout << "itsAllowedRotamers.size() = " << itsAllowedRotamers.size() << endl;
#endif
		if (_bpt < itsAllowedRotamers.size())
		{	UInt numAllowed;
			numAllowed = itsAllowedRotamers[_bpt].size();
			//cout << "numAllowed = " << numAllowed << endl;
			if (numAllowed > 0)
			{	
				UInt random = UInt(_ran.getNext()*numAllowed);
#ifdef _ALLOWED_RESIDUE_DEBUG
				cout << endl << "numAllowedRotamers = " << numAllowed << endl;
				cout << "choosing " << random << endl;
#endif
				return itsAllowedRotamers[_bpt][random];
			}
		}
		else
		{	cout << "Error reported from allowedResidue::chooseRotamerFromLibrary()" << endl;
			cout << "Branchpoint out of range: " << _bpt << endl;
		}
	}
	return -1;
}

int allowedResidue::chooseBranchpoint(ran& _ran) const
{	UInt numAllowed = itsAllowedRotamers.size();
	//cout << endl << "itsAllowedRotamers.size() = " << numAllowed << endl;
	if (numAllowed > 0)
	{	return int(_ran.getNext()*numAllowed);
	}
	return -1;
}

void allowedResidue::listAllowedRotamers() const
{	for (UInt i=0; i<itsAllowedRotamers.size(); i++)
	{	cout << "Branchpoint " << i << ": ";
		for (UInt j=0; j<itsAllowedRotamers[i].size(); j++)
		{	cout << itsAllowedRotamers[i][j] << " ";
		}
		cout << endl;
	}
	cout << endl;
}

void allowedResidue::deleteRotamer(const UInt _bpt, const UInt _index)
{	if (_bpt < itsAllowedRotamers.size())
	{	if (_index < itsAllowedRotamers[_bpt].size())
		{	iterINT firstINT;
			firstINT = itsAllowedRotamers[_bpt].begin();
			itsAllowedRotamers[_bpt].erase(firstINT + _index);
		}
		else
		{	cout << "Invalid rotamer index passed to allowedResidue::deleteRotamer" << endl;
		}
	}
	else 
	{	cout << "Invalid branch passed to allowedResidue::deleteRotamer" << endl;
	}
}

int allowedResidue::getRotamerFromIndex(const UInt _bpt, const UInt _index) const
{	if (_bpt < itsAllowedRotamers.size())
	{	if (_index < itsAllowedRotamers[_bpt].size())
		{	return itsAllowedRotamers[_bpt][_index];
		}
	}
	return -1;
}

void allowedResidue::setCurrentRotamerIndex(const UInt _rotamerIndex)
{
	itsCurrentRotamerIndex[0] = _rotamerIndex;
}

void allowedResidue::setCurrentRotamerIndex(vector<UInt> _rotamerIndex)
{
	itsCurrentRotamerIndex = _rotamerIndex;
}

void allowedResidue::setCurrentPolarHRotamerIndex(const UInt _rotamerIndex)
{
	itsCurrentPolarHRotamerIndex = _rotamerIndex;
}




// CONFORMATION CODE - not needed at moment.  redundant with rotamers ... vik 3/25/01
// contains *incorrect* treatment of rotamers and branchpoints

/*void allowedResidue::listAllowedConformations()
{
	for (UInt i = 0; i < itsAllowedConformations.size(); i++)
	{
		for (UInt j = 0; j < itsAllowedConformations[i].size(); j++)
		{
			cout << itsAllowedConformations[i][j] << " ";
		}
		cout << endl;
	}
	return;
}*/

/*void allowedResidue::setConformationNotAllowed( UIntVec _conformation )
{
	bool found = false;
	int index = -1;

	for (UInt i = 0; i < itsAllowedConformations.size(); i++)
	{
		if (itsAllowedConformations[i] == _conformation)
		{
			index = (int)i;
			found = true;
		}
	}
	if (found)
	{
		vector <UIntVec> tempConformation = itsAllowedConformations;
		itsAllowedConformations.resize(0);
		for (int i = 0; i < index; i ++)
		{	
			itsAllowedConformations.push_back(tempConformation[i]);
		}
		for (int i = index + 1; i < (int)tempConformation.size(); i++)
		{
			itsAllowedConformations.push_back(tempConformation[i]);
		}
	}
	return;
}*/


/*void allowedResidue::buildAllowedConformations(UInt _rotamerModulus)
{
    itsAllowedConformations.resize(0);
    UInt count = 1;
	UInt bpts = residue::getNumBpt(itsIdentity);
    for (UInt i = 0; i < bpts; i ++) {count = count * _rotamerModulus;}
    UIntVec tempBptVec; tempBptVec.resize(0);
    for (UInt i = 0; i < bpts; i++) {tempBptVec.push_back(1);} // create initial conformation of all 1's

    for (UInt i = 0; i < count; i++)
    {
        itsAllowedConformations.push_back(tempBptVec);
		//for (UInt k = 0; k < tempBptVec.size(); k++) cout << tempBptVec[k] << " ";
		//cout << endl;
        UInt position = tempBptVec.size()-1;
        tempBptVec = incrementConf(tempBptVec, position, _rotamerModulus);
    }
    return;
}*/

/*UIntVec allowedResidue::incrementConf(UIntVec _tmpVec, UInt _position, UInt _rotamerModulus) // recursive companion function for building conformation library
{
    UInt test = _tmpVec[_position] + 1;
    if (test > _rotamerModulus)  // if incrementConf causes rollover
    {
        _tmpVec[_position] = 1;  // reset at new value
        if (_position != 0)       // if not first position
        {
            _tmpVec = incrementConf(_tmpVec, _position - 1, _rotamerModulus); // then incrementConf previous index
        }
        else                     // if indeed first position
        {
            test = _tmpVec[_position - 1] + 1; // test if first index is maxed ... shouldn't be, but it doesn't hurt to check
            if (test < _rotamerModulus)
                _tmpVec[_position -1] ++;
        }
    }
    else
    {
        _tmpVec[_position]++;
    }
    return _tmpVec;
}*/
	
