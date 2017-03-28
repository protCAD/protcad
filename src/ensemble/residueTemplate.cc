#include "residueTemplate.h"

amberVDW residueTemplate::itsAmberVDW(0);
amberElec residueTemplate::itsAmberElec(0);
aaBaseline residueTemplate::itsAABaseline(0);
UInt residueTemplate::howManyTemplates = 0;
helixPropensity residueTemplate::itsHelixPropensity;
residueTemplate::residueTemplate()
{
#ifdef RESIDUE_TEMPLATE_DEBUG
	cout << "residueTemplate constructor called" << endl;
#endif
	initialize();
	itsRotamerLibs.resize(0);
	howManyTemplates++;
}

residueTemplate::residueTemplate(const residueTemplate& _rhs)
{	
		typeString = _rhs.typeString;
		typeIndex = _rhs.typeIndex;
		atomNameList = _rhs.atomNameList;
		connectivity = _rhs.connectivity;
		mainChain = _rhs.mainChain;
		branchPoints = _rhs.branchPoints;
		isMainChain = _rhs.isMainChain;
		atomList = _rhs.atomList;
		itsAtomEnergyTypeIndex = _rhs.itsAtomEnergyTypeIndex;	
		itsRotamerLibs = _rhs.itsRotamerLibs;
		chiDefinitions = _rhs.chiDefinitions;
                hasPolarHRotamers = _rhs.hasPolarHRotamers;		
		howManyTemplates++;
}

residueTemplate::~residueTemplate()
{
		howManyTemplates--;
}

void residueTemplate::initialize()
{	
#ifdef RESIDUE_TEMPLATE_DEBUG
	cout << "called residueTemplate::initialize()" << endl;
#endif
	typeString = "UNK";
	typeIndex = 0;
	atomList.resize(0);
	atomNameList.resize(0);
	connectivity.resize(0);
	mainChain.resize(0);
	isMainChain.resize(0);
	branchPoints.resize(0);
	chiDefinitions.resize(0);
	chiDefinitionsInitialized = false;
	itsAtomEnergyTypeIndex.resize(0);
	itsBondingPattern.resize(0);
        hasPolarHRotamers = false;
}

void residueTemplate::reset()
{	initialize();
	for(UInt i=0; i<itsRotamerLibs.size(); i++)
	{	if(itsRotamerLibs[i] != 0)
		{	delete itsRotamerLibs[i];
		}
	}
	itsRotamerLibs.resize(0);
}

int residueTemplate::getAtomIndexOf(const string& _name) const
{	for(UInt i=0; i<atomNameList.size(); i++)
	{	if(_name == atomNameList[i])
		{	return i;
		}
	}
	// no fit
	cout << "Invalid Atom Name Found: " << _name << endl;
	return -1;
}
	
void residueTemplate::initializeChiDefinitions()
{	if(chiDefinitions.size() != branchPoints.size())
	{	chiDefinitions.resize(branchPoints.size());
		for(UInt i=0; i<chiDefinitions.size(); i++)
		{	chiDefinitions[i].resize(0);
		}
	}
	//cout << "in residueTemplate::initializeChiDefinitions()" << endl;
	//cout << "chiDefinitions.size()= " << chiDefinitions.size() << endl;
	//for (UInt i=0; i<chiDefinitions.size(); i++)
	//{	cout << i << " size = " << chiDefinitions[i].size() << endl;
	//}
	chiDefinitionsInitialized = true;
}

void residueTemplate::initializeHasPolarHRotamers()
{
	switch (typeIndex)
	{
		case 15: // Ser
			hasPolarHRotamers=true;
			break;
		case 16: // Thr
			hasPolarHRotamers=true;
			break;
		case 18: // Tyr
			hasPolarHRotamers=false;
			break;
		default: hasPolarHRotamers=false;
	}
}

void residueTemplate::addChiDefinitions(const StrVec& _strVect)
{	if(!chiDefinitionsInitialized)
	{	initializeChiDefinitions();
	}

	int bpt = 0;
	// is the first string an unsigned number?
	//cout << "size of strvec = " << _strVect.size() << endl;
	//for (UInt i=0; i< _strVect.size(); i++)
	//{	cout << _strVect[i] << " ";
	//}
	//cout << endl;
	for(UInt i=0; i<_strVect[0].size(); i++)
	{	if(_strVect[0][i] < '0' || _strVect[0][i] > '9')
		{	bpt = -1;
			break;
		}
	}

	// if the first string is indeed an unsigned number
	// by construction, the number is the branch point from
	// which rotamers extended
	if(bpt != -1)
	{	
		sscanf(_strVect[0].c_str(), "%i", &bpt);
	}

	// if the first string is an unsigned number
	// we should have at least five strings
	// else we should have at least four strings
	// because four atoms are needed to define a chi
	if((bpt == -1 && _strVect.size() < 4) ||
		(bpt != -1 && (bpt > (int)chiDefinitions.size() || _strVect.size() < 5)))
	{	cout << typeString << " : Invalid chi definition format" << endl;
		cout << "Ignored" << endl;
		return;
	}
	
	// depending on whether there is a branch point definition
	// the position of the first string that defines an atom
	// is different
	UInt startIndex = 0;
	// if the first string is a branchpoints index
	//cout << "bpt = " << bpt << endl;
	if(bpt != -1)
	{	startIndex = 1;
		// index in the data file starts from 1
		bpt--;
	}
	// else by construction it means the branchpoint is 0
	else
	{	bpt = 0;
	}

	// chi definition is stored in a compact way, taking advantage
	// of the chi's that extend from a single branchpoint have
	// overlapping atoms
	int tempInt = 0;
	for(UInt i=startIndex; i<_strVect.size(); i++)
	{	tempInt = getAtomIndexOf(_strVect[i]);
		if(tempInt == -1)
		{	cout << typeString << ": Invalid chi input encountered !!!" << endl;
			cout << "Undefined atom name found : " << _strVect[i] << endl;
			exit(1);
		}
		chiDefinitions[bpt].push_back(tempInt);
	}

}

bool residueTemplate::chiDefinitionsNonempty() const
{	for(UInt i=0; i<chiDefinitions.size(); i++)
	{	if(chiDefinitions[i].size() >= 4)
		{	return true;
		}
	}
	return false;
}

bool residueTemplate::isValidChi(const UInt _bpt, const UInt _index) const
{	if(_bpt < chiDefinitions.size())
	{	if(chiDefinitions[_bpt].size() >= 4)
		{	if(_index + 3 < chiDefinitions[_bpt].size())
			{	return true;
			}
			else return false;
		}
		else return false;
	}
	else return false;
}

UIntVec residueTemplate::getAtomsOfChi(const UInt _bpt, const UInt _index) const
{	UIntVec quad(4);
	if(!isValidChi(_bpt, _index))
	{	return quad;
	}
	for(UInt i=0; i<4; i++)
	{	quad[i] = chiDefinitions[_bpt][_index + i];
	}
	return quad;
}

UIntVec residueTemplate::getAtomsOfPolarHChi() const
{
	UIntVec quad(4,0);
	if(getHasPolarHRotamers())
	{
		switch (typeIndex)
		{
			case 15: // Ser
			{	quad[0] = 1; //Calpha
				quad[1] = 4; //Cbeta
				quad[2] = 5; //Ogamma
				quad[3] = 10; //Hgamma
			}
			break;
			case 16: // Thr
			{	quad[0] = 1; //Calpha
				quad[1] = 4; //Cbeta
				quad[2] = 5; //Ogamma1
				quad[3] = 10; //Hgamma1
			}
			break;
			case 18: // Tyr
			{	quad[0] = 8; //CE1
				quad[1] = 10; //CZ
				quad[2] = 11; //OH
				quad[3] = 20; //HH
			}
			break;
			default: // just so it doesn't crash
			{	quad[0]=0;
				quad[1]=1;
				quad[2]=2;
				quad[3]=3;
			}
		}
	}
	return quad;
}

void residueTemplate::setHasPolarHRotamers(const bool _hasPolarHRotamers)
{
	hasPolarHRotamers = _hasPolarHRotamers;
}

bool residueTemplate::getHasPolarHRotamers() const
{
	return hasPolarHRotamers;
}

/*
void residueTemplate::assignAtomEnergyTypes()
{	for(UInt i=0; i<atomNameList.size(); i++)
	{	itsAtomEnergyTypeIndex.push_back(
			itsAtomEnergyType.getAtomEnergyTypeIndex(typeString, atomNameList[i]));
	}
} 
*/
void residueTemplate::addAtomTypeDefinitions(const StrVec& _strVect)
{	
#ifdef ATOM_TYPE_DEBUG
	//cout << "Output from residueTemplate::addAtomTypeDefinitions" << endl;
	for (UInt i=0; i< _strVect.size(); i++)
	{	//cout << _strVect[i] << " " ;
	}
	//cout << endl;
#endif
	itsAtomEnergyTypeDefinitionStrings.push_back(_strVect);
	convertAtomTypeStringsToIndices(_strVect);
}

void residueTemplate::convertAtomTypeStringsToIndices(const StrVec& _strVect)
{
#ifdef ATOM_TYPE_DEBUG
	//cout << "Output from residueTemplate::convertAtomTypeStringsToIndices" << endl;
//	cout << "  Atom name list: ";
	for (UInt i=0; i< atomNameList.size(); i++)
	{
	//	cout << atomNameList[i] << " ";
	}
	//cout << endl;
	//cout << "Residue name is : " <<  getName() << endl;
#endif
	for (UInt i=0; i< atomNameList.size(); i++)
	{	if (atomNameList[i] == _strVect[0])
		{	// we found the atom name, now parse the type
			vector< int > tempVector;
			int amberAllAtomInt = itsAmberVDW.getIndexFromNameString(_strVect[1]);
			tempVector.push_back(amberAllAtomInt);
// This is where the atom types are actually set
#ifdef ATOM_TYPE_DEBUG
		//	cout << _strVect[0] << " " << _strVect[1] << " : " << amberAllAtomInt;
#endif
			int amberUnitedAtomInt = itsAmberVDW.getIndexFromNameString(_strVect[2]);
			tempVector.push_back(amberUnitedAtomInt);
#ifdef ATOM_TYPE_DEBUG
		//	cout << " " << _strVect[2] << " : " << amberUnitedAtomInt;
#endif
			itsAtomEnergyTypeDefinitions.push_back(tempVector);
			return;
		}
	}
#ifdef ATOM_TYPE_DEBUG
	//cout << endl;
#endif
}

double residueTemplate::getAmberElecEnergy (const int _resType1, const int _atomType1, const int _resType2, const int _atomType2, const double _distance)
{
	if (_resType1 >= 0 && _atomType1 >= 0 && _resType2 >= 0 && _atomType2 >= 0)
	{
		return itsAmberElec.getEnergy(UInt(_resType1), UInt(_atomType1), UInt(_resType2), UInt(_atomType2), _distance);
	}
	return 0.0;
}

double residueTemplate::getAmberElecSoluteEnergy (const int _resType1, const int _atomType1, const int _resType2, const int _atomType2, const double _distance, const double _dielectric)
{
	if (_resType1 >= 0 && _atomType1 >= 0 && _resType2 >= 0 && _atomType2 >= 0)
	{
		return itsAmberElec.getSoluteEnergy(UInt(_resType1), UInt(_atomType1), UInt(_resType2), UInt(_atomType2), _distance, _dielectric);
	}
	return 0.0;
}

double residueTemplate::getAmberElecSoluteEnergySQ (const int _resType1, const int _atomType1, const int _resType2, const int _atomType2, const double _distanceSquared,  const double _dielectric)
{
	if (_resType1 >= 0 && _atomType1 >= 0 && _resType2 >= 0 && _atomType2 >= 0)
	{
		return itsAmberElec.getSoluteEnergySQ(UInt(_resType1), UInt(_atomType1), UInt(_resType2), UInt(_atomType2), _distanceSquared, _dielectric);
	}
	return 0.0;
}

double residueTemplate::getAmberElecEnergySQ (const int _resType1, const int _atomType1, const int _resType2, const int _atomType2, const double _distanceSquared)
{
    if (_resType1 >= 0 && _atomType1 >= 0 && _resType2 >= 0 && _atomType2 >= 0)
    {
        return itsAmberElec.getEnergySQ(UInt(_resType1), UInt(_atomType1), UInt(_resType2), UInt(_atomType2), _distanceSquared);
    }
    return 0.0;
}
double residueTemplate::getVDWRadius(const int _type1)
{
	if( _type1 >= 0)
	{
		return itsAmberVDW.getRadius(UInt(_type1));
	}
	return 0.0;
}

double residueTemplate::getPolarizability(const int _type1)
{
    if( _type1 >= 0)
    {
        return itsAmberVDW.getPolarizability(UInt(_type1));
    }
    return 0.0;
}

double residueTemplate::getVolume(const int _type1)
{
    if( _type1 >= 0)
    {
        return itsAmberVDW.getVolume(UInt(_type1));
    }
    return 0.0;
}

double residueTemplate::getVDWEnergy(const int _type1, const int _type2, const double _distance)
{
	if( _type1 >= 0 && _type2 >= 0)
		{	//cout << _type1 << " " << _type2 << endl;
			return itsAmberVDW.getEnergy(UInt(_type1),UInt(_type2), _distance);
		}
	return 0.0;
}

bool residueTemplate::isClash(const int _type1, const int _type2, const double _distance)
{
	if( _type1 >= 0 && _type2 >= 0)
                {       //cout << _type1 << " " << _type2 << endl;
                        return itsAmberVDW.isClash(UInt(_type1),UInt(_type2), _distance);
                }
        return false;
}


double residueTemplate::getVDWEnergySQ(const int _type1, const int _type2, const double _distanceSquared)
{
	if (_type1 >= 0 && _type2 >=0)
	{
		return itsAmberVDW.getEnergySQ(UInt(_type1), UInt(_type2), _distanceSquared);
	}
    else
    {
        cout << "VDW types not found in database: " << _type1 << " " << _type2 << endl;
        return 0.0;
    }

}

double residueTemplate::getVDWWaterEnergy(const int _type1)
{
    if (_type1 >= 0)
    {
        return itsAmberVDW.getWaterEnergy(UInt(_type1));
    }
    else
    {
        cout << "VDW types not found in database: " << _type1 << endl;
        return 0.0;
    }

}

double residueTemplate::getAABaselineEnergy(const string& _name)
{
	return itsAABaseline.getEnergy(_name);
}

vector<string> residueTemplate::getAABaselineList()
{
	return itsAABaseline.list();
}

int residueTemplate::getAtomEnergyTypeDefinition(const int _type, const int _field) const
{
	if (_type >=0 && _field >=0)
	{	if (UInt(_type) < itsAtomEnergyTypeDefinitions.size())
		{	if (UInt(_field) < itsAtomEnergyTypeDefinitions[_type].size())
			{	return itsAtomEnergyTypeDefinitions[_type][_field];
			}
		}
	}	
	return -1;
}

void residueTemplate::printAtomEnergyTypeDefinitions() const
{
	cout << "residueTemplate::printAtomEnergyTypeDefinitions()" << endl;
	for (UInt i=0; i< itsAtomEnergyTypeDefinitions.size(); i++)
	{	for (UInt j=0; j< itsAtomEnergyTypeDefinitions[i].size(); j++)
		{
			cout << itsAtomEnergyTypeDefinitions[i][j] << "  ";
		}
		cout << endl;
	}	
	cout << endl;
}

/************************************************************************/
/************************************************************************/
//	END OF THE FILE
/************************************************************************/
/************************************************************************/
