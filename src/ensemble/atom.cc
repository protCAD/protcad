// ***********************************************************************
// ***********************************************************************
// 	filename: atom.cc
// 	contents: class atom implementation
// 	Last Modified: 8/12/2002 by Jeff Kearns
// ***********************************************************************
// ***********************************************************************


#include "atom.h"
#include "residue.h"
#include "point.h"
#include "unitSphere.h"


// ***********************************************************************
// ***********************************************************************
// 	static class member initialization
// ***********************************************************************
// ***********************************************************************

UInt atom::howMany = 0;
bool atom::dataBaseBuilt = false;
vector<atom::typeInfo> atom::dataBase;
typedef vector <point*>::iterator iterPOINT;
double atom::itsProbeRadius = 1.4;
// ***********************************************************************
// ***********************************************************************
// 	Constructors
// ***********************************************************************
// ***********************************************************************

void atom::atomDefaultValues()
{	itsRadius = 0;
	itsEpsilon = 0;
	itsPolarizability = 0;
	itsVolume = 0;
	for( int i=0;i<itsCoords.dim();i++) itsCoords[i] = 0.0;
	itsType = 0; // 0 means "C"
	itsName = -1; // 0 means "UNK"
	itsSerialNumber = 0;
	// indicating not set, valid serial number starts from 1
	itsOccupancy = -1.0; // indicating not set with real data
	itsChainID = ' ';
	itsLigChainID= "UNK";
	itsTempFactor = -1.0;// indicating not set with real data
	itsCharge = 0; // starting with neutral charge
	itsSolvationEnergy = 0.0;
	itsDielectric = 1.0;
	itsWaters = 0;
	itsMaxDielectric = 80.4;
	itsMinDielectric = 2.25;
	itsResType = -1; // 0 means "UNK"
	itsAtomEnergyType = -1;
	generateNewSpherePoints();
	isSilent = false;
	fullySpecified = false;
	hetatmFlag=false;
}

atom::atom()
{	atomDefaultValues();
	isSilent = false;
	howMany++;
}

atom::atom(const string& _atomType)
{	atomDefaultValues();
	defineType( _atomType );
#ifdef __ATOM_DEBUG
	cout << "Atom constructor called: "
             << "atom::atom(const string& _atomType) "
             << endl;
#endif
    itsRadius = dataBase[itsType].vdwRadius[0];
	isSilent = false;
	itsEpsilon = dataBase[itsType].vdwRadius[1];
	howMany++;
}

atom::atom(const UInt _atomType)
{       atomDefaultValues();
        defineType( _atomType );
#ifdef __ATOM_DEBUG
	cout << "Atom constructor called: "
	     << "atom::atom(const UInt ) "
	     << endl;
#endif
    itsRadius = dataBase[itsType].vdwRadius[0];
	itsEpsilon = dataBase[itsType].vdwRadius[1];
	howMany++;
}

atom::atom(const pdbAtom& _pdbAtomData)
{	//atomDefaultValues();
	defineType(_pdbAtomData.getItem(element));
#ifdef __ATOM_DEBUG
	cout << "Atom constructor called: "
	     << "atom::atom(const pdbAtom&) "
	     << endl;
#endif
	itsCoords = _pdbAtomData.getAtomCoord();
    itsRadius = dataBase[itsType].vdwRadius[0];
	itsEpsilon = dataBase[itsType].vdwRadius[1];
	itsSerialNumber = _pdbAtomData.getSerial();
	itsOccupancy = _pdbAtomData.getOccupancy();
	itsTempFactor = _pdbAtomData.getTempFactor();
	itsCharge = _pdbAtomData.getCharge();
	itsChainID = _pdbAtomData.getChainID();
	itsSolvationEnergy = 0.0;
	itsDielectric = 1.0;
    itsWaters = 0;
	itsMaxDielectric = 80.4;
	itsMinDielectric = 2.25;
	howMany++;
	hetatmFlag=false;
	generateNewSpherePoints();
}

atom::atom(const PDBAtomRecord& _theRecord, bool _hetflag)
{	//atomDefaultValues();
        

	hetatmFlag=_hetflag;
	defineType(_theRecord.getElement());
#ifdef __ATOM_DEBUG
        cout << "Atom constructor called: "
             << "atom::atom(const PDBAtomRecord&) "
             << endl;
#endif
        itsCoords = _theRecord.getAtomCoord();
        itsRadius = dataBase[itsType].vdwRadius[0];
        itsEpsilon = dataBase[itsType].vdwRadius[1];
        itsSerialNumber = _theRecord.getSerial();
        itsOccupancy = _theRecord.getOccupancy();
        itsTempFactor = _theRecord.getTempFactor();
        //itsCharge = _theRecord.getCharge();
        itsCharge = 0.0;
        itsSolvationEnergy = 0.0;
        itsDielectric = 1.0;
        itsWaters = 0;
	   itsMaxDielectric = 80.4;
	   itsMinDielectric = 2.25;
	if(hetatmFlag){
		itsLigChainID=_theRecord.getChainID();
		const char* pTheChainID  = (_theRecord.getChainID()).c_str();
       		char theChainID = *pTheChainID;
        	itsChainID = theChainID;
	}	
	else{
		const char* pTheChainID  = (_theRecord.getChainID()).c_str();
       		char theChainID = *pTheChainID;
        	itsChainID = theChainID;
	}
        isSilent = false;
        fullySpecified = false;
        howMany++;
	
	if(hetatmFlag){
		itsLigName=_theRecord.getAtomName();
		itsResName=_theRecord.getResName();
	}
	
	//cout << "Name= " << itsLigName << " Serial= " << itsSerialNumber << endl;
        generateNewSpherePoints();
}

atom::atom(const PDBAtomRecord& _theRecord)
{	//atomDefaultValues();
	defineType(_theRecord.getElement());
#ifdef __ATOM_DEBUG
	cout << "Atom constructor called: "
	     << "atom::atom(const PDBAtomRecord&) "
	     << endl;
#endif
	itsCoords = _theRecord.getAtomCoord();
    itsRadius = dataBase[itsType].vdwRadius[0];
	itsEpsilon = dataBase[itsType].vdwRadius[1];
	itsSerialNumber = _theRecord.getSerial();
	itsOccupancy = _theRecord.getOccupancy();
	itsTempFactor = _theRecord.getTempFactor();
	//itsCharge = _theRecord.getCharge();
	itsCharge = 0.0;
	itsSolvationEnergy = 0.0;
	itsDielectric = 1.0;
    itsWaters = 0;
	itsMaxDielectric = 80.4;
	itsMinDielectric = 2.25;
	const char* pTheChainID  = (_theRecord.getChainID()).c_str();
	char theChainID = *pTheChainID;
	itsChainID = theChainID;
	isSilent = false;
	fullySpecified = false;
	howMany++;
	hetatmFlag=false;
	generateNewSpherePoints();
}

// deep copy
atom::atom(const atom& _rhs)
{	
	hetatmFlag=_rhs.hetatmFlag;
	if(hetatmFlag){
		itsLigName=_rhs.itsLigName;
		itsResName=_rhs.itsResName;
	}
	else{
		//cout <<"Deep copying an atom of type Protein..."<<endl; 
		itsName = _rhs.itsName;
		itsResType = _rhs.itsResType;
	}

	itsRadius = _rhs.itsRadius;
	itsEpsilon = _rhs.itsEpsilon;
	setCoords( _rhs.itsCoords);
	itsSerialNumber = _rhs.itsSerialNumber;
	itsOccupancy = _rhs.itsOccupancy;
	itsTempFactor = _rhs.itsTempFactor;
	itsCharge = _rhs.itsCharge;
	itsType = _rhs.itsType;
	itsChainID = _rhs.itsChainID;
	itsSolvationEnergy = 0.0;
	itsDielectric = 1.0;
	itsWaters = 0;
	itsMaxDielectric = 80.4;
	itsMinDielectric = 2.25;
	isSilent = _rhs.isSilent;
	itsAtomEnergyType = _rhs.itsAtomEnergyType;
	fullySpecified = _rhs.fullySpecified;
	howMany++;
	itsSpherePointFlags = _rhs.itsSpherePointFlags;
}


// ***********************************************************************
// ***********************************************************************
// 	Destructor
// ***********************************************************************
// ***********************************************************************

atom::~atom() 
{
#ifdef __ATOM_DEBUG
	cout << "Atom destructor called" << endl;
#endif
	howMany--;
}


// ***********************************************************************
// ***********************************************************************
// 	Atom Properties Related Operations 
// ***********************************************************************
// ***********************************************************************

void atom::setRadius(double _radius)
{/*	if( _radius >= dataBase[itsType].vdwRadius[1] &&
	    _radius <= dataBase[itsType].vdwRadius[2] )
	{	itsRadius = _radius;
	}
	else
	{	cout << " the radius value is invalid " << endl;
		cout << " the radius is not reset " << endl;
	}
*/
}
void atom::setSolvationEnergy(double _solvationEnergy)
{
	itsSolvationEnergy = _solvationEnergy;
}
void atom::setDielectric(double _dielectric)
{	
	itsDielectric = _dielectric;
}

void atom::setNumberofWaters(double _waters)
{
    itsWaters = _waters;
}

void atom::setMaxDielectric(double _maxDielectric)
{	
	itsMaxDielectric = _maxDielectric;
}

void atom::setMinDielectric(double _minDielectric)
{	
	itsMinDielectric = _minDielectric;
}

void atom::setSerialNumber(const UInt _number)
{
       itsSerialNumber = _number;
}

void atom::setOccupancy(const double _occupancy)
{
       if (_occupancy >= 0.00 && _occupancy < 100.00)
       {
              itsOccupancy = _occupancy;
       }
       else
       {
             // cout << " _occupancy is invalid " << endl;
       }
}

void atom::setTempFactor(const double _tempFactor)
{      itsTempFactor = _tempFactor;
}          


void atom::setCharge(const double _charge)
{	itsCharge = _charge;
}


void atom::makeAtomSilent()
{
	isSilent = true;
    //cout << "chain: " << getChainID() << " " << getName() << " number: " << getSerialNumber() << " silenced." << endl;
}

// ***********************************************************************
// ***********************************************************************
// Inter-atomic Distance Related Operations
// ***********************************************************************
// ***********************************************************************
 
double atom::distance(atom* pOtherAtom) const
{	return CMath::distance(itsCoords,pOtherAtom->getCoords());
}

double atom::distance(const atom& otherAtom) const
{	return CMath::distance(itsCoords,otherAtom.getCoords());
}

double atom::distanceSquared(atom* pOtherAtom) const
{	return CMath::distanceSquared(itsCoords,pOtherAtom->getCoords());
}

double atom::distanceSquared(const atom& otherAtom) const
{	return CMath::distanceSquared(itsCoords,otherAtom.getCoords());
}

bool atom::inCutoffSQ(const atom* _pOtherAtom, double _cutoff, double _cutoffSquared)
{
    if (fabs(itsCoords[0] - _pOtherAtom->getX()) > _cutoff) return false;
    if (fabs(itsCoords[1] - _pOtherAtom->getY()) > _cutoff) return false;
    if (fabs(itsCoords[2] - _pOtherAtom->getZ()) > _cutoff) return false;
    if (CMath::distanceSquared(itsCoords, _pOtherAtom->getCoords()) > _cutoffSquared) return false;

	return true;
}

bool atom::inCutoff(const atom* _pOtherAtom, double _cutoff)
{
	if (fabs(itsCoords[0] - _pOtherAtom->getX()) > _cutoff) return false;
	if (fabs(itsCoords[1] - _pOtherAtom->getY()) > _cutoff) return false;
	if (fabs(itsCoords[2] - _pOtherAtom->getZ()) > _cutoff) return false;
	if (CMath::distance(itsCoords, _pOtherAtom->getCoords()) > _cutoff) return false;
	return true;
}

double atom::inCubeWithDist(const atom* _pOtherAtom, double _cutoff)
{
	double distance = 0.0;
	if (fabs(itsCoords[0] - _pOtherAtom->getX()) > _cutoff) return 0.0;
	if (fabs(itsCoords[1] - _pOtherAtom->getY()) > _cutoff) return 0.0;
	if (fabs(itsCoords[2] - _pOtherAtom->getZ()) > _cutoff) return 0.0;
	distance = CMath::distance(itsCoords, _pOtherAtom->getCoords());
	return distance;
}

double atom::inCubeWithDistSQ(const atom* _pOtherAtom, double _cutoff)
{
	double distance = 0.0;
    if (fabs(itsCoords[0] - _pOtherAtom->getX()) > _cutoff) return 999.0;
    if (fabs(itsCoords[1] - _pOtherAtom->getY()) > _cutoff) return 999.0;
    if (fabs(itsCoords[2] - _pOtherAtom->getZ()) > _cutoff) return 999.0;
	distance = CMath::distanceSquared(itsCoords, _pOtherAtom->getCoords());
	return distance;
}

bool atom::inCube(const atom* _pOtherAtom, double _cutoff)
{
	if (fabs(itsCoords[0] - _pOtherAtom->getX()) > _cutoff) return false;
	if (fabs(itsCoords[1] - _pOtherAtom->getY()) > _cutoff) return false;
	if (fabs(itsCoords[2] - _pOtherAtom->getZ()) > _cutoff) return false;
	return true;
}

bool atom::inCutoff(const atom& _otherAtom, double _cutoff)
{
    if (fabs(itsCoords[0] - _otherAtom.getX()) > _cutoff) return false;
    if (fabs(itsCoords[1] - _otherAtom.getY()) > _cutoff) return false;
    if (fabs(itsCoords[2] - _otherAtom.getZ()) > _cutoff) return false;
    if (CMath::distance(itsCoords, _otherAtom.getCoords()) > _cutoff) return false;
    return true;
}



// ***********************************************************************
// ***********************************************************************
// Atom Type Related Operations
// ***********************************************************************
// ***********************************************************************

string atom::getType() const
{	if( itsType < dataBase.size() )
	{	return dataBase[itsType].typeName;
	}
	else
	{	cout << "Error: itsType: [" << itsType << "]"
                     << " incompatible with dataBase "
		     << endl;
		cout << "Error reported by: atom::getType() "
		     << endl;
		return "";
	}
}

void atom::setType(const string& _atomType)
{	for(UInt i=0 ; i<dataBase.size() ; i++)
	{	if( _atomType == dataBase[i].typeName )
		{	itsType=i;
			return;
		}
	}
	cout << "Error: _atomType: [" << _atomType << "]"
	     << " not found in the type list " << endl
	     << " Atom's type is not set or changed "
	     << endl;
	cout << "Error reported by: atom::setType(const string&) "
	     << endl;
}

void atom::setType(const UInt _atomType)
{	if( _atomType <= dataBase.size() )
	{	itsType = _atomType;
	}
	else
	{	cout << " not found in type list " << endl
		     << " Atom's type is not set or changed "
		     << endl;
		cout << "Error reported by: atom::setType(const UInt ) "
		     << endl;
	}
}

void atom::defineType(const UInt _atomType)
{	if( !dataBaseBuilt )
	{	buildDataBase();
	}

	setType( _atomType );
}

void atom::defineType(const string& _atomType)
{	if(!dataBaseBuilt)
	{	buildDataBase();
	}

	setType( _atomType );
}

void atom::buildDataBase()
{	
	string evname = "PROTCADDIR";
	string path = getEnvironmentVariable(evname);

	path += "/data/";
	string filebase="ATOMTYPE";
	
	string iFile = path + filebase;
	ifstream dataBaseFile( iFile.c_str() );
	dataBase.resize(0);

	if( ! dataBaseFile )
	{	cout << "Error: unable to open input file: "
		     << iFile << endl;
		cout << "Error reported by: atom::buildDataBase() " << endl;
		exit (1);
	}

	atom::typeInfo typeBuf;
	string tempString;
	double tempDouble;

	for(;;)
	{	if( dataBaseFile >> tempString )
		{	typeBuf.typeName = tempString;
		}
		else
		{	break;
		}

		for(UInt i=0;i<typeBuf.vdwRadius.size();i++)
		{	if( dataBaseFile >> tempDouble )
			{	typeBuf.vdwRadius[i] = tempDouble;
			}
			else
			{	cout << "Error: file: ATOMTYPE is corrupted, exit forced ";
				cout << endl;
				cout << "Error reported by: atom::buildDataBase() ";
				cout << endl;
				exit (1);
			}
		}

		dataBase.push_back(typeBuf);
	}

	dataBaseBuilt = true;

#ifdef __ATOM_DEBUG
	cout << " The Atom Type List : " << endl;
	
	for(UInt i = 0 ; i<dataBase.size() ; i++)
	{	cout<< i << "\t:\t" << dataBase[i].typeName << "\t"
			<< dataBase[i].vdwRadius[0] << "\t"
			<< dataBase[i].vdwRadius[1] << "\t"
			<< endl;
	}
#endif

}

// ********************************************************************
// ********************************************************************
// Atom Name Operations
// ********************************************************************
// ********************************************************************

string atom::getName() const
{
	string tempstring;
	
	if(hetatmFlag){
		return itsLigName;
	}
	else{
		if (itsName >=0 && itsName < (int)residue::getAtomNameBaseSize(itsResType))
			tempstring = residue::getAtomNameBaseItem(itsResType,itsName);
	
		else
			tempstring = "UNK";
	}
	return tempstring;
}

// ********************************************************************
// ********************************************************************
// Linkage to the residue it belongs to
// ********************************************************************
// ********************************************************************

string atom::getResType() const
{	
	if(hetatmFlag){
		//cout << "in getResType..."<< endl;
		//cout << "itsResName= " << itsResName<< endl;

		return itsResName;
	}
	else{
		return residue::getDataBaseItem( itsResType );
	}
}


// spherepoint and surface area operations

void atom::generateNewSpherePoints()
{
	for (UInt i = 0; i < spherePoint::getSphereSize(); i++)
	{
		itsSpherePointFlags.push_back(true);
	}
	return;
}

void atom::initializeSpherePoints()
{
	for (UInt i = 0; i < itsSpherePointFlags.size(); i++)
	{
		itsSpherePointFlags[i] = true;
	}
	return;
}

void atom::removeSpherePoint(UInt _index)
{
	if (_index >=0 && _index < itsSpherePointFlags.size())
	{
		itsSpherePointFlags[_index] = false;
	}
	return;
}

double atom::calculateTotalSASA()
{
    double radius = itsRadius + itsProbeRadius;
    double surfaceArea = 4 * 3.1415927 * radius * radius;

    return surfaceArea;
}

double atom::calculateExposedSASA()
{
    UInt totalNumberOfPoints = spherePoint::getSphereSize();
    UInt numberOfExistingPoints = countTrueSpherePoints();

    if (totalNumberOfPoints == 0)
    {
        cout << "ERROR in calculateExposedSASA ... division by zero!" << endl;
        cout << "zot!  this means that the unitSphere point size is zero ... " << endl;
        return -1.0;
    }
    //cout << numberOfExistingPoints << " " << totalNumberOfPoints << " ";
    double ratioExposed = (double)numberOfExistingPoints / (double)totalNumberOfPoints;
    double exposedSASA = ratioExposed * this->calculateTotalSASA();

    return exposedSASA;
}

double atom::calculateBuriedSASA()
{
    UInt totalNumberOfPoints = spherePoint::getSphereSize();
    UInt numberOfExistingPoints = countTrueSpherePoints();

    if (totalNumberOfPoints == 0)
    {
        cout << "ERROR in calculateBuriedSASA ... division by zero!" << endl;
        cout << "wuhhhrg!  this means that the unitSphere point size is zero ... " << endl;
        return -1.0;
    }

    double ratioExposed = (double)numberOfExistingPoints / (double)totalNumberOfPoints;
    double buriedSASA = (1-ratioExposed) * this->calculateTotalSASA();

    return buriedSASA;
}

UInt atom::countTrueSpherePoints()
{
	UInt count = 0;

	for (UInt i = 0; i < itsSpherePointFlags.size(); i ++)
	{
		if (itsSpherePointFlags[i] == true)
		{
			count ++;
		}
	}
	return count;
}

void atom::removeSpherePoints(atom* _pOtherAtom)
{
    // return if no spherepoints left
    if ( countTrueSpherePoints() == 0) return;


	double thisX = itsCoords[0];
	double thisY = itsCoords[1];
	double thisZ = itsCoords[2];

	double otherX = _pOtherAtom->itsCoords[0];
	double otherY = _pOtherAtom->itsCoords[1];
	double otherZ = _pOtherAtom->itsCoords[2];

	double minSeparation = itsRadius + itsProbeRadius + _pOtherAtom->itsRadius + _pOtherAtom->itsProbeRadius;

	if (minSeparation < fabs(thisX-otherX) || minSeparation < fabs(thisY-otherY) || minSeparation < fabs(thisZ-otherZ)) return;

	for (UInt i = 0; i < itsSpherePointFlags.size(); i ++)
	{
		if (itsSpherePointFlags[i])
		{
			double sphereX = spherePoint::getX(i)*(itsRadius+itsProbeRadius)+thisX;
			double sphereY = spherePoint::getY(i)*(itsRadius+itsProbeRadius)+thisY;
			double sphereZ = spherePoint::getZ(i)*(itsRadius+itsProbeRadius)+thisZ;

			double theDistance = sqrt( (sphereX-otherX)*(sphereX-otherX) + (sphereY-otherY)*(sphereY-otherY) + (sphereZ-otherZ)*(sphereZ-otherZ) );
			if (theDistance < (_pOtherAtom->itsRadius + _pOtherAtom->itsProbeRadius) )
			{
				removeSpherePoint(i);
			}
		}
	}
    return;
}


void atom::setAtomicRadius(double _radius)
{
    itsRadius = _radius;
    return;
}


// ********************************************************************
// ********************************************************************
// END OF FILE 
// ********************************************************************
// ********************************************************************
