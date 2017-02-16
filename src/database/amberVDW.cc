#include "amberVDW.h"

double amberVDW::itsScaleFactor = 1.0;
double amberVDW::itsRadiusScaleFactor = 1.0;
bool amberVDW::linearRepulsionDampening = false;
double amberVDW::itsAttractionScaleFactor = 1.0;
double amberVDW::itsRepulsionScaleFactor = 1.0;

amberVDW::amberVDW()
{	// default constructor
	// read the file until the end.....
#ifdef AMBERVDW_DEBUG
	cout << "amberVDW::amberVDW() called" << endl;
#endif
	itsFileName = "amberVDW.frc";
	R_ref.resize(0);
	EPS.resize(0);
    Pol_ref.resize(0);
    Vol_ref.resize(0);
	buildDataBase();
	//cout << " amberVDW database is built " << endl;
#ifdef AMBERVDW_DEBUG
	for (UInt i=0; i< R_ref.size(); i++)
	{	cout << amberAtomTypeNames[i] << "   " << R_ref[i];
		cout << "    " << EPS[i] << endl;
	}
#endif
}

amberVDW::amberVDW(int _Dummy)
{	
#ifdef AMBERVDW_DEBUG
	cout << "amberVDW::amberVDW(int) called" << endl;
#endif
	// another constructor
	// read the file until the end.....
	itsFileName = "amberVDW.frc";
	R_ref.resize(0);
	EPS.resize(0);
    Pol_ref.resize(0);
    Vol_ref.resize(0);
	buildDataBase();
	//cout << " amberVDW database is built " << endl;
#ifdef AMBERVDW_DEBUG
	cout << "AmberVDWAtomTypeNames: " << endl;
	for (UInt i=0; i< R_ref.size(); i++)
	{	cout << i << "  "<< amberAtomTypeNames[i] << "   " << R_ref[i];
		cout << "    " << EPS[i] << endl;
	}
#endif
}

//deep copy constructor
amberVDW::amberVDW(const amberVDW& _otherAmberVDW)
{
#ifdef AMBERVDW_DEBUG
	cout << "amberVDW deep copy constructor called " << endl;
#endif
	itsFileName = _otherAmberVDW.itsFileName;
	R_ref = _otherAmberVDW.R_ref;
	EPS = _otherAmberVDW.EPS;
    Pol_ref = _otherAmberVDW.Pol_ref;
    Vol_ref = _otherAmberVDW.Vol_ref;
	amberAtomTypeNames = _otherAmberVDW.amberAtomTypeNames;
}

amberVDW::~amberVDW()
{
}

bool amberVDW::isClash(const UInt _type1, const UInt _type2, const double _distance)
{
	double R_ref_pair = itsRadiusScaleFactor * (R_ref[_type1] + R_ref[_type2]);
	if (_distance < R_ref_pair) return true;
	return false;
}

double amberVDW::getRadius(const UInt _type1)
{
    double radius  = (R_ref[_type1]) * itsRadiusScaleFactor;
    return radius;
}

double amberVDW::getPolarizability(const UInt _type1)
{
    double polarizability  = Pol_ref[_type1];
    return polarizability;
}

double amberVDW::getVolume(const UInt _type1)
{
    double volume  = Vol_ref[_type1];
    return volume;
}

double amberVDW::getEnergySQ(const UInt _type1, const UInt _type2, const double _distanceSquared) const //optimized to avoid distance calculation if possible
{
    double energy = 0.0;
    double R_ref_pair = 0.0;
    double EPS_pair = 0.0;
    if (_type1 < R_ref.size())
    {
        if (_type2 < R_ref.size())
        {       
			R_ref_pair  = itsRadiusScaleFactor * (R_ref[_type1] + R_ref[_type2]);
            if (EPS[_type1] == EPS[_type2])
				EPS_pair = EPS[_type1]; // save a sqrt operation
			else
				EPS_pair = sqrt( EPS[_type1] * EPS[_type2]);
            if (linearRepulsionDampening) // see Kuhlman & Baker PNAS v97 p10383  (2000) - online supplementary materials
			{
				double distance = sqrt(_distanceSquared);
                if ( (pow(R_ref_pair,2)/_distanceSquared) < 1.12)
                    energy = EPS_pair * (( itsRepulsionScaleFactor * pow(R_ref_pair,12)/pow(_distanceSquared,6)) - ( 2 * itsAttractionScaleFactor * pow(R_ref_pair,6)/pow(_distanceSquared,3)) );
                else
                    energy = 10 - 11.2 * (distance / R_ref_pair );
			}
            else 
			{
				//cout << "R_ref_pair used in getEnergySQ is " << R_ref_pair << endl;
				energy = EPS_pair * (( itsRepulsionScaleFactor * pow(R_ref_pair,12)/pow(_distanceSquared, 6)) - (2 * itsAttractionScaleFactor * pow(R_ref_pair,6)/pow(_distanceSquared,3)));
				//cout << "distancesquared: " << _distanceSquared << " " << energy << endl;
			}
        }
    }
    energy *= itsScaleFactor;
    return energy;
}   


double amberVDW::getEnergy(const UInt _type1, const UInt _type2, const double _distance) const
{
	double energy = 0.0;
	double R_ref_pair = 0.0;
	double EPS_pair = 0.0;
	if (_type1 < R_ref.size())
	{
		if (_type2 < R_ref.size())
		{		
			R_ref_pair  = itsRadiusScaleFactor * (R_ref[_type1] + R_ref[_type2]);
			if (EPS[_type1] == EPS[_type2])
				EPS_pair = EPS[_type1];
			else
				EPS_pair = sqrt( EPS[_type1] * EPS[_type2]);
			if (linearRepulsionDampening) // see Kuhlman & Baker PNAS v97 p10383  (2000) - online supplementary materials
				if ( (R_ref_pair/_distance) < 1.12)
					energy = EPS_pair * ( itsRepulsionScaleFactor * pow( (R_ref_pair/_distance),12) - (2 * itsAttractionScaleFactor * pow((R_ref_pair/_distance),6)));
				else
					energy = 10 - 11.2 * (_distance / R_ref_pair );
			else energy = EPS_pair * ( itsRepulsionScaleFactor * pow( (R_ref_pair/_distance),12) - (2 * itsAttractionScaleFactor * pow((R_ref_pair/_distance),6)));
			//cout << "R_ref_pair used in getEnergy is " << R_ref_pair << endl;
		}
	}
	energy *= itsScaleFactor;
	return energy;
}

double amberVDW::getWaterEnergy(const UInt _type1) const
{
    double energy = 0.0;
    double EPS_pair = 0.0;
    UInt waterType = 52;
    if (_type1 < EPS.size())
    {
        if (EPS[_type1] == EPS[waterType])
            EPS_pair = EPS[_type1];
        else
            EPS_pair = sqrt( EPS[_type1] * EPS[waterType]);
        energy = EPS_pair * ( itsRepulsionScaleFactor - (2 * itsAttractionScaleFactor));
    }
    return energy;
}

void amberVDW::buildDataBase()
{	
	string evname = "PROTCADDIR";
	string path = getEnvironmentVariable(evname);

	path += "/data/";
	string iFile = path + itsFileName;
	ifstream inFile;
	string currentLine;
	StrVec parsedStrings;
	parsedStrings.resize(0);
	
	// read the amber.frc data file
	inFile.open(iFile.c_str());
	if (!inFile)
	{	cout << "Error: unable to open input file: ";
		cout << iFile << endl;
		exit (1);
	}

	while (getline(inFile, currentLine, '\n'))
	{	// ignore the comment line
		// comment line should start with #
		if(currentLine[0] != '#' && currentLine[0] != '@'
			&& currentLine[0] != '!' && currentLine[0] != '>')
		{	
			parsedStrings=Parse::parse(currentLine);
            if (parsedStrings.size() == 7) convertToDataElements(parsedStrings);
			parsedStrings.resize(0);
		}
	}

	inFile.close();
	inFile.clear();
}

// where specific information about parsed data is intepreted
void amberVDW::convertToDataElements(const StrVec& _parsedStrings)
{
	double tmpDouble;
	amberAtomTypeNames.push_back(_parsedStrings[2]);
	sscanf(_parsedStrings[3].c_str(), "%lf", &tmpDouble);
	R_ref.push_back(tmpDouble);
	tmpDouble = 0.0;
	sscanf(_parsedStrings[4].c_str(), "%lf", &tmpDouble);
	EPS.push_back(tmpDouble);
    tmpDouble = 0.0;
    sscanf(_parsedStrings[5].c_str(), "%lf", &tmpDouble);
    Pol_ref.push_back(tmpDouble);
    tmpDouble = 0.0;
    sscanf(_parsedStrings[6].c_str(), "%lf", &tmpDouble);
    Vol_ref.push_back(tmpDouble);
}

int amberVDW::getIndexFromNameString(string _name)
{
	for (UInt i=0; i<amberAtomTypeNames.size(); i++)
	{	if (amberAtomTypeNames[i] == _name)
		{
#ifdef AMBERVDW_DEBUG
			cout <<  "In amberVDW::getIndexFromNameString(string _name)" << endl;
			cout <<  "   found atom type " << _name << " at position " << i << endl;
#endif
			return i;
		}
	}
#ifdef AMBERVDW_DEBUG
	cout << "In amberVDW::getIndexFromNameString(string _name)" << endl;
	cout << "    couldnt find atom type: " << _name << endl;
#endif
	// If we're not able to find the atom name here, warn
	// the user, and return a -1, which is out of range normally.
	return -1;
}
