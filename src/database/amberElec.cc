#include "amberElec.h"

double amberElec::itsScaleFactor = 1.0;
bool amberElec::distanceDependance = true;
double amberElec::itsDielectricConstant = 12.0;
bool amberElec::highEnergyCutOff = true;
amberElec::amberElec()
{
#ifdef AMBERELEC_DEBUG
	cout << "amberElec::amberElec() called" << endl;
#endif
	resNames.resize(0);
	charges.resize(0);
	atomNames.resize(0);
}

amberElec::amberElec(int _dummy)
{
#ifdef AMBERELEC_DEBUG
	cout << "amberElec::amberElec(int) called" << endl;
#endif
	resNames.resize(0);
	charges.resize(0);
	atomNames.resize(0);
}

amberElec::amberElec(const amberElec& _other)
{
#ifdef AMBERELEC_DEBUG
	cout << "amberElec deep copy constuctor called" << endl;
#endif
	charges = _other.charges;
	atomNames = _other.atomNames;
	resNames = _other.resNames;
	itsFileName = _other.itsFileName;
}

amberElec::~amberElec()
{
}

double amberElec::getEnergy(const UInt _resType1, const UInt _atomType1, const UInt _resType2, const UInt _atomType2, const double _distance) const
{
	double energy = 0.0;
	if (_resType1 < resNames.size() && _resType2 < resNames.size())
	{
		if (_atomType1 < atomNames[_resType1].size() && _atomType2 < atomNames[_resType2].size())
		{
			//cout << "charge 1: " <<charges[_resType1][_atomType1]  << " charge 2: " <<charges[_resType2][_atomType2] << " distance: " << _distance << endl;
			energy = KC * (charges[_resType1][_atomType1] * charges[_resType2][_atomType2]) / _distance;

			if (distanceDependance)
				energy /= 4 * _distance;
			else if (itsDielectricConstant > 0)
				energy /= itsDielectricConstant;
			else
				cout << "ERROR - dielectric constant must be positive & nonzero" << endl;

			energy *= itsScaleFactor;
			if (highEnergyCutOff && energy > 10.0) energy = 10.0;
			return energy;
		}
		else cout << "error - atom indicies out of range: r1: " << _resType1 << " r2: " << _resType2 << " a1: " << _atomType1 << " a2: " << _atomType2 << endl;
	}
	else  cout << "error - residue indicies out of range: r1: " << _resType1 << " r2: " << _resType2 << " a1: " << _atomType1 << " a2: " << _atomType2 << endl;
	return 0.0;
}

double amberElec::getSoluteEnergy(const UInt _resType1, const UInt _atomType1, const UInt _resType2, const UInt _atomType2, const double _distance, const double _dielectric) const
{
	double energy = 0.0;
	if (_resType1 < resNames.size() && _resType2 < resNames.size())
	{
		if (_atomType1 < atomNames[_resType1].size() && _atomType2 < atomNames[_resType2].size())
		{
			energy = (KC * (charges[_resType1][_atomType1] * charges[_resType2][_atomType2]) / _distance);

			if (_dielectric > 0)
			{
				energy /= _dielectric;
			}
			else
			{
				cout << "ERROR - dielectric constant must be positive & nonzero" << endl;
			}
			energy *= itsScaleFactor;
			return energy;
		}
		else cout << "error - atom indicies out of range: r1: " << _resType1 << " r2: " << _resType2 << " a1: " << _atomType1 << " a2: " << _atomType2 << endl;
	}
	else  cout << "error - residue indicies out of range: r1: " << _resType1 << " r2: " << _resType2 << " a1: " << _atomType1 << " a2: " << _atomType2 << endl;
	return 0.0;
}

double amberElec::getSoluteEnergySQ(const UInt _resType1, const UInt _atomType1, const UInt _resType2, const UInt _atomType2, const double _distanceSquared, const double _dielectric) const
{
	double energy = 0.0;
	if (_resType1 < resNames.size() && _resType2 < resNames.size())
	{
		if (_atomType1 < atomNames[_resType1].size() && _atomType2 < atomNames[_resType2].size())
		{
			energy = (KC * (charges[_resType1][_atomType1] * charges[_resType2][_atomType2]) / sqrt(_distanceSquared)) / _dielectric;
			energy *= itsScaleFactor;
			return energy;
		}
		else cout << "error - atom indicies out of range: r1: " << _resType1 << " r2: " << _resType2 << " a1: " << _atomType1 << " a2: " << _atomType2 << endl;
	}
	else  cout << "error - residue indicies out of range: r1: " << _resType1 << " r2: " << _resType2 << " a1: " << _atomType1 << " a2: " << _atomType2 << endl;
	return 0.0;
}


double amberElec::getEnergySQ(const UInt _resType1, const UInt _atomType1, const UInt _resType2, const UInt _atomType2, const double _distanceSquared) const
{
	if (distanceDependance)
	{
		double energy = itsScaleFactor * KC * (charges[_resType1][_atomType1] * charges[_resType2][_atomType2]) / (4 * _distanceSquared);
		//cout << "charge 1: " <<charges[_resType1][_atomType1]  << " charge 2: " <<charges[_resType2][_atomType2] << " distance squared: " << _distanceSquared << endl;
		if (highEnergyCutOff && energy > 10.0) energy = 10.0;
		return energy;
	}
	double distance = sqrt(_distanceSquared);
	return getEnergy(_resType1, _atomType1, _resType2, _atomType2, distance);
}


double amberElec::getItsCharge(const UInt _resType, const UInt _atomType) const
{
	if (_resType >= 0 && _resType < resNames.size())
	{
		if (_atomType >= 0 && _atomType < atomNames[_resType].size())
		{
			return charges[_resType][_atomType];
		}
		else cout << " no valid charge assignment for atom " << getItsAtomName(_resType, _atomType) << " (" << _atomType << ") in res " << resNames[_resType] << endl;
	}
	else cout << " no valid residue type " <<  _resType << endl;
	return -666.66;
}

string amberElec::getItsAtomName(const UInt _resType, const UInt _atomType) const
{
	if (_resType >= 0 && _resType < resNames.size())
	{
		if (_atomType >=0 && _atomType < atomNames[_resType].size())
		{
			return atomNames[_resType][_atomType];
		}
		else cout << " no atom name for atom " << _atomType << " in res " << resNames[_resType] << endl;
	}
	else cout << " no valid residue type " <<  _resType << endl;
	return "UNK";
}

void amberElec::buildElectrostatics()
{
	itsFileName = "amber.prep";
	buildDataBase();
	//cout << " AMBER all atom electrostatics force field built successfully\n";
	return;
}


void amberElec::buildDataBase()
{
	string evname = "PROTCADDIR";
	string path = getEnvironmentVariable(evname);

	path += "/data/";
	string iFile = path + itsFileName;
	ifstream inFile;
	string currentLine;
	StrVec parsedStrings;
	parsedStrings.resize(0);

	inFile.open(iFile.c_str());
	if (!inFile)
	{
		cout << "Error: unable to open input file: " ;
		cout << iFile << endl;
		exit (1);
	}
	while (getline (inFile, currentLine))
	{
		parsedStrings=Parse::parse(currentLine);
		convertToDataElements(parsedStrings);
		parsedStrings.resize(0);
	}
	//cout << "file read in\n";
	inFile.close();
	inFile.clear();
	orderDataElements();
	return;
}

void amberElec::convertToDataElements(const StrVec& _parsedStrings)
{
	double tmpDbl;
	if (_parsedStrings.size() > 1){
		if (_parsedStrings[1].size() == 3){
			if (_parsedStrings[1][0] == 'I' && _parsedStrings[1][1] == 'N' && _parsedStrings[1][2] == 'T') // #Xaa is the mark of a new residue in the list
			{
				string tmpStr;
				string tmpChar;
				tmpChar.resize(1);
				tmpStr.resize(0);
				for (UInt i = 0; i < _parsedStrings[0].size(); i++) // loop starts after the # mark
				{
					tmpChar = _parsedStrings[0][i];
					tmpStr.append(tmpChar);
				}
				resNames.push_back(tmpStr);
			}
		}
	}
	if (_parsedStrings.size() == 11)
	{
		if (_parsedStrings[2][0] != 'D'){
			UInt i = _parsedStrings.size() - 1; // take last column for the charge
			tmpDbl = 0.0;
			sscanf(_parsedStrings[i].c_str(), "%lf", &tmpDbl);
			if (resNames.size() != atomNames.size()) // if there are fewer atom vectors than their are residues
			{
				vector <string> atomList;
				vector <double> chargeList;
				atomList.push_back(_parsedStrings[1]);
				atomNames.push_back(atomList);
				chargeList.push_back(tmpDbl);
				charges.push_back(chargeList);
			}
			else
			{
				UInt size = atomNames.size();
				atomNames[size - 1].push_back(_parsedStrings[1]);
				charges[size - 1].push_back(tmpDbl);
			}
		}
	}
	return;
}

void amberElec::orderDataElements()  // NOTE:  this requires that the residue type base has already been built
{
	if (residue::dataBase.size() == 0)
	{
		cout << "residue templates not built, cannot order electrostatics elements." << endl;
		return;
	}

	vector <double> tmpDblVec;
	vector <string> tmpStrVec;
	tmpDblVec.resize(0);
	tmpStrVec.resize(0);

	// order the residues


	for (UInt i = 0; i < residue::dataBase.size(); i++)
	{
		bool resFlag = false;
		for (UInt j = i; j < charges.size(); j ++)
		{
			if (resNames[j] == residue::getDataBaseItem(i))
			{
				resFlag = true;
				resNames[j] = resNames[i];
				resNames[i] = residue::getDataBaseItem(i);

				tmpDblVec = charges[j];  	charges[j] = charges[i];  		charges[i] = tmpDblVec;
				tmpStrVec = atomNames[j];  	atomNames[j] = atomNames[i];  	atomNames[i] = tmpStrVec;

				// sort atoms and charges within resiude i

				tmpDblVec = charges[i];
				tmpStrVec = atomNames[i];

				charges[i].resize(0);
				atomNames[i].resize(0);

				for (UInt k = 0; k < residue::getAtomNameBaseSize(i); k++)
				{
					bool atomFlag = false;
					for (UInt l = 0; l < tmpStrVec.size(); l++)
					{
						if (tmpStrVec[l] == residue::getAtomNameBaseItem(i,k)) // if atom name in database matches charge list
						{
							atomNames[i].push_back(tmpStrVec[l]);
							charges[i].push_back(tmpDblVec[l]);
							atomFlag = true;
						}
					}
					if (atomFlag == false) // if atom name not found in charge list, add atom with charge of zero
					{
						atomNames[i].push_back(residue::getAtomNameBaseItem(i,k));
						charges[i].push_back(0.00);
						cout << residue::getDataBaseItem(i) << " no charge for atom " << residue::getAtomNameBaseItem(i,k) << endl;
					}
				}
				if (residue::getAtomNameBaseSize(i) != atomNames[i].size())
					cout << resNames[i] << " charges not properly built!" << endl;
			}
		}
		if (resFlag != true)
			cout << residue::getDataBaseItem(i) << " not found in the charge database." << endl;
	}
	if (resNames.size() != residue::dataBase.size())
		//cout << "Not all residues found during electrostatics database building." << endl;

	return;
}

void amberElec::dummy()
{
	return;
}
