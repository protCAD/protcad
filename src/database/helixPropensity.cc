#include "helixPropensity.h"

helixPropensity::helixPropensity()
{
// default constructor
	itsFileName = "helixPropensity.frc";
	itsEnergyMap.resize(0);
	itsHelixPropensityScaleFactor = 1.0;
}

helixPropensity::~helixPropensity()
{
}

double helixPropensity::getEnergy(const UInt _resType)
{
	if (_resType >=0 && _resType < itsEnergyMap.size())
	{
		return itsHelixPropensityScaleFactor*itsEnergyMap[_resType];
	}
	else
	{
		cout << "ERROR in helixPropensity::getEnergy(..)" << endl;
		cout << "\t_resType " << _resType << " passed into function is out of range." << endl;
	}
	return 0;
}

void helixPropensity::buildDatabase()
{
	//cout << "Building Helix Propensity Database" << endl;

	string evname = "PROTCADDIR";
	string path = getEnvironmentVariable(evname);
	path += "/data/";
	string iFile = path + itsFileName;
	ifstream inFile;

	inFile.open(iFile.c_str());
	if (!inFile)
	{   
		cout << "Error: unable to open input file: ";
		cout << iFile << endl;
		exit (1);
	}
	string currentLine;
	while (getline (inFile, currentLine, '\n'))
	{  
		// ignore the comment lines
		// comment line should start with #
		if(currentLine[0] != '#')
        {
			StrVec parsedStrings;  
			parsedStrings = Parse::parse(currentLine);
			addToDatabase(parsedStrings);
        }
    }
    inFile.close();
    inFile.clear();
	return;
}

void helixPropensity::addToDatabase(const StrVec _parsedStrings)
{
	double energy;
	string residueType;
	residueType = _parsedStrings[0];
	sscanf(_parsedStrings[1].c_str(), "%lf", &energy);
	// if this is first addition, initialize energy map
	if (itsEnergyMap.size() == 0)
	{
		UInt size = residue::getDataBaseSize();

		for (UInt i = 0; i < size; i++)
		{
			itsEnergyMap.push_back(0.0);
		}
	}

	for (UInt i = 0; i < itsEnergyMap.size(); i++)
	{
		string tempString = residue::getDataBaseItem(i);
		if (tempString == residueType)
		{
			itsEnergyMap[i] = energy;
			//cout << residueType << " " << energy << endl;	
		} 
	}
	return;
}

void helixPropensity::printEnergyMap()
{
	for (UInt i = 0; i < itsEnergyMap.size(); i ++)
	{
		cout << residue::getDataBaseItem(i) << " " << itsEnergyMap[i] << endl;
	}
	return;
}



