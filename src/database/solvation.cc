#include "solvation.h"

double solvation::itsScaleFactor = 1.0;

solvation::solvation()
{
	dataFileName = "sixAtomSolvation.frc";
	itsParams.resize(0);
	buildSolvationDataBase();
	printParamTable();
}

solvation::~solvation()
{
}

void solvation::buildSolvationDataBase()
{
	string evname = "PROTCADDIR";
	string path = getEnvironmentVariable(evname);
	//cout << path << endl;
	path = path + "/data/";
	//cout << path << " " << dataFileName << endl;
	string iFile = path + dataFileName;
	//cout << iFile << endl;
	ifstream inFile;
	string currentLine;
	StrVec parsedStrings;
	parsedStrings.resize(0);

	inFile.open(iFile.c_str());
	if (!inFile)
	{
		cout << "ERROR:  unable to open solvation data file!" << endl;
		cout << iFile << endl;
		exit(1);
	}

	while (getline (inFile, currentLine, '\n'))
	{
		if (currentLine[0] != '#') // skip comment lines
		{
			parsedStrings = Parse::parse(currentLine);
			convertToDataElements(parsedStrings);
		}
	}

//	cout << " solvation file read in" << endl;
	inFile.close();
	return;
}

void solvation::convertToDataElements(StrVec& _parsedStrings)
{
	vector <double> tempDoubleVec;
	tempDoubleVec.resize(0);

	atomTypeList.push_back(_parsedStrings[0]);

	for (UInt i = 1; i < _parsedStrings.size(); i++)
	{
		double tmpDbl;
		sscanf(_parsedStrings[i].c_str(), "%lf", &tmpDbl);
		tempDoubleVec.push_back(tmpDbl);
	}

	itsParams.push_back(tempDoubleVec);
	return;
}


void solvation::printParamTable()
{
	//cout << "SOLVATION PARAMETER TABLE ..." << endl;
	for (UInt i = 0; i < atomTypeList.size(); i++)
	{		
	//	cout << atomTypeList[i] << " ";
		for (UInt j = 0; j < itsParams[i].size(); j++)
		{
			//cout << itsParams[i][j] << " ";
		}
	//	cout << endl;
	}
	return;
}

int solvation::getIndexFromNameString(string _atomName)
{
	for (int i = 0; i < (int)atomTypeList.size(); i++)
	{
		if (atomTypeList[i] == _atomName)
		{
			return i;
		}
	}
	cout << "atom type " << _atomName << " not found." << endl;
	return -1;
}

double solvation::getSolvationEnergy(double _surfaceArea, UInt _atomType, UInt _paramSet)
{
	if (_atomType >= 0 && _atomType < atomTypeList.size())
	{
		if (_paramSet >=0 && _paramSet < itsParams[_atomType].size())
		{
			return itsScaleFactor*_surfaceArea*itsParams[_atomType][_paramSet] / 1000;
		}
	}
	else
	{
		cout << "ERROR in getSolvationEnergy -- params out of range!" << endl;
		exit(1);
	}
	return 0;
}



