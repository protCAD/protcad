#include "referenceSolvation.h"

referenceSolvation::referenceSolvation()
{
	dataFileName = "referenceResidueSolvation.frc";
	buildReferenceSolvationTable();
	printTable();
}

referenceSolvation::~referenceSolvation()
{
}

void referenceSolvation::buildReferenceSolvationTable();
{
	string path = getEnvironmentVariable("PROTCADDIR");
	path = path + "/data/";
	string iFile = path + dataFileName;
	ifstream inFile;
	string currentLine;	
	
	StrVec parsedStrings;
	parsedStrings.resize(0);
	
	inFile.open(iFile.c_str());
	if (!inFile)
	{
		cout << "ERROR:  unable to open referenceSolvation data file!" << endl;
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

	//cout << " reference solvation energy file read in" << endl;
	inFile.close();
	return;
}

void referenceSolvation::convertToDataElements(StrVec& _parsedStrings)
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


