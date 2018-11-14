#include "aaBaseline.h"

double aaBaseline::itsScaleFactor = 1.0;


aaBaseline::aaBaseline()
{	// default constructor
	// read the file until the end.....
	itsFileName = "resBaseline.frc";
	baselineData.resize(0);
	residueNameStrings.resize(0);
	buildDataBase();
	//cout << " aaBaseline database is built " << endl;
/*
	for (UInt i=0; i< R_ref.size(); i++)
	{	cout << amberAtomTypeNames[i] << "   " << R_ref[i];
		cout << "    " << EPS[i] << endl;
	}
*/
}

aaBaseline::aaBaseline(int _dummy)
{	// default constructor
	// read the file until the end.....
	itsFileName = "resBaseline.frc";
	baselineData.resize(0);
	residueNameStrings.resize(0);
	buildDataBase();
	//cout << " aaBaseline database is built " << endl;
/*
	for (UInt i=0; i< R_ref.size(); i++)
	{	cout << amberAtomTypeNames[i] << "   " << R_ref[i];
		cout << "    " << EPS[i] << endl;
	}
*/
}

aaBaseline::aaBaseline(const aaBaseline& _other)
{
	baselineData = _other.baselineData;
	residueNameStrings = _other.residueNameStrings;
	itsFileName = _other.itsFileName;
}

aaBaseline::~aaBaseline()
{
}

double aaBaseline::getEnergy(const string& _name) const
{
	double energy = 0.0;
	for (UInt i=0; i< baselineData.size(); i++)
	{
		if (_name == residueNameStrings[i])
		{
			energy = baselineData[i];
		}
		energy *= itsScaleFactor;
		return energy;
	}
	return energy;
}

vector <string> aaBaseline::list() const
{
    return residueNameStrings;
}

void aaBaseline::buildDataBase()
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
	{	cout << "Error: unable to open input file: ";
		cout << iFile << endl;
		exit (1);
	}

	while (getline(inFile, currentLine, '\n'))
	{	// ignore the comment line
		// comment line should start with #
		if(currentLine[0] != '#' && currentLine[0] != '@'
			&& currentLine[0] != '!' && currentLine[0] != '>'
			&& currentLine[0] != '\n')
		{	parsedStrings=Parse::parse(currentLine);
			convertToDataElements(parsedStrings);
			parsedStrings.resize(0);
		}
	}
	inFile.close();
	inFile.clear();
}

// where specific information about parsed data is intepreted
void aaBaseline::convertToDataElements(const StrVec& _parsedStrings)
{
	double tmpDouble;
	sscanf(_parsedStrings[1].c_str(), "%lf", &tmpDouble);
	baselineData.push_back(tmpDouble);
	residueNameStrings.push_back(_parsedStrings[0]);
}
