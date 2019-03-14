#include <iostream>
#include <string>
#include <vector>
#include "ensemble.h"
#include "PDBInterface.h"
#include "ran.h"
#include "generalio.h"

vector <string> parse(string& _currentLine);

int main (int argc, char* argv[])
{

	string inputFileName = argv[1];
	ifstream inFile;
	inFile.open(inputFileName.c_str());
	if (!inFile)
	{
		cout << "Unable to find or open file" << endl;
		exit(1);
	}

	string currentLine;
	vector <string> parsedStrings;
	parsedStrings.resize(0);

	vector <atom*> hemeAtoms;

	while (getline(inFile, currentLine, '\n'))
	{
		parsedStrings = parse(currentLine);
		if (parsedStrings[0] == "ATOM")
		{
			atom* newAtom = new atom;
			double xcoord,ycoord,zcoord;
			sscanf(parsedStrings[5].c_str(), "%lf", &xcoord);
			sscanf(parsedStrings[6].c_str(), "%lf", &ycoord);
			sscanf(parsedStrings[7].c_str(), "%lf", &zcoord);
			newAtom->setCoords(xcoord,ycoord,zcoord);
			hemeAtoms.push_back(newAtom);
		}
	}
	cout << "Number of atoms is " << hemeAtoms.size() << endl;
	for (UInt i = 0; i < hemeAtoms.size(); i ++)
	{
		cout << hemeAtoms[i]->getX() << " " << hemeAtoms[i]->getY() << " " << hemeAtoms[i]->getZ() << endl;
	}

	return 0;
}

vector <string> parse(string& _currentLine)
{
	StrVec parsedStrings;
	string tmpStrChi;
	string tmpStr;
	parsedStrings.resize(0);
	tmpStrChi.resize(1);
	tmpStr.resize(0);

	for (UInt i = 0; i < _currentLine.size(); i++)
	{
		if (_currentLine[i] != ' ' && _currentLine[i] != '\t' && _currentLine[i] != '\n' && i != (_currentLine.size()-1))
		{
			tmpStrChi[0] = _currentLine[i];
			tmpStr.append(tmpStrChi);
		}
		else if (tmpStr.size() != 0)
		{
			parsedStrings.push_back(tmpStr);
			tmpStr.resize(0);
		}
	}

	for (UInt i = 0; i < parsedStrings.size(); i++)
	{
		cout << parsedStrings[i] << " ";
	}
	cout << endl;
	return parsedStrings;
}
