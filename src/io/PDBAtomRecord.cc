#include "PDBAtomRecord.h"
#include <stdio.h>

PDBAtomRecord::PDBAtomRecord()
{
}

PDBAtomRecord::PDBAtomRecord(string& _pdbAtomLine)
{
//	cout << "Entering PDBAtomRecord(string& _pdbAtomLine)" << endl;
	atomName = "";
	resName = "";
	iCode = "";
	altLoc = "";
	resSeq = 0;
	serial = 0;
	occupancy = 0.0;
	tempFactor = 0.0;
	charge = "";
	segID = "";
	element = "";
	chainID = "";
	atomCoord.newsize(3);
	atomCoord[0] = 0.0;
	atomCoord[1] = 0.0;
	atomCoord[2] = 0.0;
	hetflag=false;
	convert(_pdbAtomLine);
}

PDBAtomRecord::PDBAtomRecord(const PDBAtomRecord& _otherRecord)
{
	atomName = _otherRecord.atomName;
	resName = _otherRecord.resName;
	iCode = _otherRecord.iCode;
	altLoc = _otherRecord.altLoc;
	resSeq = _otherRecord.resSeq;
	serial = _otherRecord.serial;
	occupancy = _otherRecord.occupancy;
	tempFactor = _otherRecord.tempFactor;
	charge = _otherRecord.charge;
	chainID = _otherRecord.chainID;
	segID = _otherRecord.segID;
	element = _otherRecord.element;
	atomCoord = _otherRecord.atomCoord;
	hetflag= _otherRecord.hetflag;
}

PDBAtomRecord::~PDBAtomRecord()
{
}

void PDBAtomRecord::convert(string& _pdbAtomLine)
{
	//cout << "Entering PDBAtomRecord::convert(string& _pdbAtomLine)" << endl;
	UInt lengthOfLine = _pdbAtomLine.length();
	//cout << "length of line: " << lengthOfLine << endl;
	int tempint;
	double tempDouble;
	//cout << _pdbAtomLine << endl;

// Get atom serial number
	string lineType= _pdbAtomLine.substr(0,6);
	if(lineType=="HETATM"){hetflag=true;}
	
	string sSerial = _pdbAtomLine.substr(6,5);
	//cout << sSerial << endl;
	sscanf(sSerial.c_str(), "%d", &tempint);
	serial = tempint;
	//cout << "Serial Number= " << serial << endl;

	string atomNameString = _pdbAtomLine.substr(12,4);
	//cout << atomNameString << ":";
	// We only want the alphanumeric portion of the atom
	// name, stripping away any leading or trailing spaces
	UInt firstAlpha = 0;
	UInt lastAlpha = 3;
	for (UInt i=0; i<4; i++)
	{ 	
		//cout << i << " ";
		if (atomNameString.substr(i,1) != " ")
		{	lastAlpha = i;
		}
	}

	for (int i=3; i>=0; i--)
	{	
		//cout << i << " ";
	 	if (atomNameString.substr(i,1) != " ")
		{	firstAlpha = i;
		}
	}
	//cout << "AtomNameString:" << atomNameString << ":" << firstAlpha;
	//cout << ":" << lastAlpha << endl;

	ASSERT(firstAlpha+((lastAlpha - firstAlpha) + 1) <=4);
	atomName = atomNameString.substr(firstAlpha,(lastAlpha-firstAlpha)+1);
    /*if (atomName == "FE1")
    {
        atomName = "F1";
    }
    if (atomName == "FE2")
    {
        atomName = "F2";
    }
    if (atomName == "FE3")
    {
        atomName = "F3";
    }
    if (atomName == "FE4")
    {
        atomName = "F4";
    }*/

	//cout << "atomName= " <<atomName << endl;

	altLoc = _pdbAtomLine.substr(16,1);
	//cout << "altLoc= " <<altLoc << endl;

        
	resName= _pdbAtomLine.substr(17,3);
    if (resName == "HIS")  // default to primary protonation state
    {
        resName = "HIE";
    }
	if (resName == "WAT")  // default to pdb water nomenclature
    {
        resName = "HOH";
    }

	chainID = _pdbAtomLine.substr(21,1);

	string sResSeq = _pdbAtomLine.substr(22,4);
	sscanf(sResSeq.c_str(), "%d", &tempint);
	resSeq = tempint;
	//cout << "resSeq= " <<resSeq << endl;

	iCode = _pdbAtomLine.substr(26,1);
	//cout << "iCode= " <<iCode << endl;

	string sX = _pdbAtomLine.substr(30,8);
	//cout << sX << endl;
	sscanf( sX.c_str(), "%lf", &tempDouble);
	atomCoord[0] = tempDouble;
	//cout << atomCoord[0] << endl;

	string sY = _pdbAtomLine.substr(38,8);
	//cout << sY << endl;
	sscanf( sY.c_str(), "%lf", &tempDouble);
	atomCoord[1] = tempDouble;
	//cout << atomCoord[1] << endl;

	string sZ = _pdbAtomLine.substr(46,8);
	//cout << sZ << endl;
	sscanf( sZ.c_str(), "%lf", &tempDouble);
	atomCoord[2] = tempDouble;
	//cout << atomCoord[2] << endl;

        if (lengthOfLine > 58)
	{
		string sOccupancy = _pdbAtomLine.substr(54,6);
		//cout << sOccupancy << " : ";
		sscanf( sOccupancy.c_str(), "%lf", &tempDouble);
		occupancy = tempDouble;
		//cout << occupancy << endl;
	}
	else
	{
		occupancy = 0.00;
	}

        if (lengthOfLine > 64)
	{
		string sTempFactor = _pdbAtomLine.substr(60,6);
		//cout << sTempFactor << " : ";
		sscanf( sTempFactor.c_str(), "%lf", &tempDouble);
		tempFactor = tempDouble;
		//cout << tempFactor << endl;
	}
	else
	{
		tempFactor = 0.00;
	}
	
        if (lengthOfLine > 75)
	{
		segID = _pdbAtomLine.substr(72,4);
		//cout << "segID = " << segID << endl;
	}
	else
	{
		segID = "    ";
	}

        if (lengthOfLine > 76)
	{
		string elementString = _pdbAtomLine.substr(76,2);
                //cout << "elementString=*" << elementString << "*\n";
                element=elementString;
		if (elementString.substr(0,1) == " ")
		{	element = elementString.substr(1,1);
		}
		if (elementString.substr(1,1) == " ")
		{	element = elementString.substr(0,1);
		}

		//cout << "element = " << element << endl;
	}
	else
	{
		element = "";
	}

	if (lengthOfLine > 79)
	{
		charge = _pdbAtomLine.substr(78,2);
		//cout << "charge = " << charge << endl;
	}
	else
	{	
		charge = "  ";
	}

}
