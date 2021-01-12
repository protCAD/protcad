#include "PDBInterface.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <stdio.h>

aaBaseline PDBInterface::itsAABaseline(0);
string chainIdString[] = {"A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V","W","X","Y","Z","a","b","c","d","e","f","g","h","i","j","k","l","m","n","o","p","q","r","s","t","u","v","w","x","y","z","1","2","3","4",
"5","6","7","8","9"};
UInt chainIDSize = sizeof(chainIdString)/sizeof(chainIdString[0]);
PDBInterface::PDBInterface()
{
}

PDBInterface::PDBInterface(const string& _filename)
{
	resNames = itsAABaseline.list();
	ifstream inFile(_filename.c_str());
	if(inFile.fail())
	{	cout << "ERROR: cannot open file: " << _filename << endl;
	}
	pItsEnsemble = 0;
	readData(inFile);
	inFile.close();
	inFile.clear();
	categorizeLines();
	parseAtomLine();
}

// _hetflag=true if you want to read in hetatms from file.
// _atomflag=true if you want to read in protein atoms from file.
PDBInterface::PDBInterface(const string& _filename, bool _hetflag, bool _atomflag)
{
	ifstream inFile(_filename.c_str());
	if(!inFile){
		cout << "Error: cannot open file: " << _filename << endl;
	}
	//cout << "\nStarting pItsEnsemble";
	pItsEnsemble=0;
	//cout << "\nStarting readData(inFile)";
	readData(inFile);
	//cout << "\nStarting categorizeLines()";
	categorizeLines();

	setHetatmFlag(_hetflag);
	setAtomFlag(_atomflag);
	
	if(itsAtomFlag){parseAtomLine();}
	if(itsHetatmFlag){parseHetatmLine();}
}

PDBInterface::PDBInterface(const string& _filename, const bool _Hflag)
{	ifstream inFile(_filename.c_str());
	if(!inFile)
	{	cout << "ERROR: cannot open file: " << _filename << endl;
	}
	pItsEnsemble = 0;
	readData(inFile);
	inFile.close();
	inFile.clear();
	categorizeLines();
	parseAtomLine(_Hflag);
}

PDBInterface::PDBInterface(const string& _filename, const bool _Hflag, const bool _HPflag, const bool _dummyBool)
{	ifstream inFile(_filename.c_str());
	if(!inFile)
	{	cout << "ERROR: cannot open file : " << _filename << endl;
	}
	pItsEnsemble = 0;
	readData(inFile);
	inFile.close();
	inFile.clear();
	categorizeLines();
	if((_Hflag) && (_HPflag))
	{
		cout << "_Hflag and _HPflag can't both be true in PDBInterface.\n";
		cout << "In PDBInterface, setting only _HPflag as true.\n";
		parseAtomLine(false,true);
	}
	else
	{	parseAtomLine(_Hflag,_HPflag); }
}


PDBInterface::PDBInterface(const PDBInterface& _other)
{
}

PDBInterface::~PDBInterface()
{
	delete pItsEnsemble;
}

void PDBInterface::readData(ifstream& _infile)
{
	string linebuffer;
	bool inLibrary;
	while (getline(_infile,linebuffer,'\n'))
	{
		if (linebuffer.length() > 16)
		{
			inLibrary = false;
			for (UInt i = 0; i < resNames.size(); i++)
			{
				if (linebuffer.compare(17,3,resNames[i]) == 0)
				{
					inLibrary = true;
					break;
				}
			}
			if (inLibrary)
			{
				if (linebuffer.compare(13,1,"H") != 0 || linebuffer.compare(12,1,"H") != 0 || linebuffer.compare(16,1," ") != 0)// do not include hydrogens as we will use internal nomenclature and inclusion criteria
				{
					theLines.push_back(linebuffer);
				}
			}
		}
	}
}

void PDBInterface::categorizeLines()
{
	string header;
	for (UInt i=0; i<theLines.size(); i++)
	{
		header = theLines[i].substr(0,6);
		if (header == "ATOM  " || header == "HETATM")
		{	atomLines.push_back(i);
			continue;
		}
	}
}

void PDBInterface::parseHetatmLine()
{
	//cout << "\nIn parseHetatmLine..." << hetatmLines.size() << " <--- number of hetatm lines\n";
		
	string tempstring,chainIDtemp="EMPTY";
        
	if (!pItsEnsemble){pItsEnsemble = new ensemble();}
        //cout << "Done making a new ensemble"<<endl;
}

void PDBInterface::parseAtomLine()
{	
	//cout << "Entering PDBInterface::parseAtomLine()" << endl;
	// First, Let's find out how many chains we have in the file
	// We're going to assume (for now) that there's at least one...

	UInt numchains = 1;
	vector<string> chainsSeen;
	string theProteinName;
	int counter = 1;
	bool iveSeenIt;
	string lastSeenChain;
	string tempstring;
	//UInt nextChainID;
	for (UInt i=0; i< atomLines.size(); i++)
	{
		tempstring = theLines[atomLines[i]];
		PDBAtomRecord currentRecord(tempstring);
		//cout << currentRecord.getChainID() << " ";
		if (counter == 1)
		{	// This is the first residue
			chainsSeen.push_back(currentRecord.getChainID());
			theProteinName = currentRecord.getSegID();
			counter++;
			continue;
		}
		iveSeenIt = false;
		for (UInt i=0; i<chainsSeen.size(); i++)
		{	if (chainsSeen[i] == currentRecord.getChainID())
			{	iveSeenIt = true;
			}
		}
		if (!iveSeenIt)
		{	
			numchains++;
			chainsSeen.push_back(currentRecord.getChainID());
			//if (currentRecord.getChainID() != " "){lastSeenChain = currentRecord.getChainID();}
		}
		counter++;
	}
	// find next chain to seperate cofactors by
	/*for (UInt i=0; i < chainIDSize; i++)
	{
		if (lastSeenChain == chainIdString[i]){nextChainID = i+1;}
	}*/

	//cout << endl;
	//cout << "Number of unique protein chains in file: " << numchains << endl;
	// OK, now build the protein, and add the chains in (empty at first)
	protein* pTheProtein = new protein(theProteinName);
	vector<chain*> vecChainPointers;
	for (UInt i=0; i<numchains; i++)
	{
		const char* pTheChainID  = chainsSeen[i].c_str();
		char theChainID = *pTheChainID;
		chain* theChain = new chain(theChainID);
		pTheProtein->add(theChain);
		vecChainPointers.push_back(theChain);
	}
	// Now, start looking for residues - we'll assume that the pdb
	// file is numbered "sanely" and that different residues will 
	// be distinguished by different residue numbers and/or insertion
	// codes
	counter = 1;
	string lastResName;
	string currentResName;
    int lastResSeq = 0;
	int currentResSeq;
	string lastICode;
	string currentICode;
	string currentChainID;
	string lastChainID;
	string currentAltLoc;
	string lastAltLoc;
	vector<UInt> resbegin;
	vector<UInt> resend;
	vector<string> resname;
	vector<string> resChainID;
	vector<UInt> resnums;
	vector<string> icodes;

	// The following function is called in order to initialize all 
	// the static stuff in the
	// residue class definition, such as the residue database, etc.
	//cout << "presetup" << endl;
	residue::setupDataBase();
	//cout << "setup" << endl;
    bool hydrogensFound = true;
	bool altLocFound = false;
	vector<bool> Hflags;
	vector<bool> altLocFlags;

	for (UInt i=0; i<atomLines.size(); i++)
	{	
		tempstring = theLines[atomLines[i]];
		PDBAtomRecord currentRecord(tempstring);
		if (counter == 1)
		{	// this is the first residue
			currentResName = currentRecord.getResName();
			lastResName = currentResName;
			currentResSeq = currentRecord.getResSeq();
			lastResSeq = currentResSeq;
			currentICode = currentRecord.getICode();
			lastICode = currentICode;
			currentChainID = currentRecord.getChainID();
			lastChainID = currentChainID;
			currentAltLoc = currentRecord.getAltLoc();
			lastAltLoc = currentAltLoc;
			resbegin.push_back(i);
			resChainID.push_back(currentChainID);
			resname.push_back(currentResName);
			resnums.push_back(currentResSeq);
			icodes.push_back(currentICode);
			counter++;
			continue;
		}
		currentResName = currentRecord.getResName();
		currentResSeq = currentRecord.getResSeq();
		currentICode = currentRecord.getICode();
		currentChainID = currentRecord.getChainID();
		currentAltLoc = currentRecord.getAltLoc();
		if (currentResSeq != lastResSeq)
		{	// this is a new residue
			//cout << "start of new res at " << i << endl;
			Hflags.push_back(hydrogensFound);
			altLocFlags.push_back(altLocFound);
            hydrogensFound = true;
			altLocFound = false;
			resend.push_back(i-1);
			resbegin.push_back(i);
			resname.push_back(currentResName);
			resnums.push_back(currentResSeq);
			resChainID.push_back(currentChainID);
			icodes.push_back(currentICode);
			lastResSeq = currentResSeq;
			lastICode = currentICode;
			lastChainID = currentChainID;

		}
		if (currentResSeq == lastResSeq && currentICode != lastICode)
		{	// this is a new residue
			//cout << "start of new res at " << i << endl;
			Hflags.push_back(hydrogensFound);
			altLocFlags.push_back(altLocFound);
            hydrogensFound = true;
			altLocFound = false;
			resend.push_back(i-1);
			resname.push_back(currentResName);
			resnums.push_back(currentResSeq);
			resbegin.push_back(i);
			lastResSeq = currentResSeq;
			lastICode = currentICode;
		}
		counter++;
		if (currentRecord.getElement() == "" || currentRecord.getElement() == " " || currentRecord.getElement() == "  " )
		{ 
			// more robust method in case the element field has
			// been left blank
			// Figure out what kind of residue we're dealing with	
			UInt theType = 9999;
			for (UInt j=0; j<residue::getDataBaseSize(); j++)
			{	if (currentRecord.getResName() == residue::getDataBaseItem(j))
				{	theType = j;
					break;
				}
			}
			// OK, now we now the residue type, let's find the
			// atom name amongst the atoms in the template
			int theAtomTypeIndex = residue::dataBase[theType].getAtomIndexOf(currentRecord.getAtomName());
			if (theAtomTypeIndex < 0)
			{
				if (currentRecord.getAtomName() != "H1" && currentRecord.getAtomName() != "H2" && currentRecord.getAtomName() != "H3" && currentRecord.getAtomName() != "OXT")
				{
					cout << "Don't understand atom name ";
					cout << currentRecord.getResName() << " " << currentRecord.getAtomName() << endl;
					cout << "Further behavior unpredictable" << endl;
				}
			}
		}
		if (currentAltLoc != lastAltLoc)
		{	altLocFound = true;
		}

	}
	resend.push_back(atomLines.size()-1);
    Hflags.push_back(hydrogensFound);
	altLocFlags.push_back(altLocFound);

	counter = 0;
	for (UInt i=0; i<atomLines.size(); i++)
	{	
		tempstring = theLines[atomLines[i]];
		PDBAtomRecord currentRecord(tempstring);
	}

	// now we're going to loop over the residues and
	// generate them, with an internal loop to add the
	// atoms
	for (UInt i=0; i< resbegin.size(); i++)
	{	
		// Figure out what kind of residue we're dealing with	
		UInt theType = 9999;
		for (UInt j=0; j<residue::getDataBaseSize(); j++)
		{	if (resname[i] == residue::getDataBaseItem(j))
			{	theType = j;
				break;
			}
		}
        residue* pTheResidue = new residue(theType,Hflags[i]);
		pTheResidue->setResNum(resnums[i]);

		// if residue is a cofactor add it to it's own chain
		/*if (pTheResidue->isCofactor() && nextChainID < chainIDSize){
			const char* pTheChainID  = chainIdString[nextChainID].c_str();
			resChainID[i] = chainIdString[nextChainID];
			char theChainID = *pTheChainID;
			chain* theChain = new chain(theChainID);
			pTheProtein->add(theChain);
			vecChainPointers.push_back(theChain);
			nextChainID++;
		}*/

		// OK, now which chain do we add this to?
		chain* pCurrentChain = 0;
		char theChainID;
		for (UInt j=0; j<vecChainPointers.size(); j++)
		{
			const char* pTheChainID  = resChainID[i].c_str();
			theChainID = *pTheChainID;
			if ( theChainID == (vecChainPointers[j])->getChainID())
			{	pCurrentChain = vecChainPointers[j];
				break;
			}
		}
		// if we've found the right chain, add the residue,
		// otherwise, fail hard
		if (pCurrentChain)
		{
			pCurrentChain->add(pTheResidue);
		}
		else
		{
			cout << "Wasn't able to find chain named: ";
			cout << "unfoundChain:" << theChainID << endl;
			//cout << "Please fix your pdb file!" << endl;
			//cout << "Further behavior unpredictable  - Stopping" << endl;
			//break;
		}
		// Now, add the atoms from the start of the residue to
		// the end of the residue
		// cout << resbegin[i] << " : " << resend[i] << endl;
		UInt numAtomsInRes = 0;
		if (altLocFlags[i])
		{
			string bestAltLoc = "";
			double bestOccupancy = 0.0;
			for (UInt j=resbegin[i]; j<=resend[i]; j++)
			{
				tempstring = theLines[atomLines[j]];
				PDBAtomRecord currentRecord(tempstring);
				if (currentRecord.getAltLoc() != " ")
				{	if (currentRecord.getOccupancy() > bestOccupancy)
					{	bestOccupancy = currentRecord.getOccupancy();
						bestAltLoc = currentRecord.getAltLoc();
					}
				}
			}
			//cout << "Best altLoc = " << bestAltLoc << endl;

			for (UInt j=resbegin[i]; j<=resend[i]; j++)
			{
				tempstring = theLines[atomLines[j]];
				PDBAtomRecord currentRecord(tempstring);
				if (currentRecord.getAltLoc() == " " ||
				    currentRecord.getAltLoc() == bestAltLoc)
				{
					pTheResidue->addAtom(currentRecord);
					numAtomsInRes++;
				}
			}
		}
		else
		{
			for (UInt j=resbegin[i]; j<=resend[i]; j++)
			{	
				tempstring = theLines[atomLines[j]];
				PDBAtomRecord currentRecord(tempstring);
				pTheResidue->addAtom(currentRecord);
				numAtomsInRes++;
			}
		}
		// At this point all the atoms should be initialized in our
		// new residue, and the residue should be in the chain, which
		// is in the protein.
	} // end loop over residues

	parseBIOMT();

	pTheProtein->finishProteinBuild();

	if (!pItsEnsemble)
	{	pItsEnsemble = new ensemble();
	}
	pItsEnsemble->add(pTheProtein);
}

void PDBInterface::parseAtomLine(const bool _Hflag)
{
	//cout << "Entering PDBInterface::parseAtomLine()" << endl;
	// First, Let's find out how many chains we have in the file
	// We're going to assume (for now) that there's at least one...
	UInt numchains = 1;
	vector<string> chainsSeen;
	string theProteinName;
	int counter = 1;
	bool iveSeenIt;
	string tempstring;
	for (UInt i=0; i< atomLines.size(); i++)
	{
		tempstring = theLines[atomLines[i]];
		PDBAtomRecord currentRecord(tempstring);
		//cout << currentRecord.getChainID() << " ";
		if (counter == 1)
		{	// This is the first residue
			chainsSeen.push_back(currentRecord.getChainID());
			theProteinName = currentRecord.getSegID();
			counter++;
			continue;
		}
		iveSeenIt = false;
		for (UInt i=0; i<chainsSeen.size(); i++)
		{	if (chainsSeen[i] == currentRecord.getChainID())
			{	iveSeenIt = true;
			}
		}
		if (!iveSeenIt)
		{
			numchains++;
			chainsSeen.push_back(currentRecord.getChainID());
		}
		counter++;
	}
	//cout << endl;
	//cout << "Number of unique protein chains in file: " << numchains << endl;

	// OK, now build the protein, and add the chains in (empty at first)
	protein* pTheProtein = new protein(theProteinName);
	vector<chain*> vecChainPointers;
	for (UInt i=0; i<numchains; i++)
	{
		const char* pTheChainID  = chainsSeen[i].c_str();
		char theChainID = *pTheChainID;
		chain* theChain = new chain(theChainID);
		pTheProtein->add(theChain);
		vecChainPointers.push_back(theChain);
	}
	// Now, start looking for residues - we'll assume that the pdb
	// file is numbered "sanely" and that different residues will 
	// be distinguished by different residue numbers and/or insertion
	// codes
	counter = 1;
	string lastResName;
	string currentResName;
    int lastResSeq = 0;
	int currentResSeq;
	string lastICode;
	string currentICode;
	string currentChainID;
	string lastChainID;
	string currentAltLoc;
	string lastAltLoc;
	vector<UInt> resbegin;
	vector<UInt> resend;
	vector<string> resname;
	vector<string> resChainID;
	vector<UInt> resnums;
	vector<string> icodes;

	// The following function is called in order to initialize all
	// the static stuff in the
	// residue class definition, such as the residue database, etc.

	residue::setupDataBase(_Hflag);
	bool hydrogensFound = false;
	bool altLocFound = false;
	vector<bool> Hflags;
	vector<bool> altLocFlags;

	for (UInt i=0; i<atomLines.size(); i++)
	{	
		tempstring = theLines[atomLines[i]];
		PDBAtomRecord currentRecord(tempstring);
		if (counter == 1)
		{	// this is the first residue
			currentResName = currentRecord.getResName();
			lastResName = currentResName;
			currentResSeq = currentRecord.getResSeq();
			lastResSeq = currentResSeq;
			currentICode = currentRecord.getICode();
			lastICode = currentICode;
			currentChainID = currentRecord.getChainID();
			lastChainID = currentChainID;
			currentAltLoc = currentRecord.getAltLoc();
			lastAltLoc = currentAltLoc;
			resbegin.push_back(i);
			resChainID.push_back(currentChainID);
			resname.push_back(currentResName);
			resnums.push_back(currentResSeq);
			icodes.push_back(currentICode);
			counter++;
			continue;
		}
		currentResName = currentRecord.getResName();
		currentResSeq = currentRecord.getResSeq();
		currentICode = currentRecord.getICode();
		currentChainID = currentRecord.getChainID();
		currentAltLoc = currentRecord.getAltLoc();
		if (currentResSeq != lastResSeq)
		{	// this is a new residue
			//cout << "start of new res at " << i << endl;
			Hflags.push_back(hydrogensFound);
			altLocFlags.push_back(altLocFound);
			hydrogensFound = false;
			altLocFound = false;
			resend.push_back(i-1);
			resbegin.push_back(i);
			resname.push_back(currentResName);
			resnums.push_back(currentResSeq);
			resChainID.push_back(currentChainID);
			icodes.push_back(currentICode);
			lastResSeq = currentResSeq;
			lastICode = currentICode;
			lastChainID = currentChainID;

		}
		if (currentResSeq == lastResSeq && currentICode != lastICode)
		{	// this is a new residue
			//cout << "start of new res at " << i << endl;
			Hflags.push_back(hydrogensFound);
			altLocFlags.push_back(altLocFound);
			hydrogensFound = false;
			altLocFound = false;
			resend.push_back(i-1);
			resname.push_back(currentResName);
			resnums.push_back(currentResSeq);
			resbegin.push_back(i);
			lastResSeq = currentResSeq;
			lastICode = currentICode;
		}
        counter++;
		// Now, check for the presence of hydrogens in each of the residues,
		// if they exist, set the value of Hflag to true
		if (currentRecord.getElement() == "H")
        {
            //cout << "Found a hydrogen at " << i << endl;
            hydrogensFound = true;
        }
		if (currentRecord.getElement() == "" ||
		    currentRecord.getElement() == " " ||
                    currentRecord.getElement() == "  " )
		{
			// more robust method in case the element field has
			// been left blank
			//cout << "PDBInterface::Null element field! "<< endl;
			//cout << "Checking against residue type database" << endl;
			// Figure out what kind of residue we're dealing with	
			UInt theType = 9999;
			for (UInt j=0; j<residue::getDataBaseSize(); j++)
			{	if (currentRecord.getResName() == residue::getDataBaseItem(j))
				{	theType = j;
					break;
				}
			}
			if (theType == 9999)
			{
				cout << "Don't know how to build a residue of type ";
				cout << currentRecord.getResName() << endl;
				cout << "Please add it to the database!" << endl;
				cout << "Further behavior unpredictable  - Stopping" << endl;
				terminate();
			}
			// OK, now we now the residue type, let's find the
			// atom name amongst the atoms in the template
			//cout << "Querying atom names" << endl;
			int theAtomTypeIndex = residue::dataBase[theType].getAtomIndexOf(currentRecord.getAtomName());
			if (theAtomTypeIndex >= 0)
			{
				//cout << "atomTypeIndex = " << theAtomTypeIndex << endl;
				string atomTypeString = (residue::dataBase[theType].atomList[theAtomTypeIndex]).getType();
				//cout << "atomTypeString = " << atomTypeString << endl;
				if (atomTypeString == "H")
				{ hydrogensFound = true;
				}
			}
			else
			{
				cout << "Don't understand atom name ";
				cout << currentRecord.getAtomName() << endl;
				cout << "Further behavior unpredictable" << endl;
			}
		}
		if (currentAltLoc != lastAltLoc)
		{	altLocFound = true;
		}

	}
	resend.push_back(atomLines.size()-1);
	Hflags.push_back(hydrogensFound);
	altLocFlags.push_back(altLocFound);
/*	
	for (UInt i=0; i<resend.size(); i++)
	{	cout << resname [i] << "  " << resbegin[i] << "  " << resend[i];
		cout << "  " << Hflags[i] << "  " << altLocFlags[i] << endl;
	}
*/
	counter = 0;
	for (UInt i=0; i<atomLines.size(); i++)
	{	
		tempstring = theLines[atomLines[i]];
		PDBAtomRecord currentRecord(tempstring);
		if (counter == 0)
		{	
		}
	}
	// now we're going to loop over the residues and
	// generate them, with an internal loop to add the
	// atoms

	for (UInt i=0; i< resbegin.size(); i++)
	{	
		// Figure out what kind of residue we're dealing with
		UInt theType = 9999;
		for (UInt j=0; j<residue::getDataBaseSize(); j++)
		{	if (resname[i] == residue::getDataBaseItem(j))
			{	theType = j;
				break;
			}
		}
		if (theType == 9999)
		{
			cout << "Don't know how to build a residue of type ";
			cout << resname[i] << endl;
			cout << "Please add it to the database!" << endl;
			cout << "Further behavior unpredictable  - Stopping" << endl;
			terminate();
		}
		//cout << residue::getDataBaseItem(theType) << "  "<< "Hflags[" << i << "] = " << Hflags[i] << endl;
		residue* pTheResidue = new residue(theType,Hflags[i]);
		pTheResidue->setResNum(resnums[i]);
		chain* pCurrentChain = 0;
		// OK, now which chain do we add this to?
		for (UInt j=0; j<vecChainPointers.size(); j++)
		{
			const char* pTheChainID  = resChainID[i].c_str();
            char theChainID = *pTheChainID;
			if ( theChainID == (vecChainPointers[j])->getChainID())
			{	pCurrentChain = vecChainPointers[j];
				break;
			}
		}
		// if we've found the right chain, add the residue,
		// otherwise, fail hard
		if (pCurrentChain)
		{
			pCurrentChain->add(pTheResidue);
		}
		else
		{
            cout << "Wasn't able to find chain";
            //cout << theChainID << endl;
			cout << "Please fix your pdb file!" << endl;
			cout << "Further behavior unpredictable  - Stopping" << endl;
			terminate();	
		}
		// Now, add the atoms from the start of the residue to
		// the end of the residue
		// cout << resbegin[i] << " : " << resend[i] << endl;
		UInt numAtomsInRes = 0;
		if (altLocFlags[i])
		{
			string bestAltLoc = "";
			double bestOccupancy = 0.0;
			for (UInt j=resbegin[i]; j<=resend[i]; j++)
			{
				tempstring = theLines[atomLines[j]];
				PDBAtomRecord currentRecord(tempstring);
				if (currentRecord.getAltLoc() != " ")
				{	if (currentRecord.getOccupancy() > bestOccupancy)
					{	bestOccupancy = currentRecord.getOccupancy();
						bestAltLoc = currentRecord.getAltLoc();
					}
				}
			}
			//cout << "Best altLoc = " << bestAltLoc << endl;

			for (UInt j=resbegin[i]; j<=resend[i]; j++)
			{
				tempstring = theLines[atomLines[j]];
				PDBAtomRecord currentRecord(tempstring);
				if (currentRecord.getAltLoc() == " " ||
				    currentRecord.getAltLoc() == bestAltLoc)
				{
					pTheResidue->addAtom(currentRecord);
					numAtomsInRes++;
				}
			}
		}
		else
		{
			for (UInt j=resbegin[i]; j<=resend[i]; j++)
			{	
				tempstring = theLines[atomLines[j]];
				PDBAtomRecord currentRecord(tempstring);
				pTheResidue->addAtom(currentRecord);
				numAtomsInRes++;
			}
		}
		// At this point all the atoms should be initialized in our
		// new residue, and the residue should be in the chain, which
		// is in the protein.
		if (numAtomsInRes < pTheResidue->getNumAtoms())
		{
			cout << "Too few atoms in residue ";
			cout << pTheResidue->getResNum() ;
			cout << " - expected: " << pTheResidue->getNumAtoms();
			cout << " found: " <<  numAtomsInRes << endl;
			pCurrentChain->fixBrokenResidue(pCurrentChain->getNumResidues()-1);
		}
		else if (numAtomsInRes > pTheResidue->getNumAtoms())
		{
			cout << "Too many atoms in residue ";
			cout << pTheResidue->getResNum() ;
			cout << " - expected: " << pTheResidue->getNumAtoms();
			cout << " found: " << numAtomsInRes << endl;
			//pCurrentChain->fixBrokenResidue(pCurrentChain->getNumResidues()-1);
		}
	} // end loop over residues

	parseBIOMT();

	pTheProtein->finishProteinBuild();

	if (!pItsEnsemble)
	{	pItsEnsemble = new ensemble();
	}
	pItsEnsemble->add(pTheProtein);
}


void PDBInterface::parseAtomLine(const bool _Hflag, const bool _HPflag)
{
	//cout << "Entering PDBInterface::parseAtomLine()" << endl;
	// First, Let's find out how many chains we have in the file
	// We're going to assume (for now) that there's at least one...
	UInt numchains = 1;
	vector<string> chainsSeen;
	string theProteinName;
	int counter = 1;
	bool iveSeenIt;
	string tempstring;
	for (UInt i=0; i< atomLines.size(); i++)
	{
		tempstring = theLines[atomLines[i]];
		PDBAtomRecord currentRecord(tempstring);
		//cout << currentRecord.getChainID() << " ";
		if (counter == 1)
		{	// This is the first residue
			chainsSeen.push_back(currentRecord.getChainID());
			theProteinName = currentRecord.getSegID();
			counter++;
			continue;
		}
		iveSeenIt = false;
		for (UInt i=0; i<chainsSeen.size(); i++)
		{	if (chainsSeen[i] == currentRecord.getChainID())
			{	iveSeenIt = true;
			}
		}
		if (!iveSeenIt)
		{
			numchains++;
			chainsSeen.push_back(currentRecord.getChainID());
		}
		counter++;
	}
	//cout << endl;
	//cout << "Number of unique protein chains in file: " << numchains << endl;

	// OK, now build the protein, and add the chains in (empty at first)
	protein* pTheProtein = new protein(theProteinName);
	vector<chain*> vecChainPointers;
	for (UInt i=0; i<numchains; i++)
	{
		const char* pTheChainID  = chainsSeen[i].c_str();
		char theChainID = *pTheChainID;
		chain* theChain = new chain(theChainID);
		pTheProtein->add(theChain);
		vecChainPointers.push_back(theChain);
	}
	// Now, start looking for residues - we'll assume that the pdb
	// file is numbered "sanely" and that different residues will 
	// be distinguished by different residue numbers and/or insertion
	// codes
	counter = 1;
	string lastResName;
	string currentResName;
    int lastResSeq = 0;
	int currentResSeq;
	string lastICode;
	string currentICode;
	string currentChainID;
	string lastChainID;
	string currentAltLoc;
	string lastAltLoc;
	vector<UInt> resbegin;
	vector<UInt> resend;
	vector<string> resname;
	vector<string> resChainID;
	vector<UInt> resnums;
	vector<string> icodes;

	// The following function is called in order to initialize all
	// the static stuff in the
	// residue class definition, such as the residue database, etc.

	residue::setupDataBase(_Hflag,_HPflag);
	bool hydrogensFound = false;
	bool altLocFound = false;
	vector<bool> Hflags;
	vector<bool> altLocFlags;

	for (UInt i=0; i<atomLines.size(); i++)
	{	
		tempstring = theLines[atomLines[i]];
		PDBAtomRecord currentRecord(tempstring);
		if (counter == 1)
		{	// this is the first residue
			currentResName = currentRecord.getResName();
			lastResName = currentResName;
			currentResSeq = currentRecord.getResSeq();
			lastResSeq = currentResSeq;
			currentICode = currentRecord.getICode();
			lastICode = currentICode;
			currentChainID = currentRecord.getChainID();
			lastChainID = currentChainID;
			currentAltLoc = currentRecord.getAltLoc();
			lastAltLoc = currentAltLoc;
			resbegin.push_back(i);
			resChainID.push_back(currentChainID);
			resname.push_back(currentResName);
			resnums.push_back(currentResSeq);
			icodes.push_back(currentICode);
			counter++;
			continue;
		}
		currentResName = currentRecord.getResName();
		currentResSeq = currentRecord.getResSeq();
		currentICode = currentRecord.getICode();
		currentChainID = currentRecord.getChainID();
		currentAltLoc = currentRecord.getAltLoc();
		if (currentResSeq != lastResSeq)
		{	// this is a new residue
			//cout << "start of new res at " << i << endl;
			Hflags.push_back(hydrogensFound);
			altLocFlags.push_back(altLocFound);
			hydrogensFound = false;
			altLocFound = false;
			resend.push_back(i-1);
			resbegin.push_back(i);
			resname.push_back(currentResName);
			resnums.push_back(currentResSeq);
			resChainID.push_back(currentChainID);
			icodes.push_back(currentICode);
			lastResSeq = currentResSeq;
			lastICode = currentICode;
			lastChainID = currentChainID;

		}
		if (currentResSeq == lastResSeq && currentICode != lastICode)
		{	// this is a new residue
			//cout << "start of new res at " << i << endl;
			Hflags.push_back(hydrogensFound);
			altLocFlags.push_back(altLocFound);
			hydrogensFound = false;
			altLocFound = false;
			resend.push_back(i-1);
			resname.push_back(currentResName);
			resnums.push_back(currentResSeq);
			resbegin.push_back(i);
			lastResSeq = currentResSeq;
			lastICode = currentICode;
		}
		counter++;
		// Now, check for the presence of hydrogens in each of the residues,
		// if they exist, set the value of Hflag to true
		if (currentRecord.getElement() == "H")
		{
			//cout << "Found a hydrogen at " << i << endl;
			hydrogensFound = true;
		}
		if (currentRecord.getElement() == "" ||
		    currentRecord.getElement() == " " ||
                    currentRecord.getElement() == "  " )
		{
			// more robust method in case the element field has
			// been left blank
			//cout << "PDBInterface::Null element field! "<< endl;
			//cout << "Checking against residue type database" << endl;
			// Figure out what kind of residue we're dealing with	
			UInt theType = 9999;
			for (UInt j=0; j<residue::getDataBaseSize(); j++)
			{	if (currentRecord.getResName() == residue::getDataBaseItem(j))
				{	theType = j;
					break;
				}
			}
			if (theType == 9999)
			{
				cout << "Don't know how to build a residue of type ";
				cout << currentRecord.getResName() << endl;
				cout << "Please add it to the database!" << endl;
				cout << "Further behavior unpredictable  - Stopping" << endl;
				terminate();
			}
			// OK, now we now the residue type, let's find the
			// atom name amongst the atoms in the template
			//cout << "Querying atom names" << endl;
			int theAtomTypeIndex = residue::dataBase[theType].getAtomIndexOf(currentRecord.getAtomName());
			if (theAtomTypeIndex >= 0)
			{
				//cout << "atomTypeIndex = " << theAtomTypeIndex << endl;
				string atomTypeString = (residue::dataBase[theType].atomList[theAtomTypeIndex]).getType();
				//cout << "atomTypeString = " << atomTypeString << endl;
				if (atomTypeString == "H")
				{ hydrogensFound = true;
				}
			}
			else
			{
				cout << "Don't understand atom name ";
				cout << currentRecord.getAtomName() << endl;
				cout << "Further behavior unpredictable" << endl;
			}
		}
		if (currentAltLoc != lastAltLoc)
		{	altLocFound = true;
		}

	}
	resend.push_back(atomLines.size()-1);
	Hflags.push_back(hydrogensFound);
	altLocFlags.push_back(altLocFound);
/*	
	for (UInt i=0; i<resend.size(); i++)
	{	cout << resname [i] << "  " << resbegin[i] << "  " << resend[i];
		cout << "  " << Hflags[i] << "  " << altLocFlags[i] << endl;
	}
*/
	counter = 0;
	for (UInt i=0; i<atomLines.size(); i++)
	{	
		tempstring = theLines[atomLines[i]];
		PDBAtomRecord currentRecord(tempstring);
		if (counter == 0)
		{	
		}
	}
	// now we're going to loop over the residues and
	// generate them, with an internal loop to add the
	// atoms

	for (UInt i=0; i< resbegin.size(); i++)
	{	
		// Figure out what kind of residue we're dealing with
		UInt theType = 9999;
		for (UInt j=0; j<residue::getDataBaseSize(); j++)
		{	if (resname[i] == residue::getDataBaseItem(j))
			{	theType = j;
				break;
			}
		}
		if (theType == 9999)
		{
			cout << "Don't know how to build a residue of type ";
			cout << resname[i] << endl;
			cout << "Please add it to the database!" << endl;
			cout << "Further behavior unpredictable  - Stopping" << endl;
			terminate();
		}
		//cout << residue::getDataBaseItem(theType) << "  "<< "Hflags[" << i << "] = " << Hflags[i] << endl;
		residue* pTheResidue = new residue(theType,Hflags[i]);
		pTheResidue->setResNum(resnums[i]);
		chain* pCurrentChain = 0;
		// OK, now which chain do we add this to?
		for (UInt j=0; j<vecChainPointers.size(); j++)
		{
			const char* pTheChainID  = resChainID[i].c_str();
            char theChainID = *pTheChainID;
			if ( theChainID == (vecChainPointers[j])->getChainID())
			{	pCurrentChain = vecChainPointers[j];
				break;
			}
		}
		// if we've found the right chain, add the residue,
		// otherwise, fail hard
		if (pCurrentChain)
		{
			pCurrentChain->add(pTheResidue);
		}
		else
		{
            cout << "Wasn't able to find chain";
            //cout << theChainID << endl;
			cout << "Please fix your pdb file!" << endl;
			cout << "Further behavior unpredictable  - Stopping" << endl;
			terminate();	
		}
		// Now, add the atoms from the start of the residue to
		// the end of the residue
		// cout << resbegin[i] << " : " << resend[i] << endl;
		UInt numAtomsInRes = 0;
		if (altLocFlags[i])
		{
			string bestAltLoc = "";
			double bestOccupancy = 0.0;
			for (UInt j=resbegin[i]; j<=resend[i]; j++)
			{
				tempstring = theLines[atomLines[j]];
				PDBAtomRecord currentRecord(tempstring);
				if (currentRecord.getAltLoc() != " ")
				{	if (currentRecord.getOccupancy() > bestOccupancy)
					{	bestOccupancy = currentRecord.getOccupancy();
						bestAltLoc = currentRecord.getAltLoc();
					}
				}
			}
			//cout << "Best altLoc = " << bestAltLoc << endl;

			for (UInt j=resbegin[i]; j<=resend[i]; j++)
			{
				tempstring = theLines[atomLines[j]];
				PDBAtomRecord currentRecord(tempstring);
				if (currentRecord.getAltLoc() == " " ||
				    currentRecord.getAltLoc() == bestAltLoc)
				{
					pTheResidue->addAtom(currentRecord);
					numAtomsInRes++;
				}
			}
		}
		else
		{
			for (UInt j=resbegin[i]; j<=resend[i]; j++)
			{	
				tempstring = theLines[atomLines[j]];
				PDBAtomRecord currentRecord(tempstring);
				pTheResidue->addAtom(currentRecord);
				numAtomsInRes++;
			}
		}
		// At this point all the atoms should be initialized in our
		// new residue, and the residue should be in the chain, which
		// is in the protein.
		if (numAtomsInRes < pTheResidue->getNumAtoms())
		{
			cout << "Too few atoms in residue ";
			cout << pTheResidue->getResNum() ;
			cout << " - expected: " << pTheResidue->getNumAtoms();
			cout << " found: " <<  numAtomsInRes << endl;
			pCurrentChain->fixBrokenResidue(pCurrentChain->getNumResidues()-1);
		}
		else if (numAtomsInRes > pTheResidue->getNumAtoms())
		{
			cout << "Too many atoms in residue ";
			cout << pTheResidue->getResNum() ;
			cout << " - expected: " << pTheResidue->getNumAtoms();
			cout << " found: " << numAtomsInRes << endl;
			//pCurrentChain->fixBrokenResidue(pCurrentChain->getNumResidues()-1);
		}
	} // end loop over residues

	parseBIOMT();

	pTheProtein->finishProteinBuild();

	if (!pItsEnsemble)
	{	pItsEnsemble = new ensemble();
	}
	pItsEnsemble->add(pTheProtein);
}

/*
int where;
where = string1.find(string2);
// const npos in the string class is a non-valid position
if  (where == string::npos)
// means i didn't find it!
*/
void PDBInterface::parseBIOMT()
{
// NOT YET FINISHED! c.s. 12/5/01
	UInt numRemarkLines = remarkLines.size();
	string biomtstring = "REMARK 350";
	bool startFound = false;
	UInt start = 0;
	bool endFound = false;
	UInt end = 0;
	for (UInt i=0; i<numRemarkLines; i++)
	{
		if (theLines[remarkLines[i]].substr(0,10) == biomtstring)
		{	// this means i've found it
			if (!startFound)
			{
				startFound = true;
				start = i;
			}
		}
		else
		{	if (startFound)
			{
				endFound = true;
				end = i-1;
				break;
			}
		}
	}
	if (startFound && endFound)
	{//	cout << "Remark 350 Lines:" << endl;
		for (UInt i=start; i<=end; i++)
		{	//cout << theLines[remarkLines[i]] << endl;
		}
		//cout << endl;
	}
	else if (!startFound and !endFound)
	{	//cout << "No Remark 350 record found" << endl;
		return;
	}
	else
	{	cout << "Error reading Remark 350 records" << endl;
		return;
	}
	vector<UInt> applyLines;
	for (UInt i=start; i<=end; i++)
	{
		UInt where;
		where = theLines[remarkLines[i]].find("APPLY");
		if  (where != string::npos)
		{	applyLines.push_back(i);
		}
	}
}
