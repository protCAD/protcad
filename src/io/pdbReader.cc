#include "protein.h"
#include "pdbReader.h"
#include <iostream>
#include <fstream>

molecule* pdbReader(const string& _pdbFile)
{	ifstream inFile(_pdbFile.c_str());
	if(!inFile)
	{	cout << "ERROR: cannot open file: " << _pdbFile << endl;
		return 0;
	}

	//Initialize the protein structure (this may need to be modified
	//in future to account for multiple chains,  etc...	
	protein* pTheProtein = new protein;
	// implementation of multi-chain proteins
	char currentChainID = ' ';
	chain* pTheChain = new chain(currentChainID);
	pTheProtein->add(pTheChain);

	// Initializing the counters, etc.
	int currentResSeq = -1;
	unsigned int skip = 0;
	unsigned int firstAtom = 1;
	// Creating dummy residue to initialize residue type-base
	residue* tempRes = new residue(0);
	delete tempRes;
	// Reset tempRes pointer to 0 for further use
	tempRes = 0;
	//unsigned int ResType;
	unsigned int numatomsInRes = 0;

	pdbAtom theData;
#ifdef PDBREADER_DEBUG
	//cout << "starting to read file.." << endl;
#endif
	while (theData.pdbGetLine(inFile))
	{	if (firstAtom)
		{	currentChainID = theData.getChainID();
		//	cout << "theData ID;" << theData.getChainID() << endl;
			pTheChain->setChainID(theData.getChainID());
		//	cout << "theChain ID:" << pTheChain->getChainID() << endl;
			firstAtom = 0;
	   	}
		
		//cout << theData.getChainID() << "  " << currentChainID << endl;
		if ( theData.getChainID() != currentChainID)
		{	pTheChain = new chain(theData.getChainID());
			pTheProtein->add(pTheChain);
			currentChainID = theData.getChainID();
		}

		//Initially, the currentResSeq is set to -1 to enable entry into loop
		//and building of first residue.  Thereafter, it's reset to the actual
		//residue number read in from the pdb file
		if ( theData.getResSeq() != currentResSeq)
		{	// Check to see that the number of atoms in the residue we just
			// built is correct.  If not, fix it!
			if (tempRes !=0)
			{  // We've actually got a residue to work with!
				if (numatomsInRes < tempRes->getNumAtoms())
				{	cout << "Residue Number " << currentResSeq << " is broken." << endl;
					cout << "TOO FEW ATOMS" << endl << endl;
					if (numatomsInRes < residue::dataBase[tempRes->getTypeIndex()].mainChain.size())
						{  cout << "Residue cannot be fixed : missing main chain atoms" << endl;
						 	exit(1);
						}
					else
						pTheChain->fixBrokenResidue(pTheChain->getNumResidues()-1);	
						tempRes = 0;
				}
				else if (numatomsInRes > tempRes->getNumAtoms())
				{
					cout << "Residue Number " << currentResSeq << " is broken." << endl;
					cout << "TOO MANY ATOMS" << endl << endl;
					pTheChain->fixBrokenResidue(pTheChain->getNumResidues()-1);	
					tempRes = 0;
				}
			}
			numatomsInRes = 0;
			skip = 0;
			unsigned int nameFound = 0;
			for (unsigned int i=0 ; i<residue::getDataBaseSize() ; i++)
			{	if ( theData.getItem(resName) == residue::getDataBaseItem(i) )
				{	tempRes = new residue(i);
					pTheChain->add(tempRes);
					nameFound = 1;
					//ResType = i;
					currentResSeq =  theData.getResSeq();
					tempRes->setResNum(currentResSeq);
					break;
				}
			}
			if (!nameFound)
			{	cout << theData.getItem(resName) << " cannot be made " << endl;
				cout << "Do you want to continue? (y/n)";
				char tempchar;
				cin >> tempchar;
				if (tempchar == 'y')
				{	skip = 1;
				}
				else
				{	exit(0);
				}
			}
			if (!skip)
			{	// add the atom to the residue
				tempRes->addAtom(theData);
				numatomsInRes++;
			}
		}
		else
		{	if (!skip)
			{	// add the atom to the residue
				tempRes->addAtom(theData);
				numatomsInRes++;
			}
		}
	}	
	if (currentResSeq == -1)
	{	cout << "Error reading pdbfile " << _pdbFile << endl;
		cout << "File appears to be empty!" << endl;
		inFile.close();
		inFile.clear();
		exit(1);
	}	
	inFile.close();
	inFile.clear();
	static_cast<protein*>(pTheProtein)->finishProteinBuild();
	return pTheProtein;
}
