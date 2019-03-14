#include <iostream>
#include <string>
#include "typedef.h"
#include "PDBInterface.h"
#include "pdbWriter.h"
#include "atomIterator.h"


int main(int argc, char* argv[])
{	
	if (argc != 2)
	{	cout << "Usage: PDBInterfacetest pdbname" << endl;
		exit (1);
	}	
		
	string fname = argv[1];
	PDBInterface* thePDB = new PDBInterface(fname);
	ensemble* theEnsemble = thePDB->getEnsemble();

	molecule* pTheMolecule = theEnsemble->getMoleculePointer(0);
	if (pTheMolecule == 0)
	{	cout << "Failed" << endl;
		return 1;
	}

	protein* pTheProtein = static_cast<protein*>(pTheMolecule);
	axis theAxis;

	theAxis = Z_axis;
	pTheProtein->rotate(0,theAxis,-160.0);

	pTheProtein->translate(0,9.0,0.0,0.0);
	//theAxis = X_axis;
	//pTheProtein->rotate(0,theAxis,5.0);
	theAxis = Z_axis;
	pTheProtein->rotate(0,theAxis,-45.0);

	string newfilename = fname + "_panel2";

	int success = pdbWriter(pTheProtein,newfilename);
	
	if (success == 1)
	{	cout << "Success!" << endl;
	}


	// OK, now lets get the chain pointer and duplicate
	chain* pTheChain = pTheProtein->getChain(0);

	chain* pTheChain2 = new chain(*pTheChain);
	chain* pTheChain3 = new chain(*pTheChain);
	chain* pTheChain4 = new chain(*pTheChain);

	char id = 'A';
	pTheChain->setChainID(id);
	id = 'B';
	pTheChain2->setChainID(id);
	id = 'C';
	pTheChain3->setChainID(id);
	id = 'D';
	pTheChain4->setChainID(id);

	pTheProtein->add(pTheChain2);
	pTheProtein->add(pTheChain3);
	pTheProtein->add(pTheChain4);

	theAxis = Z_axis;
	pTheProtein->rotate(1,theAxis,180.0);
	theAxis = X_axis;
	pTheProtein->rotate(2,theAxis,180.0);
	theAxis = Y_axis;
	pTheProtein->rotate(3,theAxis,180.0);

	newfilename = fname + "_panel3";

	success = pdbWriter(pTheProtein,newfilename);
	
	if (success == 1)
	{	cout << "Success!" << endl;
	}

	// OK, now let's coil this protein up
	double zmin = 100.0;

	atomIterator it1(pTheProtein);

	atom* tempAtom = 0;
	for (;!(it1.last());it1++)
	{	tempAtom = it1.getAtomPointer();
		if ( (tempAtom->getCoords())[2] < zmin )
		{	zmin = (tempAtom->getCoords())[2];
		}
	}
	double pitch = 189.0;
	double twopi = 2.0 * 3.14159;
	double theta;
	double sth;
	double cth;
	it1.initialize();
	for (;!(it1.last());it1++)
	{	tempAtom = it1.getAtomPointer();
		dblVec tempCoords = tempAtom->getCoords();
		theta = (-(tempCoords[2]-zmin)/pitch)*twopi;
		sth = sin(theta);
		cth = cos(theta);
		dblVec newCoords(tempCoords);
		newCoords[0] = tempCoords[0]*cth - tempCoords[1]*sth;
		newCoords[1] = tempCoords[0]*sth + tempCoords[1]*cth;
		newCoords[2] = tempCoords[2];
		tempAtom->setCoords(newCoords[0],newCoords[1],newCoords[2]);
	}

	newfilename = fname + "_panel4";

	success = pdbWriter(pTheProtein,newfilename);
	
	if (success == 1)
	{	cout << "Success!" << endl;
	}

	return 0;
}
