#include <iostream>
#include <string>
#include <vector>
#include "ensemble.h"
#include "PDBInterface.h"
#include "ran.h"
#include "generalio.h"

int main (int argc, char* argv[])
{
	string inputFile = argv[1];
	// read in prot structure
	PDBInterface* thePDB = new PDBInterface(inputFile);
	ensemble* theEnsemble = thePDB->getEnsemblePointer();
	molecule* theMol = theEnsemble->getMoleculePointer(0);
	protein* prot = static_cast<protein*>(theMol);
	/*
	int chainIndexi, resIndexi;
	UInt chainIndex, resIndex;

	sscanf(argv[2], "%d", &chainIndexi);
	sscanf(argv[3], "%d", &resIndexi);
	chainIndex = (UInt)chainIndexi;
	resIndex = (UInt)resIndexi;

	residue* temp = prot->superimposeGLY(chainIndex,resIndex);
	atom* hca1 = temp->getAtom(5);
	atom* hca2 = temp->getAtom(6);

	dblVec HCA1 = hca1->getCoords();
	dblVec HCA2 = hca2->getCoords();

	cout << HCA1[0] << " " << HCA1[1] << " " << HCA1[2] << endl;
	cout << HCA2[0] << " " << HCA2[1] << " " << HCA2[2] << endl;
	*/
	cout << "HBONDENERGY " << prot->calculateHCA_O_hBondEnergy() << endl;

	//delete temp;
	return 0;
}
