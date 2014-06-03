#include <iostream>
#include <string>
#include "ensemble.h"
#include "PDBInterface.h"

int main (int argc, char* argv[])
{
	string infile = argv[1];
	// read in protein
	PDBInterface* thePDB = new PDBInterface(infile);
	ensemble* theEnsemble = thePDB->getEnsemblePointer();
	molecule* pMol = theEnsemble->getMoleculePointer(0);
	protein* prot = static_cast<protein*>(pMol);
	prot->silenceMessages();
	prot->activateAllForRepacking(0);
//	prot->setCanonicalHelixRotamersOnly(0);
	pdbWriter(prot, "out1.pdb");
	prot->saveCurrentState();
	prot->optimizeRotamers();
	pdbWriter(prot, "out2.pdb");
	//prot->optimizeSmallRotations(2,30.0);
	//pdbWriter(prot, "out3.pdb");
	prot->undoState();
	pdbWriter(prot, "out4.pdb");
	return 0;
}



