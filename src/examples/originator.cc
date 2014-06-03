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

	// center on backbone centroid
	dblVec center = prot->getBackBoneCentroid();
	center = center * -1.0;
	prot->translate(center);
	string outfile = argv[2];
	pdbWriter(prot, outfile);
	return 0;
}



