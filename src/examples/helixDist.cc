#include <iostream>
#include <string>
#include "ensemble.h"
#include "PDBInterface.h"
int main (int argc, char* argv[])
{
    	if (argc != 2)
    	{
		cout << "helixDist  <infilename.pdb>" << endl;
    	exit(1);
    	}
	
    	string infile = argv[1];
    	PDBInterface* thePDB = new PDBInterface(infile);
    	ensemble* theEnsemble = thePDB->getEnsemblePointer();
    	molecule* pMol = theEnsemble->getMoleculePointer(0);

    	protein* bundle = static_cast<protein*>(pMol);

					
		dblVec coords1 = bundle->getCoords(0,10,"CB");
		dblVec coords2 = bundle->getCoords(0,59,"CB");

		coords1[2] = 0.0;
		coords2[2] = 0.0;

		double dist = CMath::distance(coords1, coords2);

		cout << " dist " << dist <<endl;

	return 0;
}

