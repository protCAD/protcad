#include "typedef.h"
#include "ensemble.h"
#include "PDBInterface.h"
#include "ran.h"
#include <sstream>
#include<cstdlib>
#include<ctime>



int main (int argc, char* argv[])
{

        enum aminoAcid {A, R, N, D, C, Q, E, G, H, I, L, K, M, F, P, S, T, W, Y, V, X};

	string inputFileName = argv[1];

       	// read in prot structure
        PDBInterface* thePDB = new PDBInterface(inputFileName);
        ensemble* theEnsemble = thePDB->getEnsemblePointer();
        molecule* theMol = theEnsemble->getMoleculePointer(0);
        protein* prot = static_cast<protein*>(theMol);

	prot->silenceMessages();
	prot->activateAllForRepacking(0);


	for (UInt i = 0; i < prot->getNumResidues(0); i ++)
	{
		UInt ID;
		prot->setPhi(0,i, 180.0);
		prot->setPsi(0,i, 180.0);
		sscanf (argv[3 + i], "%u", &ID);
		prot->mutate(0,i,ID);
	}

	string outFile = argv[2];

	pdbWriter (prot, outFile);
	return 0;
}





