#include <iostream>
#include <fstream>
#include <string>
#include "ensemble.h"
#include "PDBInterface.h"
#include "pdbWriter.h"
#include <vector>
int main ()
{
	enum aminoAcid {A,R,N,D,C,Q,E,G,H,I,L,K,M,F,P,S,T,W,Y,V};

	string GFP_Infile = "TMGFP-globalOptSeq434254_superfolder.pdb";
	PDBInterface* thePDB = new PDBInterface(GFP_Infile);
	ensemble* theEnsemble = thePDB->getEnsemblePointer();
	molecule* pMol = theEnsemble->getMoleculePointer(0);
	protein* GFP = static_cast<protein*>(pMol);

	//here's the residue positions targeted for mutation
	int resPosToOptimize[49] =
		{15,  17,  19,  21,  28,  30,  32,
		 34,  36,  41,  43,  45,  47,  93,
		 95,  97,  99, 101, 105, 107, 109,
		111, 113, 118, 120, 122, 124, 126,
		147, 149, 151, 153, 162, 164, 166,
		178, 180, 182, 184, 186, 200, 202,
		204, 206, 208, 219, 221, 223, 225};

    int randLocNumber;



	for(int i = 1; i<= 5000; i++)
    {
        randLocNumber = (rand()%(49)); //between 1 and 49
        GFP->activateForRepacking(0, GFP->getIndexFromResNum(0,(UInt)resPosToOptimize[randLocNumber]));
		GFP->protein::optimizeRotamers(0, GFP->getIndexFromResNum(0,(UInt)resPosToOptimize[randLocNumber]));
		cout << endl << GFP->intraEnergy;
	}



	//write the PDB file
	pdbWriter(GFP,"TMGFP-globalOptSeq434254_superfolder_optimized.pdb");

	return 0;
}


