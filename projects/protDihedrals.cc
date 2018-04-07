//*******************************************************************************************************
//*******************************************************************************************************
//**************************************                       ******************************************
//**************************************   getDihedrals 1.0    ******************************************
//**************************************                       ******************************************
//*******************************************************************************************************
//*******************************************************************************************************


//--Program setup----------------------------------------------------------------------------------------
#include <iostream>
#include <string>
#include "ensemble.h"
#include "PDBInterface.h"
int main (int argc, char* argv[])
{
    if (argc !=2)
	{
    cout << "getDihedrals <inFile.pdb>" << endl;
	exit(1);
	}

	enum aminoAcid {A, R, N, D, Dh, C, Q, E, Eh, G, H, I, L, K, M, F, P, S, T, W, Y, V, dA, dR, dN, dD, dDh, dC, dQ, dE, dEh, dH, dI, dL, dK, dM, dF, dP, dS, dT, dW, dY, dV};
	string infile = argv[1];
	PDBInterface* thePDB = new PDBInterface(infile);
	ensemble* theEnsemble = thePDB->getEnsemblePointer();
	molecule* pMol = theEnsemble->getMoleculePointer(0);
    protein* bundle = static_cast<protein*>(pMol);

	//cout << endl << "chain res phi psi" << endl;


	//--Search sequence of inFile for phis >= 0 --------------------------------------------------------
	UInt chainNum = bundle->getNumChains();	
	for (UInt i = 0; i < chainNum; i ++)
	{
		UInt resNum = bundle->getNumResidues(i);
		if (resNum > 1)
		{
			for (UInt j = 1; j < resNum-1; j ++)
			{
				double phi = bundle->getPhi(i,j);
				double psi = bundle->getPsi(i,j);
				cout << phi << " " << psi << endl;
			}
		}
	}
	return 0;
}

