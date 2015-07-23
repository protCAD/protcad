//*******************************************************************************************************
//*******************************************************************************************************
//**************************************                       ******************************************
//**************************************   structFinder 1.0    ******************************************
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
	cout << "StructFinder <inFile.pdb>" << endl;
	exit(1);
	}

	enum aminoAcid {A, R, N, D, Dh, C, Q, E, Eh, G, H, I, L, K, M, F, P, S, T, W, Y, V, dA, dR, dN, dD, dDh, dC, dQ, dE, dEh, dH, dI, dL, dK, dM, dF, dP, dS, dT, dW, dY, dV};
	string infile = argv[1];
	PDBInterface* thePDB = new PDBInterface(infile);
	ensemble* theEnsemble = thePDB->getEnsemblePointer();
	molecule* pMol = theEnsemble->getMoleculePointer(0);
	protein* bundle = static_cast<protein*>(pMol);
	bundle->silenceMessages();
	residue::setCutoffDistance(8.0);
	rotamer::setScaleFactor(0.0);
	amberVDW::setScaleFactor(1.0);
	amberVDW::setRadiusScaleFactor(0.95);
	amberVDW::setLinearRepulsionDampeningOff();
	amberElec::setScaleFactor(0.0);
	
	cout << endl << "\t*StructFinder*" << endl;
	cout << endl <<"Generate phi psi angles for modeling" << endl;
	cout << endl << "chain res phi psi" << endl;


	//--Search sequence of inFile for phis >= 0 --------------------------------------------------------
	UInt chainNum = bundle->getNumChains();	
	for (UInt i = 0; i < chainNum; i ++)
	{
		UInt resNum = bundle->getNumResidues(i);
		for (UInt j = 0; j < resNum; j ++)
		{	
			double phi = bundle->getPhi(i,j);
			double psi = bundle->getPsi(i,j);
			cout << i+1 << " " << j+1 << " " << phi << " " << psi << endl;
		}
	}
	return 0;
}

