//*******************************************************************************************************
//*******************************************************************************************************
//**************************************                       ******************************************
//**************************************     acidMutator 1.0   ******************************************
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
	if (argc !=3)
	{
	cout << "acidMutator <inFile.pdb> <outFile.pdb>" << endl;
	exit(1);
	}

    enum aminoAcid {A,R,N,D,Dh,C,Cx,Q,E,Eh,Hd,He,Hn,Hp,I,L,K,M,F,P,O,S,T,W,Y,V,G,dA,dR,dN,dD,dDh,dC,dCx,dQ,dE,dEh,dHd,dHe,dHn,dHp,dI,dL,dK,dM,dF,dP,dO,dS,dT,dW,dY,dV,HC};
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


	//--Mutate sequence of inFile-----------------------------------------------------------------------
	UInt chainNum = bundle->getNumChains();
    int resID[] = {A, Dh, Eh, dDh, dEh};
	for (UInt i = 0; i < chainNum; i ++)
	{
		UInt resNum = bundle->getNumResidues(i);
		for (UInt j = 0; j < resNum; j ++)
		{	
			UInt restype = bundle->getTypeFromResNum(i, j);
			bundle->activateForRepacking(i, j);
            if (restype != G)
			{
				bundle->mutate(i, j, resID[0]);
			}
		}
	}


	//--Write to file-----------------------------------------------------------------------------------
	cout << endl << "Mutated!!" << endl << endl;
	string outFile = argv[2];
	pdbWriter(bundle, outFile);
	return 0;
}

