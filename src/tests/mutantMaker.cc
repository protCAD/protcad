#include <iostream>
#include <string>
#include "ensemble.h"
#include "PDBInterface.h"

//--Program setup--------------------------------------------------------------
int main (int argc, char* argv[])
{
    if (argc !=3)
    {
	cout << "mutantMaker <inFile.pdb> <outFile.pdb>" << endl;
	exit(1);
    }

    enum aminoAcid {A, R, N, D, C, Q, E, G, H, I, L, K, M, F, P, S, T, W, Y, V};
    string infile = argv[1];
    PDBInterface* thePDB = new PDBInterface(infile);
    ensemble* theEnsemble = thePDB->getEnsemblePointer();
    molecule* pMol = theEnsemble->getMoleculePointer(0);
    protein* bundle = static_cast<protein*>(pMol);
    bundle->silenceMessages();
    residue::setCutoffDistance(8.0);
    pmf::setScaleFactor(0.0);
    rotamer::setScaleFactor(1.0);
    microEnvironment::setScaleFactor(0.0);
    amberVDW::setScaleFactor(1.0);
    amberVDW::setRadiusScaleFactor(0.95);
    amberVDW::setLinearRepulsionDampeningOff();
    amberElec::setScaleFactor(0.0);
    solvation::setItsScaleFactor(0.0);

//--Mutate sequence of inFile---------------------------------------------------
    UInt chainNum = bundle->getNumChains();
    int resID[] = {L,E,K,S,I,D,D,L,E,D,E,L,Y,A,Q,K,L,K,Y,K,A,I,S,E,E,L,D,H,A,L,N,D,M,T,S,I};
    int arraySize = sizeof(resID)/sizeof(resID[0]);
    for (UInt i = 0; i < chainNum; i++)
    {
	    for (int j = 0; j < arraySize; j ++)
	    {	
			bundle->activateForRepacking(i,j);
			bundle->mutate(i, j, (UInt)resID[j]);	
	    }
    }
    
//--Write to file---------------------------------------
    string outFile = argv[2];
    pdbWriter(bundle, outFile);
    return 0;
}

