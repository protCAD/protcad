#include <iostream>
#include <algorithm>
#include <time.h>
#include "ensemble.h"
#include "PDBInterface.h"

	//--Setup up conditions and data metrics--------------------------------------------------------------------------------------------------------------	

	//--Program setup
int main (int argc, char* argv[])
{
	if (argc !=2)
	{
		cout << "counter <inFile.pdb> <residue to count>" << endl;
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
	amberElec::setScaleFactor(1.0);
	solvation::setItsScaleFactor(0.0);
	srand (time(NULL));


	//--Setup and get protein data
	UInt chainNum = bundle->getNumChains();
	UInt resNum = bundle->getNumResidues(0);
	int count = 0;
	for (UInt h = 0; h < chainNum; h++)
	{
		for (UInt i = 0; i < resNum; i++)
		{
			int restype1 = bundle->getTypeFromResNum(h, i);
			if (restype1 == M) 
			{
				count = count + 1;
			}
		}
	}
	cout << count << endl;
	return 0;
}

