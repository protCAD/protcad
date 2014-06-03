#include <iostream>
#include <string>
#include <vector>
#include "ensemble.h"
#include "PDBInterface.h"
#include "ran.h"
#include "generalio.h"

static double PI=3.14159265;

int main(int argc, char* argv[])
{
	enum aminoAcid {A, R, N, D, C, Q, E, G, H, I, L, K, M, F, P, S, T, W, Y, V};

	if (argc < 3)
	{
		cout << "program input.pdb  output.pdb" << endl;	
		exit(1);
	}

	// get name of input file
	string inputFileName = argv[1];
	string outputFileName = argv[2];

	// convert inFile to protein object 
	PDBInterface* thePDB = new PDBInterface(inputFileName);
	ensemble* theEnsemble = thePDB->getEnsemblePointer();
	molecule* pMol = theEnsemble->getMoleculePointer(0);
	protein* myProt = static_cast<protein*>(pMol);

	myProt->activateAllForRepacking(0);

	myProt->mutate(0, i, resType);
	myProt->mutate(0, i, resType);
	myProt->mutate(0, i, resType);
	myProt->mutate(0, i, resType);
	myProt->mutate(0, i, resType);
	myProt->mutate(0, i, resType);
	myProt->mutate(0, i, resType);
	myProt->mutate(0, i, resType);
	myProt->mutate(0, i, resType);
	myProt->mutate(0, i, resType);
	
	myProt->mutate(0, i, resType);
	myProt->mutate(0, i, resType);
	myProt->mutate(0, i, resType);
	myProt->mutate(0, i, resType);
	myProt->mutate(0, i, resType);
	myProt->mutate(0, i, resType);
	myProt->mutate(0, i, resType);
	myProt->mutate(0, i, resType);
	myProt->mutate(0, i, resType);
	myProt->mutate(0, i, resType);

	myProt->mutate(0, i, resType);
	myProt->mutate(0, i, resType);
	myProt->mutate(0, i, resType);
	myProt->mutate(0, i, resType);
	myProt->mutate(0, i, resType);
	myProt->mutate(0, i, resType);
	myProt->mutate(0, i, resType);
	myProt->mutate(0, i, resType);
	myProt->mutate(0, i, resType);
	myProt->mutate(0, i, resType);

	myProt->mutate(0, i, resType);
	myProt->mutate(0, i, resType);
	myProt->mutate(0, i, resType);
	myProt->mutate(0, i, resType);
	myProt->mutate(0, i, resType);
	myProt->mutate(0, i, resType);
	myProt->mutate(0, i, resType);
	myProt->mutate(0, i, resType);
	myProt->mutate(0, i, resType);
	myProt->mutate(0, i, resType);

	myProt->mutate(0, i, resType);
	myProt->mutate(0, i, resType);
	myProt->mutate(0, i, resType);
	myProt->mutate(0, i, resType);
	myProt->mutate(0, i, resType);
	myProt->mutate(0, i, resType);
	myProt->mutate(0, i, resType);
	myProt->mutate(0, i, resType);
	myProt->mutate(0, i, resType);




	// set up the energy function
	residue::setCutoffDistance(10.0);
	pmf::setScaleFactor(0.0);
	rotamer::setScaleFactor(1.0);
	microEnvironment::setScaleFactor(0.0);
	amberVDW::setScaleFactor(1.0);
	amberVDW::setRadiusScaleFactor(1.0);
	amberVDW::setLinearRepulsionDampeningOff();
	amberElec::setScaleFactor(0.0);
	solvation::setItsScaleFactor(0.0);

	pdbWriter(myProt, outputFileName);

	return 0;

}

sscanf(argv[2], "%lf", &radiusMin);
