#include <iostream>
#include <string>
#include "ensemble.h"
#include "PDBInterface.h"

int main(int argc, char* argv[])
{
	if (argc !=6)
	{
		cout << "dimerBuilder <inFile.pdb> <radius> <faceangle> <pitch> <outFile.pdb>" << endl;
		exit(1);
	}

	//read in the input file
	string inFile = argv[1];
	PDBInterface* thePDB = new PDBInterface(inFile);
    ensemble* theEnsemble = thePDB->getEnsemblePointer();
    molecule* pMol = theEnsemble->getMoleculePointer(0);
	protein* bundle = static_cast<protein*>(pMol);

	double radius, face, pitch;

	sscanf(argv[2], "%lf", &radius);
	sscanf(argv[3], "%lf", &face);
	sscanf(argv[4], "%lf", &pitch);
	
	//anti parallel
	//bundle->rotate(Z_axis, face);
	//bundle->rotate(1, Z_axis, 180);
	//bundle->translate(0.0,radius,0.0);
	//bundle->rotate(1, Z_axis, 180.0);

	//parallel
	bundle->rotate(Z_axis, face);
	bundle->translate(0.0,radius,0.0);
	bundle->rotate(1, Z_axis, 180.0);

	bundle->coilcoil(pitch);

	string outFile = argv[5];
	pdbWriter(bundle, outFile);

	return 0;
}
