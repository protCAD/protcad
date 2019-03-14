//*******************************************************************************************************
//*******************************************************************************************************
//************************************                  *************************************************
//************************************   protMover 1.0  *************************************************
//************************************                  *************************************************
//*******************************************************************************************************
//***************************** -move protein on all six axes- ******************************************
//*******************************************************************************************************

////  Just specify a infile, translation distances (A) and rotation angles for each axes, and outfile.

//--Program setup----------------------------------------------------------------------------------------

#include <iostream>
#include <string>
#include "ensemble.h"
#include "PDBInterface.h"

int main(int argc, char* argv[])
{
	if (argc !=9)
	{
		cout << "protMover <inFile.pdb> <translate-x> <translate-y> <translate-z> <rotate-x> <rotate-y> <rotate-z> <outFile.pdb>" << endl;
		exit(1);
	}

	string inFile = argv[1];
	PDBInterface* thePDB = new PDBInterface(inFile);
    ensemble* theEnsemble = thePDB->getEnsemblePointer();
    molecule* pMol = theEnsemble->getMoleculePointer(0);
	protein* bundle = static_cast<protein*>(pMol);



	double transx, transy, transz, rotx, roty, rotz;

	sscanf(argv[2], "%lf", &transx);
	sscanf(argv[3], "%lf", &transy);
	sscanf(argv[4], "%lf", &transz);
	sscanf(argv[5], "%lf", &rotx);
	sscanf(argv[6], "%lf", &roty);
	sscanf(argv[7], "%lf", &rotz);
	
	if (rotx != 0)
	{
		bundle->rotate(X_axis, rotx);
	}
	if (roty != 0)
	{
		bundle->rotate(Y_axis, roty);
	}
	if (rotz != 0)
	{
		bundle->rotate(Z_axis, rotz);
	}

	bundle->translate(transx, transy, transz);

	string outFile = argv[8];
	pdbWriter(bundle, outFile);

	cout << endl << "Moved!!" << endl << endl;

	return 0;
}
