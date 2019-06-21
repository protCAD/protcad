//*******************************************************************************************************
//*******************************************************************************************************
//**************************************                       ******************************************
//**************************************   protDihedrals 1.0    ******************************************
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
    cout << "protDihedrals <inFile.pdb>" << endl;
	exit(1);
	}
	
	string backboneTypes[] = {"-π","-α","-ρ","-β","β","ρ","α","π","-πi","-αi","-ρi","-βi","βi","ρi","αi","πi"};
	string infile = argv[1];
	PDBInterface* thePDB = new PDBInterface(infile);
	ensemble* theEnsemble = thePDB->getEnsemblePointer();
	molecule* pMol = theEnsemble->getMoleculePointer(0);
    protein* bundle = static_cast<protein*>(pMol);

	cout << endl << "phi psi backbonetype" << endl;

	UInt chainNum = bundle->getNumChains();
	for (UInt i = 0; i < chainNum; i ++)
	{
		UInt resNum = bundle->getNumResidues(i);
		for (UInt j = 1; j < resNum-2; j ++)
		{
			UInt backboneType = bundle->getBackboneSequenceType(i,j);
			cout << i << " " << j << " " << bundle->getPhi(i,j) << " " << bundle->getPsi(i,j) << " " << backboneTypes[backboneType] << endl;
		}
	}
	return 0;
}

