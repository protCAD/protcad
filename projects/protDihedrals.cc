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
	
	string backboneTypes[] = {"-γ","-π","-α","-ρ","-β","β","ρ","α","π","γ","-γi","-πi","-αi","-ρi","-βi","βi","ρi","αi","πi","γi"};
	string infile = argv[1];
	PDBInterface* thePDB = new PDBInterface(infile);
	ensemble* theEnsemble = thePDB->getEnsemblePointer();
	molecule* pMol = theEnsemble->getMoleculePointer(0);
    protein* bundle = static_cast<protein*>(pMol);

	cout << "type phi psi" << endl;

	UInt chainNum = bundle->getNumChains();
	for (UInt i = 0; i < chainNum; i ++)
	{
		UInt resNum = bundle->getNumResidues(i);
		for (UInt j = 1; j < resNum-1; j ++)
		{
			UInt backboneType = bundle->getBackboneSequenceType(i,j);
			double rpt = bundle->getResiduesPerTurn(i,j);
			cout << rpt << " " << backboneTypes[backboneType] << " " << bundle->getPhi(i,j) << " " << bundle->getPsi(i,j) << endl;
		}
	}
	return 0;
}

