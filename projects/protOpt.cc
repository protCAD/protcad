//*******************************************************************************************************
//*******************************************************************************************************
//*************************************                      ********************************************
//*************************************       protOpt        ********************************************
//*************************************                      ********************************************
//*******************************************************************************************************
//******** -sidechain and backbone optimization with a burial-based scaling of electrostatics- **********
//*******************************************************************************************************

/////// Just specify a infile and outfile, it will optimize to a generally effective minimum.

#include "ensemble.h"
#include "PDBInterface.h"
#include <sstream>
#include <time.h>
int main (int argc, char* argv[])
{
	if (argc !=3)
	{   cout << "protOpt <inFile.pdb> <outFile.pdb>" << endl;
		exit(1); }

	string infile = argv[1];
	string outFile = argv[2];
	PDBInterface* thePDB = new PDBInterface(infile);
	ensemble* theEnsemble = thePDB->getEnsemblePointer();
	molecule* pMol = theEnsemble->getMoleculePointer(0);
	protein* _prot = static_cast<protein*>(pMol);

	bool backbone = true;
	double startEnergy = _prot->protEnergy();
	time_t start,end;
	time (&start);
	_prot->protOpt(backbone);
	time (&end);
	double endEnergy = _prot->protEnergy();
	cout << _prot->protEnergy() << " " << (endEnergy-startEnergy)/difftime(end,start) << " ";
	pdbWriter(_prot, outFile);

	return 0;
}

