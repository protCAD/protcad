//*******************************************************************************************************
//*******************************************************************************************************
//*****************************                         *************************************************
//*****************************        protTest         *************************************************
//*****************************                         *************************************************
//*******************************************************************************************************
//*******************************************************************************************************

#include <iostream>
#include <string>
#include "ensemble.h"
#include "PDBInterface.h"
#include <time.h>

//--Program setup-------------------------------------------------------------
int main (int argc, char* argv[])
{	

	if (argc !=2)
	{
		cout << "protTest <inFile.pdb>" << endl;
		exit(1);
	}

	//clock_t start, end;
	//double cpu_time_used;
	string infile = argv[1];
	
	PDBInterface* thePDB = new PDBInterface(infile);
	ensemble* theEnsemble = thePDB->getEnsemblePointer();
	molecule* pMol = theEnsemble->getMoleculePointer(0);
	protein* bundle = static_cast<protein*>(pMol);
	
	residue::setElectroSolvationScaleFactor(0.0);
	residue::setHydroSolvationScaleFactor(0.0);
	residue::setPolarizableElec(true);
	amberElec::setScaleFactor(1.0);
	amberVDW::setScaleFactor(0.0);
	residue::setTemperature(300);
	
	bundle->translateChain(1,-1,0,0);
	for (UInt i = 0; i < 102; i++)
	{
		//start = clock();
		double E = bundle->getSoluteEnergy(0, 1, 3, 1, 2, 4);
		//double E = bundle->getBackboneHBondEnergy(1, 2, 0, 1);
		//end = clock();
		//cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
		cout << E << endl; //" " << cpu_time_used << endl;
		bundle->translateChain(1,0.04,0,0);
		//bundle->rotateChain(0, Z_axis, -1);
	}
	
	return 0;
}
