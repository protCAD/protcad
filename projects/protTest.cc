//*******************************************************************************************************
//*******************************************************************************************************
//*****************************                         *************************************************
//*****************************        protEnergy       *************************************************
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
		cout << "protEnergy <inFile.pdb>" << endl;
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
	residue::setPolarizableElec(false);
	amberElec::setScaleFactor(1.0);
	amberVDW::setScaleFactor(0.0);
	residue::setTemperature(300);
	

	for (UInt i = 0; i < 200; i++)
	{
		bundle->getSoluteEnergy(0, 1, 3, 1, 2, 4);
		bundle->translateChain(1,0.04,0,0);
	}
	
	return 0;
}
