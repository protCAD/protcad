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
	string infile = argv[1];

	PDBInterface* thePDB = new PDBInterface(infile);
	ensemble* theEnsemble = thePDB->getEnsemblePointer();
	molecule* pMol = theEnsemble->getMoleculePointer(0);
	protein* prot = static_cast<protein*>(pMol);
	
	amberElec::setScaleFactor(1.0);
	amberVDW::setScaleFactor(1.0);
	residue::setElectroSolvationScaleFactor(1.0);
	residue::setHydroSolvationScaleFactor(1.0);
	residue::setPolarizableElec(false);
	residue::setEntropyFactor(0.0);
	
#ifdef __CUDA__
	prot->loadDeviceMemAll();
	cout << infile << " " << prot->getNumClashesCU() << " clashes " << prot->protEnergyCU() << " kcal/mol at 300K" << endl;
#else
	cout << infile << " " << prot->getNumHardClashes() << " clashes " << prot->protEnergy() << " kcal/mol at 300K" << endl;
#endif
	//pdbWriter(prot, infile);
	return 0;
}
