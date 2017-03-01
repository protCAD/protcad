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

//--Program setup-------------------------------------------------------------
int main (int argc, char* argv[])
{
	if (argc !=2)
	{cout << "protEnergy <inFile.pdb> "<< endl;
		exit(1);}

	string infile = argv[1];
	PDBInterface* thePDB = new PDBInterface(infile);
	ensemble* theEnsemble = thePDB->getEnsemblePointer();
	molecule* pMol = theEnsemble->getMoleculePointer(0);
	protein* bundle = static_cast<protein*>(pMol);
	residue::setCutoffDistance(8.0);
	residue::setTemperature(300);
	residue::setElectroSolvationScaleFactor(1.0);
	residue::setHydroSolvationScaleFactor(1.0);
	amberElec::setScaleFactor(1.0);
	amberVDW::setScaleFactor(1.0);

	double Energy = bundle->protEnergy();
	cout << "Energy: " << Energy << " " << infile << endl;
	return 0;
}
