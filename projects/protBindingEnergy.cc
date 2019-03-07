//*******************************************************************************************************
//*******************************************************************************************************
//**********************************                            *****************************************
//**********************************    protBindingEnergy 1.0   *****************************************
//**********************************                            *****************************************
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
	//--Running parameters
	if (argc !=2)
	{
        cout << "protBindingEnergy <inFile.pdb>" << endl;
		exit(1);
	}
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
    vector <double> Energies = bundle->protLigandBindingEnergy(0,12);
    cout <<  infile << " " << Energies[0] << " " << Energies[1] << endl;
	return 0;
}
