//*******************************************************************************************************
//*******************************************************************************************************
//**********************************                            *****************************************
//**********************************  foldingBindingEnergy 1.0  *****************************************
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
		cout << "bindingEnergy <inFile.pdb>" << endl;
		exit(1);
	}
    clock_t t;
	string infile = argv[1];
	PDBInterface* thePDB = new PDBInterface(infile);
	ensemble* theEnsemble = thePDB->getEnsemblePointer();
	molecule* pMol = theEnsemble->getMoleculePointer(0);
	protein* bundle = static_cast<protein*>(pMol);
	bundle->silenceMessages();
	residue::setCutoffDistance(9.0);
	rotamer::setScaleFactor(0.0);
	amberVDW::setScaleFactor(1.0);
	amberVDW::setRadiusScaleFactor(1.0);
	amberVDW::setLinearRepulsionDampeningOff();
    amberElec::setScaleFactor(1.0);

    //vector <double> Energy = bundle->chainFoldingBindingEnergy(false);
    //double interEnergy = bundle->interSoluteEnergy(true, 0, 1);
   // cout << Energy[1] << " " << Energy[0];
    t=clock();
    //bundle->buildResidueMatrices();
    t=clock()-t;
    cout << "Time to run: " << ((float)t)/CLOCKS_PER_SEC << endl << endl;
	return 0;
}
