//*******************************************************************************************************
//*******************************************************************************************************
//**********************************                            *****************************************
//**********************************  		solvation 1.0 		*****************************************
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
	string infile = argv[1];
	PDBInterface* thePDB = new PDBInterface(infile);
	ensemble* theEnsemble = thePDB->getEnsemblePointer();
	molecule* pMol = theEnsemble->getMoleculePointer(0);
	protein* bundle = static_cast<protein*>(pMol);
	bundle->silenceMessages();
	residue::setCutoffDistance(9.0);
	pmf::setScaleFactor(0.0);
	rotamer::setScaleFactor(0.0);
	microEnvironment::setScaleFactor(0.0);
	amberVDW::setScaleFactor(1.0);
	amberVDW::setRadiusScaleFactor(1.0);
	amberVDW::setLinearRepulsionDampeningOff();
	amberElec::setScaleFactor(1.0);
	solvation::setItsScaleFactor(0.0);

	double Energy = bundle->intraSoluteEnergy(true);
	pdbWriter(bundle, infile);
	return 0;
}
