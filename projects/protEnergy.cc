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
    clock_t t;
    string infile = argv[1];
    PDBInterface* thePDB = new PDBInterface(infile);
    ensemble* theEnsemble = thePDB->getEnsemblePointer();
    molecule* pMol = theEnsemble->getMoleculePointer(0);
    protein* bundle = static_cast<protein*>(pMol);
    residue::setCutoffDistance(9.0);
    residue::setTemperature(300);
    residue::setElectroSolvationScaleFactor(1.0);
    residue::setHydroSolvationScaleFactor(1.0);
    amberElec::setScaleFactor(0.0);
    amberVDW::setScaleFactor(1.0);
    t=clock();
    double Energy = bundle->protEnergy();
    //cout << bundle->getNumAtoms(0,0) << " " << infile << " ";
    //bundle->mutateWBC(0,0,23);
    t=clock()-t;
    cout << Energy << " " << infile << " ";
    pdbWriter(bundle, infile);
	return 0;
}
