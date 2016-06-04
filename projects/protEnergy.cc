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
    t=clock();
    double Energy = bundle->protEnergy();
    double symEnergy = bundle->intraSoluteEnergy(true, 0);
    t=clock()-t;
    cout << "Time " << ((float)t)/CLOCKS_PER_SEC << " Energy " << Energy << " SymEnergy " << symEnergy << " TotalSym " << symEnergy+symEnergy+symEnergy+symEnergy+symEnergy+symEnergy << endl;
    pdbWriter(bundle, infile);
	return 0;
}
