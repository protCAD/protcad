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
    for (int i = 0; i<20; i++ )
    {

        PDBInterface* thePDB = new PDBInterface(infile);
        ensemble* theEnsemble = thePDB->getEnsemblePointer();
        molecule* pMol = theEnsemble->getMoleculePointer(0);
        protein* bundle = static_cast<protein*>(pMol);
        double cutoff = 1+ (i*0.5);
        residue::setCutoffDistance(cutoff);
        amberVDW::setScaleFactor(1.0);
        amberVDW::setRadiusScaleFactor(1.0);
        amberElec::setScaleFactor(0.0);
        t=clock();
        double Energy = bundle->intraSoluteEnergy(false);
        t=clock()-t;
        cout << "VDW" << " " << cutoff << " " << ((float)t)/CLOCKS_PER_SEC << " " << Energy << endl;
        delete thePDB;
    }
    for (int i = 0; i<40; i++ )
    {
        PDBInterface* thePDB = new PDBInterface(infile);
        ensemble* theEnsemble = thePDB->getEnsemblePointer();
        molecule* pMol = theEnsemble->getMoleculePointer(0);
        protein* bundle = static_cast<protein*>(pMol);
        double cutoff = 1+ (i*0.5);
        residue::setCutoffDistance(cutoff);
        amberVDW::setScaleFactor(0.0);
        amberVDW::setRadiusScaleFactor(0.0);
        amberElec::setScaleFactor(1.0);
        t=clock();
        double Energy = bundle->intraSoluteEnergy(true);
        t=clock()-t;
        cout << "Elec" << " " << cutoff << " " << ((float)t)/CLOCKS_PER_SEC << " " << Energy << endl;
        delete thePDB;
    }


    //pdbWriter(bundle, infile);
	return 0;
}
