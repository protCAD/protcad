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
	string infile = argv[1];
	PDBInterface* thePDB = new PDBInterface(infile);
	ensemble* theEnsemble = thePDB->getEnsemblePointer();
	molecule* pMol = theEnsemble->getMoleculePointer(0);
	protein* bundle = static_cast<protein*>(pMol);
    residue::setCutoffDistance(8.0);
	residue::setTemperature(300);
    residue::setElectroSolvationScaleFactor(0.0);
    residue::setHydroSolvationScaleFactor(0.0);
    amberElec::setScaleFactor(0.0);
    amberVDW::setScaleFactor(1.0);

    cout << bundle->protEnergy() << endl;
    string outFile = "test.pdb";
    pdbWriter(bundle, outFile);
    /*double solventEnergy=0.0;
    UInt chainNum = bundle->getNumChains();
    for (UInt i = 0; i < chainNum; i ++)
    {
        UInt resNum = bundle->getNumResidues(i);
        for (UInt j = 0; j < resNum; j ++)
        {
            solventEnergy += bundle->getSolvationEnergy(i,j);
        }
    }
    cout << solventEnergy*4184*1000/1.60217662e-19/6.022140857e23 << endl;*/
	return 0;
}
