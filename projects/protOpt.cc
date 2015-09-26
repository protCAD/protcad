//*******************************************************************************************************
//*******************************************************************************************************
//*************************************                      ********************************************
//*************************************    protOptSolvent    ********************************************
//*************************************                      ********************************************
//*******************************************************************************************************
//******** -sidechain and backbone optimization with a burial-based scaling of electrostatics- **********
//*******************************************************************************************************

/////// Just specify a infile and outfile, it will optimize to a generally effective minimum.

#include "ensemble.h"
#include "PDBInterface.h"
#include <time.h>
int main (int argc, char* argv[])
{
    if (argc !=3)
    {   cout << "protOpt <inFile.pdb> <outFile.pdb>" << endl;
        exit(1); }

    string infile = argv[1];
    string outFile = argv[2];
	PDBInterface* thePDB = new PDBInterface(infile);
	ensemble* theEnsemble = thePDB->getEnsemblePointer();
	molecule* pMol = theEnsemble->getMoleculePointer(0);
	protein* _prot = static_cast<protein*>(pMol);

    clock_t t;
    t=clock();
    _prot->protOpt(100, false);
    t=clock()-t;
    cout << "Time for protOpt: " << ((float)t)/CLOCKS_PER_SEC << " Energy: " << _prot->protEnergy() << endl;
    pdbWriter(_prot, outFile);

    return 0;
}

