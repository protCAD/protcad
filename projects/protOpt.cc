//*******************************************************************************************************
//*******************************************************************************************************
//*************************************                      ********************************************
//*************************************       protOpt        ********************************************
//*************************************                      ********************************************
//*******************************************************************************************************
//******** -sidechain and backbone optimization with a burial-based scaling of electrostatics- **********
//*******************************************************************************************************

/////// Just specify a infile and outfile, it will optimize to a generally effective minimum.

#include "ensemble.h"
#include "PDBInterface.h"
#include <sstream>
#include <time.h>
int main (int argc, char* argv[])
{
    if (argc !=3)
    {   cout << "protOpt <inFile.pdb> <outFile.pdb>" << endl;
        exit(1); }

    clock_t t;
    string infile = argv[1];
    string outFile = argv[2];
    PDBInterface* thePDB = new PDBInterface(infile);
    ensemble* theEnsemble = thePDB->getEnsemblePointer();
    molecule* pMol = theEnsemble->getMoleculePointer(0);
    protein* _prot = static_cast<protein*>(pMol);

    bool backbone = false;
    bool homoSymmetric = true;
    UInt _frozenResidues[] = {3,5,6,13};

    UIntVec frozenResidues;
    UInt frozenResiduesSize = sizeof(_frozenResidues)/sizeof(_frozenResidues[0]);
    for (UInt i = 0; i < frozenResiduesSize; i++)
    {
        frozenResidues.push_back(_frozenResidues[i]);
    }
    if (homoSymmetric)
    {
        for (UInt i = 1; i < _prot->getNumChains(); i++)
        {
            _prot->symmetryLinkChainAtoB(i, 0);
        }
    }

    t=clock();
    _prot->protOpt(backbone, frozenResidues, 0);
    t=clock()-t;
    cout << "Time: " << ((float)t)/CLOCKS_PER_SEC << " Energy: " << _prot->protEnergy() << endl;
    pdbWriter(_prot, outFile);

    return 0;
}

