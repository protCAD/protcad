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
#include <sstream>
#include <time.h>
int main (int argc, char* argv[])
{
    if (argc !=2)
    {   cout << "protOpt <inFile.pdb>" << endl;
        exit(1); }

    string infile = argv[1];
    string outFile;
    for (UInt i = 0; i < 100; i++)
    {
        PDBInterface* thePDB = new PDBInterface(infile);
        ensemble* theEnsemble = thePDB->getEnsemblePointer();
        molecule* pMol = theEnsemble->getMoleculePointer(0);
        protein* _prot = static_cast<protein*>(pMol);

        clock_t t;
        t=clock();
        _prot->protOpt(true);
        t=clock()-t;
        cout << i+1 << " " << ((float)t)/CLOCKS_PER_SEC << " " << _prot->protEnergy() << endl;
        stringstream convert;
        string countstr;
        convert << i+1, countstr = convert.str();
        outFile = countstr + ".pdb";
        pdbWriter(_prot, outFile);
        delete thePDB;
    }

    return 0;
}

