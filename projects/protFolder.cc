//*******************************************************************************************************
//*******************************************************************************************************
//*************************************                      ********************************************
//*************************************       protFolder     ********************************************
//*************************************                      ********************************************
//*******************************************************************************************************
//******** -sidechain and backbone optimization with a burial-based scaling of electrostatics- **********
//*******************************************************************************************************

/////// Just specify a infile and outfile, it will optimize to a generally effective minimum.

#include "ensemble.h"
#include "PDBInterface.h"
#include <sstream>
#include <time.h>
#include <unistd.h>
int main (int argc, char* argv[])
{
    if (argc !=2)
    {   cout << "protFolder <inFile.pdb>" << endl;
        exit(1); }

    stringstream convert;
    string startstr;
    srand (getpid());
    UInt name = rand() % 100000000;
    convert << name, startstr = convert.str();
    string foldModel = startstr + "_fold.pdb";
    string infile = argv[1];
    PDBInterface* thePDB = new PDBInterface(infile);
    ensemble* theEnsemble = thePDB->getEnsemblePointer();
    molecule* pMol = theEnsemble->getMoleculePointer(0);
    protein* _prot = static_cast<protein*>(pMol);
    double Energy = _prot->protEnergy();
    double newEnergy = Energy;
    UInt nobetter = 0, plateau = 100;
    while (true)
    {
        _prot->protSampling(plateau);
        _prot->protOpt(false);
        newEnergy = _prot->protEnergy();
        if (newEnergy < Energy){
            pdbWriter(_prot, foldModel);
            cout << startstr << " " << newEnergy << endl;
            Energy = newEnergy;
            nobetter = 0;
        }
    }
    return 0;
}
