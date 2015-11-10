//*******************************************************************************************************
//*******************************************************************************************************
//*************************************                      ********************************************
//*************************************     protInverter     ********************************************
//*************************************                      ********************************************
//*******************************************************************************************************
//********************** -Flip chirality of amino acids and invert dihedrals- ***************************
//*******************************************************************************************************

/////// Just specify a infile and outfile, it will optimize to a generally effective minimum.

#include "ensemble.h"
#include "PDBInterface.h"
int main (int argc, char* argv[])
{
    if (argc !=4)
    {   cout << "protInverter <inFile.pdb> <outFile.pdb>" << endl;
        exit(1); }

    string infile = argv[1];
    string infile2 = argv[2];
    string outFile = argv[3];
    PDBInterface* thePDB = new PDBInterface(infile);
    ensemble* theEnsemble = thePDB->getEnsemblePointer();
    molecule* pMol = theEnsemble->getMoleculePointer(0);
    protein* bundle = static_cast<protein*>(pMol);

    PDBInterface* thePDB2 = new PDBInterface(infile2);
    ensemble* theEnsemble2 = thePDB2->getEnsemblePointer();
    molecule* pMol2 = theEnsemble2->getMoleculePointer(0);
    protein* bundle2 = static_cast<protein*>(pMol2);

    UInt chainNum = bundle->getNumChains();
    for (UInt i = 0; i < chainNum; i ++)
    {
        UInt resNum = bundle->getNumResidues(i);
        for (UInt j = 0; j < resNum; j ++)
        {
            double phi = bundle->getPhi(i,j);
            double psi = bundle->getPsi(i,j);
            if (j < 4 || j > 29)
            {

                if (j != 0)
                {
                    bundle2->setPhi(i, j, phi*-1);
                }
                if (j != resNum)
                {
                    bundle2->setPsi(i, j, psi*-1);
                }
            }
        }
    }
    bundle2->protOpt(false);
    pdbWriter(bundle2, outFile);
    return 0;
}

