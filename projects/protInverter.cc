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
#include <sstream>

enum aminoAcid {A,R,N,D,Dh,C,Cx,Cf,Q,E,Eh,Hd,He,Hn,Hp,I,L,K,M,F,P,O,S,T,W,Y,V,G,dA,dR,dN,dD,dDh,dC,dCx,dQ,dE,dEh,dHd,dHe,dHn,dHp,dI,dL,dK,dM,dF,dP,dO,dS,dT,dAT,dW,dY,dV,Hce,Pch,Csf};

int main (int argc, char* argv[])
{
    if (argc !=3)
    {   cout << "protInverter <inFile.pdb> <outFile.pdb>" << endl;
        exit(1); }

    string infile = argv[1];
    string outFile = argv[2];
    PDBInterface* thePDB = new PDBInterface(infile);
    ensemble* theEnsemble = thePDB->getEnsemblePointer();
    molecule* pMol = theEnsemble->getMoleculePointer(0);
    protein* bundle = static_cast<protein*>(pMol);
    //UIntVec rots = bundle->getAllowedRotamers(0,0,L,0);
    //UInt count = 0;
    //for (UInt h = 0; h < rots.size(); h++)
    //{
        UInt chainNum = bundle->getNumChains();
        for (UInt i = 0; i < chainNum; i ++)
        {
            UInt resNum = bundle->getNumResidues(i);
            for (UInt j = 0; j < resNum; j ++)
            {
                UInt restype = bundle->getTypeFromResNum(i,j);
                if (restype == L)
                {
                    bundle->mutateWBC(i,j,A);
                }
                if (restype == dL)
                {
                    bundle->mutateWBC(i,j,dA);
                }
            }
        }
       /* count++;
        stringstream convert;
        string countstr;
        convert << count, countstr = convert.str();
        string outFile = countstr + ".rot.pdb";
        pdbWriter(bundle, outFile);
    }*/
    //bundle->protOpt(false);
    pdbWriter(bundle, outFile);
    return 0;
}

