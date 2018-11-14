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
    UInt restype, chainNum = bundle->getNumChains();
    vector < vector <double> > currentRot;
    for (UInt i = 0; i < chainNum; i ++)
    {
        UInt resNum = bundle->getNumResidues(i);
        for (UInt j = 0; j < resNum; j ++)
        {
            restype = bundle->getTypeFromResNum(i,j);
            currentRot = bundle->getSidechainDihedrals(i,j);
            if (restype < 27)
            {
                UInt restype = bundle->getTypeFromResNum(i,j);
                bundle->mutateWBC(i,j,restype+27);
            }
            if (restype > 27 && restype < 55)
            {
                UInt restype = bundle->getTypeFromResNum(i,j);
                bundle->mutateWBC(i,j,restype-27);
            }
            bundle->setSidechainDihedralAngles(i, j, currentRot);
            if (j == 0){
               double psi = bundle->getAngle(i,j,1);
               bundle->setDihedral(i,j,psi*-1,1,0);
            }
            else if (j == resNum-1){
                double phi = bundle->getAngle(i,j,0);
                bundle->setDihedral(i,j,phi*-1,0,0);
            }
            else{
                double phi = bundle->getAngle(i,j,0);
                double psi = bundle->getAngle(i,j,1);
                bundle->setDihedral(i,j,phi*-1,0,0);
                bundle->setDihedral(i,j,psi*-1,1,0);
            }
        }
    }
    pdbWriter(bundle, outFile);
    return 0;
}

