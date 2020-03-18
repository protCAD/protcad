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

enum aminoAcid {A,R,N,D,Dh,C,Cx,Cf,Q,E,Eh,Hd,He,Hp,I,L,K,M,F,P,O,S,T,W,Y,V,G,dA,dR,dN,dD,dDh,dC,dCx,dCf,dQ,dE,dEh,dHd,dHe,dHp,dI,dL,dK,dM,dF,dP,dO,dS,dT,dW,dY,dV,SF4,HEM,NI2,CLN,CO2,MG2,OH,OXY,CLD,HIS};

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
            for (UInt i = 0; i < currentRot.size(); i++)
            {
				for (UInt j = 0; j < currentRot[i].size(); j++)
				{
					currentRot[i][j] = currentRot[i][j]*-1;
				}
            }
            if (restype < G)
            {
                bundle->mutateWBC(i,j,restype+(G+1));
            }
            if (restype > G && restype < SF4)
            {
                bundle->mutateWBC(i,j,restype-(G+1));
            }
            bundle->setSidechainDihedralAngles(i, j, (currentRot));
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

