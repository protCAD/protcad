//*******************************************************************************************************
//*******************************************************************************************************
//***********************************                    ************************************************
//***********************************  protMutator 1.5  *************************************************
//***********************************                    ************************************************
//*******************************************************************************************************
//**************   -point mutations, then backbone and sidechain optimization-   ************************
//*******************************************************************************************************

/////// Just specify a infile and preferred outfile name.

//--Program setup----------------------------------------------------------------------------------------
#include <iostream>
#include <string>
#include <time.h>
#include <sstream>
#include "ensemble.h"
#include "PDBInterface.h"

double getAverageDielectric(protein* _bundle, UInt _resIndex);
int main (int argc, char* argv[])
{
	//--Running parameters
	if (argc !=2)
	{
        cout << "protRotamer <inFile.pdb>" << endl;
		exit(1);
	}
    string infile = argv[1];
    enum aminoAcid {A,R,N,D,Dh,C,Cx,Cf,Q,E,Eh,Hd,He,Hn,Hp,I,L,K,M,F,P,O,S,T,W,Y,V,G,dA,dR,dN,dD,dDh,dC,dCx,dQ,dE,dEh,dHd,dHe,dHn,dHp,dI,dL,dK,dM,dF,dP,dO,dS,dT,dAT,dW,dY,dV,Hcd,Pch,Csf};
    //string aminoAcidString[] = {"A","R","N","D","Dh","C","Cx","Cf","Q","E","Eh","Hd","He","Hn","Hp","I","L","K","M","F","P","O","S","T","W","Y","V","G","dA","dR","dN","dD","dDh","dC","dCx","dQ","dE","dEh","dHd","dHe","dHn","dHp","dI","dL","dK","dM","dF","dP","dO","dS","dT","dAT","dW","dY","dV","Hcd","Pch","Csf"};
    PDBInterface* thePDB = new PDBInterface(infile);
    ensemble* theEnsemble = thePDB->getEnsemblePointer();
    molecule* pMol = theEnsemble->getMoleculePointer(0);
    protein* bundle = static_cast<protein*>(pMol);
    bundle->silenceMessages();
    residue::setCutoffDistance(9.0);
    rotamer::setScaleFactor(0.0);
    amberVDW::setScaleFactor(1.0);
    amberVDW::setRadiusScaleFactor(1.0);
    amberVDW::setLinearRepulsionDampeningOff();
    amberElec::setScaleFactor(1.0);
    srand (time(NULL));
	
    //--parameters
    int chains[] = {1};
    int chainsSize = sizeof(chains)/sizeof(chains[0]);
    int residues[] = {95};//61,30,9,129/46,47,61,92,95/,87,91,110,150,152{1,3,5,7,9,13,15,17,19,21,28,30,32,34,36,54,57,59,61,70,72,74,76,78,86,88,90,92,94,108,110,112,114,116,127,129,131,133,151,153,155,157,159,161,166,169,172,174,181,183,185};
    int residuesSize = sizeof(residues)/sizeof(residues[0]);
    //int count = 0;
    //--variables
    UInt restype;
    bool better;
    double bestE = bundle->intraSoluteEnergy(true), chi;
    string outFile = "opt." + infile;
    cout << "pdb " << "energy " << endl;
    UInt randres, randchain, randrot;
    UIntVec currentRot;
    UIntVec allowedRots = bundle->getAllowedRotamers(0,47,Hcd,0);


    for (UInt h = 0; h < 1000; h++)
    {
        randres = residues[rand() % residuesSize];
        randchain = chains[rand() % chainsSize];
        restype = bundle->getTypeFromResNum(randchain, randres);
        if (restype == Hcd)
        {
            better = false;
            currentRot = bundle->getCurrentRotamer(randchain, randres);
            randrot = allowedRots[rand() % allowedRots.size()];
            bundle->setRotamerWBC(randchain, randres, 0, randrot);
            for (UInt l = 0; l < 10; l++)
            {
                chi = l*36;
                bundle->setChi(randchain,randres,1,0,chi);
                double Energy = bundle->intraSoluteEnergy(true);
                if (Energy < bestE)
                {
                    bestE = Energy;
                    better = true;
                    pdbWriter(bundle, outFile);
                    cout << randchain << " " << randres << " " << Energy << " improved!" << endl;
                    break;
                }
            }
            if (!better)
            {
                bundle->setChi(randchain,randres,1,0,(chi*-1));
                bundle->setRotamerWBC(randchain, randres, 0, currentRot[0]);
            }
        }
    }
	return 0;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////
//////// functions //////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////

double getAverageDielectric(protein* _bundle, UInt _resIndex)
{
	UInt count = 0;
	double totalDielectric = 0.0, dielectric;
	for (UInt i = 0; i < _bundle->getNumChains(); i++)
	{
		for (UInt j = 0; j < _bundle->getNumAtoms(i, _resIndex); j++)
		{
			dielectric = _bundle->getDielectric(i,_resIndex,j);
			totalDielectric += dielectric;
			count++;
		}
	}
	return totalDielectric/count;
} 

