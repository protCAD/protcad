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
    enum aminoAcid {A,R,N,D,Dh,C,Cx,Q,E,Eh,Hd,He,Hn,Hp,I,L,K,M,F,P,O,S,T,W,Y,V,G,dA,dR,dN,dD,dDh,dC,dCx,dQ,dE,dEh,dHd,dHe,dHn,dHp,dI,dL,dK,dM,dF,dP,dO,dS,dT,dW,dY,dV,Hce,Pch};
    string aminoAcidString[] = {"A","R","N","D","Dh","C","Cx","Q","E","Eh","Hd", "He","Hn","Hp","I","L","K","M","F","P","O","S","T","W","Y", "V","G","dA","dR","dN","dD","dDh","dC","dCx","dQ","dE","dEh","dHd","dHe","dHn","dHp","dI","dL","dK","dM","dF","dP","dO","dS","dT","dW","dY","dV","Hce","Pch"};
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
    int chains[] = {0};
	int chainsSize = sizeof(chains)/sizeof(chains[0]);
    int residues[] = {10};//61,30,9,129/46,47,61,92,95/,87,91,110,150,152{1,3,5,7,9,13,15,17,19,21,28,30,32,34,36,54,57,59,61,70,72,74,76,78,86,88,90,92,94,108,110,112,114,116,127,129,131,133,151,153,155,157,159,161,166,169,172,174,181,183,185};
    int residuesSize = sizeof(residues)/sizeof(residues[0]);
    int count = 0;

    //--variables
    UInt restype;
    double Energy;
    UIntVec allowedRots;
    cout << "pdb " << "energy " << endl;

    for (int i = 0; i < chainsSize; i++)
    {
        for (int j = 0; j < residuesSize; j++)
        {
            restype = bundle->getTypeFromResNum(chains[i], residues[j]);
            allowedRots = bundle->getAllowedRotamers(chains[i], residues[j], restype, 0);
            for (UInt m = 0; m < allowedRots.size(); m ++)
            {
               bundle->setRotamerWBC(chains[i], residues[j], 0, allowedRots[m]);
               Energy = bundle->protEnergy();
               count++;
               stringstream convert;
               string countstr;
               convert << count, countstr = convert.str();
               string outFile = countstr + ".pdb";
               pdbWriter(bundle, outFile);
               cout << count << " " << Energy << " " << endl;
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

