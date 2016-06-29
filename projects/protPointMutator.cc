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
#include <unistd.h>
#include "ensemble.h"
#include "PDBInterface.h"

double getAverageDielectric(protein* _bundle, UInt _resIndex);
int main (int argc, char* argv[])
{
	//--Running parameters
	if (argc !=2)
	{
        cout << "protMutator <inFile.pdb>" << endl;
		exit(1);
	}
    string infile = argv[1];
    enum aminoAcid {A,R,N,D,Dh,C,Cx,Cf,Q,E,Eh,Hd,He,Hn,Hp,I,L,K,M,F,P,O,S,T,W,Y,V,G,dA,dR,dN,dD,dDh,dC,dCx,dQ,dE,dEh,dHd,dHe,dHn,dHp,dI,dL,dK,dM,dF,dP,dO,dS,dT,dAT,dW,dY,dV,Hce,Pch,Csf};
    string aminoAcidString[] = {"A","R","N","D","Dh","C","Cx","Cf","Q","E","Eh","Hd","He","Hn","Hp","I","L","K","M","F","P","O","S","T","W","Y","V","G","dA","dR","dN","dD","dDh","dC","dCx","dQ","dE","dEh","dHd","dHe","dHn","dHp","dI","dL","dK","dM","dF","dP","dO","dS","dT","dAT","dW","dY","dV","Hce","Pch"};
    PDBInterface* thePDB = new PDBInterface(infile);
    ensemble* theEnsemble = thePDB->getEnsemblePointer();
    molecule* pMol = theEnsemble->getMoleculePointer(0);
    protein* bundle = static_cast<protein*>(pMol);
    bundle->silenceMessages();
    srand (getpid());
    UInt numres = bundle->getNumResidues(0);
    delete thePDB;
	
    //--parameters
    //int chains[] = {0};
    //int chainsSize = sizeof(chains)/sizeof(chains[0]);
    //int residues[] = {1,3,5,7,9,13,15,17,19,21,28,30,32,34,36,54,57,59,61,70,72,74,76,78,86,88,90,92,94,108,110,112,114,116,127,129,131,133,151,153,155,157,159,161,166,169,172,174,181,183,185};//{46,47,61,92,95};//61,30,9,129/46,47,61,92,95/,87,91,110,150,152{1,3,5,7,9,13,15,17,19,21,28,30,32,34,36,54,57,59,61,70,72,74,76,78,86,88,90,92,94,108,110,112,114,116,127,129,131,133,151,153,155,157,159,161,166,169,172,174,181,183,185};
    //int residuesSize = sizeof(residues)/sizeof(residues[0]);
    int resID[] = {Csf};
    UInt resIDsize = sizeof(resID)/sizeof(resID[0]);
    double Energy, bestE;

	//--Mutations
    for (UInt h = 0; h < resIDsize; h++)
    {
        for (UInt i = 0; i < numres; i++)
        {
            PDBInterface* thePDB = new PDBInterface(infile);
            ensemble* theEnsemble = thePDB->getEnsemblePointer();
            molecule* pMol = theEnsemble->getMoleculePointer(0);
            protein* bundle = static_cast<protein*>(pMol);
            bestE = 1E100;
            bundle->activateForRepacking(0, i);
            bundle->mutate(0, i, resID[h]);
            cout << "test1" << endl;
            UIntVec allowedRots = bundle->getAllowedRotamers(0, i, resID[h], 0);
            for (UInt j = 0; j < allowedRots.size(); j++)
            {
                cout << "test2" << endl;
                bundle->setRotamerWBC(0, i, 0, allowedRots[j]);
                cout << "test3" << endl;
                Energy = bundle->intraSoluteEnergy(true);
                cout << "test4" << endl;
                if (Energy < bestE)
                {
                    bestE = Energy;
                    stringstream convert;
                    string countstr;
                    convert << i, countstr = convert.str();
                    cout << "test3.5" << endl;
                    string outFile = countstr + ".good.pdb";
                    cout << "test4.5" << endl;
                    pdbWriter(bundle, outFile);
                    cout << "test5" << endl;
                }
            }
            cout << i << " " << aminoAcidString[resID[h]] << " " << bestE << endl;
            delete thePDB;
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

