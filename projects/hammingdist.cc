//*******************************************************************************************************
//*******************************************************************************************************
//*****************************                         *************************************************
//*****************************        hammingdistance  *************************************************
//*****************************                         *************************************************
//*******************************************************************************************************
//*******************************************************************************************************

#include <iostream>
#include <string>
#include "ensemble.h"
#include "PDBInterface.h"

//--Program setup-------------------------------------------------------------
int main (int argc, char* argv[])
{
        if (argc !=3)
        {
                cout << "hammingdist <wildtype.pdb> <mutant.pdb>" << endl;
                exit(1);
        }
	enum aminoAcid {A,R,N,D,Dh,C,Cx,Cf,Q,E,Eh,Hd,He,Hp,I,L,K,M,F,P,O,S,T,W,Y,V,G,dA,dR,dN,dD,dDh,dC,dCx,dCf,dQ,dE,dEh,dHd,dHe,dHp,dI,dL,dK,dM,dF,dP,dO,dS,dT,dW,dY,dV,Sf4,Saf,Hem};
	string aminoAcidString[] = {"A","R","N","D","D","C","C","C","Q","E","E","H","H","H","I","L","K","M","F","P","O","S","T","W","Y","V","G","dA","dR","dN","dD","dDh","dC","dCx","dCf","dQ","dE","dEh","dHd","dHe","dHp","dI","dL","dK","dM","dF","dP","dO","dS","dT","dW","dY","dV","Csf","Sf4","Saf","Hem"};
        string infile1 = argv[1];
	string infile2 = argv[2];
        PDBInterface* thePDB1 = new PDBInterface(infile1);
        ensemble* theEnsemble1 = thePDB1->getEnsemblePointer();
        molecule* pMol1 = theEnsemble1->getMoleculePointer(0);
        protein* bundle1 = static_cast<protein*>(pMol1);
	PDBInterface* thePDB2 = new PDBInterface(infile2);
        ensemble* theEnsemble2 = thePDB2->getEnsemblePointer();
        molecule* pMol2 = theEnsemble2->getMoleculePointer(0);
        protein* bundle2 = static_cast<protein*>(pMol2);
	vector<string>str1,str2;
	//this loop will get the number of chains in the pdb and pulls out the residues from each chains for bundle1 (i.e., pdb1)
	for (UInt i=0; i < bundle1->getNumChains();i++){
		UInt numRes = bundle1->getNumResidues(i);
		for (UInt j = 0; j < numRes; j++)
		{
			UInt restype = bundle1->getTypeFromResNum(i,j);
			str1.push_back(aminoAcidString[restype]);
		}
	}
	//this gets the residues for second pdb to compare the hammingdistance of
	for (UInt i=0; i < bundle2->getNumChains();i++){
		UInt numRes = bundle2->getNumResidues(i);
		for (UInt j = 0; j < numRes; j++)
		{
			UInt restype = bundle2->getTypeFromResNum(i,j);
			 str2.push_back(aminoAcidString[restype]);
		}
	}
	cout << "Hamming distance of: "<<argv[1]<<" and "<< argv[2]<<" is: "<<bundle1->getHammingDistance(str1,str2)<<endl; //calculates the hammingdistance and prints the score
}
