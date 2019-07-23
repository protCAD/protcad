//*******************************************************************************************************
//*******************************************************************************************************
//*************************************                      ********************************************
//*************************************      protSequence     ********************************************
//*************************************                      ********************************************
//*******************************************************************************************************
//********************************* -get sequence from pdb file- ****************************************
//*******************************************************************************************************

//--Program setup----------------------------------------------------------------------------------------
#include <iostream>
#include <string>
#include "ensemble.h"
#include "PDBInterface.h"
int main (int argc, char* argv[])
{
	//--Running parameters
    if (argc !=2)
	{
		cout << "protSequence <inFile.pdb>" << endl;
		exit(1);
	}
	enum aminoAcid {A,R,N,D,Dh,C,Cx,Cf,Q,E,Eh,Hd,He,Hp,I,L,K,M,F,P,O,S,T,W,Y,V,G,dA,dR,dN,dD,dDh,dC,dCx,dCf,dQ,dE,dEh,dHd,dHe,dHp,dI,dL,dK,dM,dF,dP,dO,dS,dT,dW,dY,dV,Sf4,Saf,Hem};
	string aminoAcidString[] = {"A","R","N","D","Dh","C","Cx","Cf","Q","E","Eh","Hd","He","Hp","I","L","K","M","F","P","O","S","T","W","Y","V","G","dA","dR","dN","dD","dDh","dC","dCx","dCf","dQ","dE","dEh","dHd","dHe","dHp","dI","dL","dK","dM","dF","dP","dO","dS","dT","dW","dY","dV","Csf","Sf4","Saf","Hem"};
	string backboneSeq[] =   { "C", "L", "P", "T","E","Y","A","I",  "D",  "Q",  "R",  "F", "H", "W", "K", "S"};
	string backboneTypes[] = {"-π","-α","-ρ","-β","β","ρ","α","π","-πi","-αi","-ρi","-βi","βi","ρi","αi","πi"};
	string infile = argv[1];
	PDBInterface* thePDB = new PDBInterface(infile);
	ensemble* theEnsemble = thePDB->getEnsemblePointer();
	molecule* pMol = theEnsemble->getMoleculePointer(0);
	protein* _prot = static_cast<protein*>(pMol);

	UInt numChains = _prot->getNumChains();
	for (UInt i = 0; i < numChains; i++)
	{	
		cout << ">" << infile << " amino acid sequence chain: "<< i << endl;
		UInt numRes = _prot->getNumResidues(i);
		for (UInt j = 0; j < numRes; j++)
		{
			UInt restype = _prot->getTypeFromResNum(i,j);
			cout << aminoAcidString[restype];
		}
		cout << endl;
	}
	for (UInt i = 0; i < numChains; i++)
	{
		cout << endl << ">" << infile << " backbone structure sequence chain:"<< i << endl;
		UInt numRes = _prot->getNumResidues(i);
		for (UInt j = 1; j < numRes-2; j++)
		{
			UInt backboneType = _prot->getBackboneSequenceType(i,j);
			cout << backboneSeq[backboneType];
		}
	}
	return 0;
}



