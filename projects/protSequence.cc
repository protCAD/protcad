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
	string aminoAcidString[] = {"A","R","N","D","D","C","C","C","Q","E","E","H","H","H","I","L","K","M","F","P","O","S","T","W","Y","V","G","dA","dR","dN","dD","dDh","dC","dCx","dCf","dQ","dE","dEh","dHd","dHe","dHp","dI","dL","dK","dM","dF","dP","dO","dS","dT","dW","dY","dV","Csf","Sf4","Saf","Hem"};
	string backboneSeq[] =   { "G", "C", "L", "P", "T","E","Y","A","I","M",  "N",  "D",  "Q",  "R",  "F", "H", "W", "K", "S", "V"};
	string backboneTypes[] = {"-γ","-π","-α","-ρ","-β","β","ρ","α","π","γ","-γi","-πi","-αi","-ρi","-βi","βi","ρi","αi","πi","γi"};
	string infile = argv[1];
	PDBInterface* thePDB = new PDBInterface(infile);
	ensemble* theEnsemble = thePDB->getEnsemblePointer();
	molecule* pMol = theEnsemble->getMoleculePointer(0);
	protein* _prot = static_cast<protein*>(pMol);
	
	//UInt numChains = _prot->getNumChains();
	//for (UInt i = 0; i < numChains; i++)
	//{	
	fstream aa;
	aa.open ("aa.faa", fstream::in | fstream::out | fstream::app);
	aa << endl << ">" << infile << " aaseq:" << endl;
	UInt numRes = _prot->getNumResidues(0);
	for (UInt j = 0; j < numRes; j++)
	{
		UInt restype = _prot->getTypeFromResNum(0,j);
		aa << aminoAcidString[restype];
	}
	aa.close();
	//}
	//for (UInt i = 0; i < numChains; i++)
	//{
	/*fstream bb;
	bb.open ("bb.faa", fstream::in | fstream::out | fstream::app);
	bb << endl << ">" << infile << " bbseq:" << endl;
	numRes = _prot->getNumResidues(0);
	for (UInt j = 0; j < numRes; j++)
	{
		UInt backboneType = _prot->getBackboneSequenceType(0,j);
		bb << backboneSeq[backboneType];
	}
	bb << endl;
	/*for (UInt j = 0; j < numRes; j++)
	{
		UInt backboneType = _prot->getBackboneSequenceType(0,j);;
		bb << backboneTypes[backboneType] << " ";
	}
	bb.close();*/
	//}
	return 0;
}



