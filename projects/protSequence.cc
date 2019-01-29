//*******************************************************************************************************
//*******************************************************************************************************
//*************************************                      ********************************************
//*************************************      getSequence     ********************************************
//*************************************                      ********************************************
//*******************************************************************************************************
//********************************* -get sequence from pdb file- ****************************************
//*******************************************************************************************************

//--Program setup----------------------------------------------------------------------------------------
#include <iostream>
#include <string>
#include <time.h>
#include <sstream>
#include "ensemble.h"
#include "PDBInterface.h"
int main (int argc, char* argv[])
{
	//--Running parameters
    if (argc !=2)
	{
        cout << "getSequence <inFile.pdb>" << endl;
		exit(1);
	}
	enum aminoAcid {A,R,N,D,Dh,C,Cx,Cf,Q,E,Eh,Hd,He,Hp,I,L,K,M,F,P,O,S,T,W,Y,V,G,dA,dR,dN,dD,dDh,dC,dCx,dCf,dQ,dE,dEh,dHd,dHe,dHp,dI,dL,dK,dM,dF,dP,dO,dS,dT,dW,dY,dV,Csf,Sf4,Hca,Eoc,Oec,Hem};
	string aminoAcidString[] = {"A","R","N","D","D","C","C","C","Q","E","E","H","H","H","I","L","K","M","F","P","O","S","T","W","Y","V","G","dA","dR","dN","dD","dD","dC","dC","dC","dQ","dE","dE","dH","dH","dH","dI","dL","dK","dM","dF","dP","dO","dS","dT","dW","dY","dV","Csf","Sf4","Hca","Eoc","Oec","Hem"};
	string backboneSeq[] =   {"", "M", "C", "L", "P", "B","E","Y","A","I","G",  "N",  "D",  "Q",  "R",  "F", "H", "W", "K", "S", "T"};
	string backboneTypes[] = {"","-γ","-π","-α","-ρ","-β","β","ρ","α","π","γ","-γl","-πl","-αl","-ρl","-βl","βl","ρl","αl","πl","γl"};
	string infile = argv[1];
	PDBInterface* thePDB = new PDBInterface(infile);
	ensemble* theEnsemble = thePDB->getEnsemblePointer();
	molecule* pMol = theEnsemble->getMoleculePointer(0);
	protein* _prot = static_cast<protein*>(pMol);

    UInt numChains = _prot->getNumChains();
    for (UInt i = 0; i < numChains; i++)
	{	cout << ">" << infile << ":" << _prot->getChainID(i) << endl;
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
		cout << ">" << infile << ":" << _prot->getChainID(i) << endl;
        UInt numRes = _prot->getNumResidues(i);
        for (UInt j = 0; j < numRes; j++)
        {
            UInt backboneType = _prot->getBackboneSequenceType(i,j);
			cout << backboneSeq[backboneType];
        }
        cout << endl;
    }
	cout << endl;
	pdbWriter(_prot,infile);
	return 0;
}



