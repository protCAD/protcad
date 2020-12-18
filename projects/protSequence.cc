//*******************************************************************************************************
//*******************************************************************************************************
//*************************************                      ********************************************
//*************************************      protSequence     *******************************************
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
	string aminoAcidString[] = {"A","R","N","D","D","C","C","C","Q","E","E","H","H","H","I","L","K","M","F","P","O","S","T","W","Y","V","G","A","R","N","D","D","C","C","C","Q","E","E","H","H","H","I","L","K","M","F","P","O","S","T","W","Y","V","A","R","N","D","D","C","C","C","Q","E","E","H","H","H","I","L","K","M","F","P","O","S","T","W","Y","V","G","A","R","N","D","D","C","C","C","Q","E","E","H","H","H","I","L","K","M","F","P","O","S","T","W","Y","V","A","R","N","D","D","C","C","C","Q","E","E","H","H","H","I","L","K","M","F","P","O","S","T","W","Y","V","G","A","R","N","D","D","C","C","C","Q","E","E","H","H","H","I","L","K","M","F","P","O","S","T","W","Y","V"};
	string backboneSeq[] =   { "m", "c", "l", "p", "b","t","y","a","i","g",  "n",  "d",  "q",  "r",  "f", "h", "w", "k", "s", "v"};
	string backboneTypes[] = {"-γ","-π","-α","-ρ","-β","β","ρ","α","π","γ","-γi","-πi","-αi","-ρi","-βi","βi","ρi","αi","πi","γi"};
	string infile = argv[1];
	PDBInterface* thePDB = new PDBInterface(infile);
	ensemble* theEnsemble = thePDB->getEnsemblePointer();
	molecule* pMol = theEnsemble->getMoleculePointer(0);
	protein* _prot = static_cast<protein*>(pMol);
	
	fstream aa, bb;
	aa.open ("aa.faa", fstream::in | fstream::out | fstream::app);
	bb.open ("bb.faa", fstream::in | fstream::out | fstream::app);
	UInt numChains = _prot->getNumChains();
	for (UInt i = 0; i < numChains; i++)
	{
		aa << ">" << infile << "_" << _prot->getChainID(i) << endl;
		bb << ">" << infile << "_" << _prot->getChainID(i) << endl;
		UInt numRes = _prot->getNumResidues(i);
		for (UInt j = 0; j < numRes; j++)
		{
			if (!_prot->isCofactor(i,j)){
				UInt restype = _prot->getTypeFromResNum(i,j);
				aa << aminoAcidString[restype];
				UInt backboneType = _prot->getBackboneSequenceType(i,j);
				bb << backboneSeq[backboneType];
			}
		}
		aa << endl;
		bb << endl;
	}
	aa.close();
	bb.close();
	return 0;
}