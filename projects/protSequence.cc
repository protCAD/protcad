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
    if (argc !=3)
	{
		cout << "protSequence <includeBBseq(t/f)> <inFile.pdb>" << endl;
		exit(1);
	}
	enum aminoAcid {A,R,N,D,Dh,C,Cx,Cf,Q,E,Eh,Hd,He,Hp,I,L,K,M,F,P,O,S,T,W,Y,V,G,dA,dR,dN,dD,dDh,dC,dCx,dCf,dQ,dE,dEh,dHd,dHe,dHp,dI,dL,dK,dM,dF,dP,dO,dS,dT,dW,dY,dV,Sf4,Saf,Hem};
	string aminoAcidString[] = {"A","R","N","D","D","C","C","C","Q","E","E","H","H","H","I","L","K","M","F","P","O","S","T","W","Y","V","G","A","R","N","D","D","C","C","C","Q","E","E","H","H","H","I","L","K","M","F","P","O","S","T","W","Y","V","A","R","N","D","D","C","C","C","Q","E","E","H","H","H","I","L","K","M","F","P","O","S","T","W","Y","V","G","A","R","N","D","D","C","C","C","Q","E","E","H","H","H","I","L","K","M","F","P","O","S","T","W","Y","V","A","R","N","D","D","C","C","C","Q","E","E","H","H","H","I","L","K","M","F","P","O","S","T","W","Y","V","G","A","R","N","D","D","C","C","C","Q","E","E","H","H","H","I","L","K","M","F","P","O","S","T","W","Y","V"};
	string backboneSeq[] =   { "m", "c", "l", "p", "e","t","y","a","i","g",  "n",  "d",  "q",  "r",  "f", "h", "w", "k", "s", "v","-"};
	string backboneTypes[] = {"-γ","-π","-α","-ρ","-β","β","ρ","α","π","γ","-γi","-πi","-αi","-ρi","-βi","βi","ρi","αi","πi","γi","-"};
	string bbseq = argv[1];
	string infile = argv[2];
	PDBInterface* thePDB = new PDBInterface(infile);
	ensemble* theEnsemble = thePDB->getEnsemblePointer();
	molecule* pMol = theEnsemble->getMoleculePointer(0);
	protein* _prot = static_cast<protein*>(pMol);
	string pdb = infile.substr(3,infile.length()-7);
	bool BB = false;
	if (bbseq == "t"){BB = true;}

	fstream aa, bb;
	aa.open ("aa.faa", fstream::in | fstream::out | fstream::app);
	if (BB){bb.open ("bb.faa", fstream::in | fstream::out | fstream::app);}

	UInt numChains = _prot->getNumChains();
	for (UInt i = 0; i < numChains; i++)
	{
		bool write = true;
		bool firstbreak = true;
		UInt numRes = _prot->getNumResidues(i);
		for (UInt j = 0; j < numRes; j++)
		{
			if (!_prot->isCofactor(i,j)){
				if (write){
					aa << ">" << pdb << _prot->getChainID(i) << " " << _prot->getResNum(i,0) << endl;
					if (BB){bb << ">" << pdb << _prot->getChainID(i) << " " << _prot->getResNum(i,0) << endl;}
					write = false;
				}
				UInt restype = _prot->getTypeFromResNum(i,j);
				aa << aminoAcidString[restype];
				if (BB){
					UInt backboneType = _prot->getBackboneSequenceType(i,j);
					bb << backboneSeq[backboneType];
					if (backboneType == 20 && j > 0 && j < numRes-1 && !_prot->isCofactor(i,j+1)){
						if (firstbreak){
							bb << backboneSeq[backboneType];
							aa << backboneSeq[backboneType];
							firstbreak = false;
						} else {firstbreak = true;}
					}
				}
			}
		}
		aa << endl;
		if (BB){bb << endl;}
	}
	aa.close();
	if (BB){bb.close();}
	return 0;
}