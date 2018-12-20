//*******************************************************************************************************
//*******************************************************************************************************
//**************************************                       ******************************************
//**************************************   getDihedrals 1.0    ******************************************
//**************************************                       ******************************************
//*******************************************************************************************************
//*******************************************************************************************************


//--Program setup----------------------------------------------------------------------------------------
#include <iostream>
#include <string>
#include "ensemble.h"
#include "PDBInterface.h"
int main (int argc, char* argv[])
{
    if (argc !=2)
	{
    cout << "getDihedrals <inFile.pdb>" << endl;
	exit(1);
	}

	enum aminoAcid {A,R,N,D,Dh,C,Cx,Cf,Q,E,Eh,Hd,He,Hp,I,L,K,M,F,P,O,S,T,W,Y,V,G,dA,dR,dN,dD,dDh,dC,dCx,dCf,dQ,dE,dEh,dHd,dHe,dHp,dI,dL,dK,dM,dF,dP,dO,dS,dT,dW,dY,dV,Csf,Sf4,Hca,Eoc,Oec,Hem};
	string aminoAcidString[] = {"A","R","N","D","D","C","C","C","Q","E","E","H","H","H","I","L","K","M","F","P","O","S","T","W","Y","V","G","dA","dR","dN","dD","dD","dC","dC","dC","dQ","dE","dE","dH","dH","dH","dI","dL","dK","dM","dF","dP","dO","dS","dT","dW","dY","dV","Csf","Sf4","Hca","Eoc","Oec","Hem"};
	string infile = argv[1];
	PDBInterface* thePDB = new PDBInterface(infile);
	ensemble* theEnsemble = thePDB->getEnsemblePointer();
	molecule* pMol = theEnsemble->getMoleculePointer(0);
    protein* bundle = static_cast<protein*>(pMol);

	//cout << endl << "phi psi RPT" << endl;


	//--Search sequence of inFile for phis >= 0 --------------------------------------------------------
	//UInt chainNum = bundle->getNumChains();	
	//for (UInt i = 0; i < chainNum; i ++)
	//{
		UInt resNum = bundle->getNumResidues(0);
        //cout << "NA " << bundle->getPsi(i,0) << " NA" << endl;
        for (UInt j = 1; j < resNum-2; j ++)
        {
            UInt restype1 = bundle->getTypeFromResNum(0,j);
            UInt restype2 = bundle->getTypeFromResNum(0,j+1);
            cout << aminoAcidString[restype1] << aminoAcidString[restype2] << " " << bundle->getBackboneSequenceType(0,j) << bundle->getBackboneSequenceType(0,j+1) << endl;
        }
       // cout << bundle->getPhi(i,resNum-1) << " NA " << "NA" << endl;
	//}
	return 0;
}

