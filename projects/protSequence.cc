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
	string aminoAcidString[] = {"A","R","N","D","Dh","C","Cx","Cf","Q","E","Eh","Hd","He","Hp","I","L","K","M","F","P","O","S","T","W","Y","V","G","dA","dR","dN","dD","dDh","dC","dCx","dCf","dQ","dE","dEh","dHd","dHe","dHp","dI","dL","dK","dM","dF","dP","dO","dS","dT","dW","dY","dV","Csf","Sf4","Hca","Eoc","Oec","Hem"};
	string backboneTypes[] = {"Z","J","G","C","D","B","F","P","H","A","K","L","M","N","O","A","Q","R","S","T"};
	string infile = argv[1];
	PDBInterface* thePDB = new PDBInterface(infile);
	ensemble* theEnsemble = thePDB->getEnsemblePointer();
	molecule* pMol = theEnsemble->getMoleculePointer(0);
	protein* _prot = static_cast<protein*>(pMol);
	_prot->silenceMessages();
	residue::setCutoffDistance(9.0);
	rotamer::setScaleFactor(0.0);
	amberVDW::setScaleFactor(1.0);
	amberVDW::setRadiusScaleFactor(1.0);
	amberVDW::setLinearRepulsionDampeningOff();
	amberElec::setScaleFactor(1.0);

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
	cout << endl;
	
	for (UInt i = 0; i < numChains; i++)
	{
		cout << ">" << infile << ":" << _prot->getChainID(i) << endl;
        UInt numRes = _prot->getNumResidues(i);
        for (UInt j = 0; j < numRes; j++)
        {
            UInt backboneType = _prot->getBackboneSequenceType(i,j);
			cout << backboneTypes[backboneType];
        }
        cout << endl;
    }
	cout << endl;
	pdbWriter(_prot,infile);
	return 0;
}



