//*******************************************************************************************************
//*******************************************************************************************************
//**************************************                         ****************************************
//**************************************   registerEnergy 1.0    ****************************************
//**************************************                         ****************************************
//*******************************************************************************************************
//*******************************************************************************************************


//--Program setup----------------------------------------------------------------------------------------
#include <iostream>
#include <string>
#include "ensemble.h"
#include "PDBInterface.h"

void randomizeSideChain(protein* _prot, UInt _chainIndex, UInt _resIndex);

int main (int argc, char* argv[])
{
	//--Program setup
	if (argc !=3)
	{
	cout << "mutantMaker <inFile.pdb> <outFile.pdb>" << endl;
	exit(1);
	}
	enum aminoAcid {A,R,N,D,Dh,C,Cx,Q,E,Eh,Hd,He,Hn,Hp,I,L,K,M,F,P,O,S,T,W,Y,V,G,dA,dR,dN,dD,dDh,dC,dCx,dQ,dE,dEh,dHd,dHe,dHn,dHp,dI,dL,dK,dM,dF,dP,dO,dS,dT,dW,dY,dV};
	string infile = argv[1];
	PDBInterface* thePDB = new PDBInterface(infile);
	ensemble* theEnsemble = thePDB->getEnsemblePointer();
	molecule* pMol = theEnsemble->getMoleculePointer(0);
	protein* bundle = static_cast<protein*>(pMol);
	bundle->silenceMessages();
	residue::setCutoffDistance(9.0);
	pmf::setScaleFactor(0.0);
	rotamer::setScaleFactor(0.0);
	microEnvironment::setScaleFactor(0.0);
	amberVDW::setScaleFactor(1.0);
	amberVDW::setRadiusScaleFactor(1.0);
	amberVDW::setLinearRepulsionDampeningOff();
	amberElec::setScaleFactor(1.0);
	solvation::setItsScaleFactor(0.0);
	

	//--h3lix
	//UInt resID_A[] = {P,K,G,P,K,G,P,K,G,K,O,G,P,D,G,D,O,G,D,O,G,D,O,G,P,K,G,P,K,G};
	//UInt resID_B[] = {P,D,G,D,O,G,D,O,G,D,O,G,P,D,G,K,O,G,P,D,G,P,D,G,P,D,G,D,O,G};
	//UInt resID_C[] = {K,O,G,P,D,G,P,D,G,P,K,G,K,O,G,P,K,G,K,O,G,K,O,G,K,O,G,K,O,G};
	vector < vector < UInt > > proteinSequence;
	vector < UInt > chainSequence;
	


	for (UInt i = 0; i < chainNum; i ++)
	{
		resNum = bundle->getNumResidues(i);
		for (UInt j = 0; j < resNum; j++)
		{
			bundle->activateForRepacking(i, j);
			bundle->mutateWBC(i, j, resID_A[j]);
			randomizeSideChain(bundle, i, j);
		}
	}

	pdbWriter(bundle, outFile);

	return 0;
}


void randomizeSideChain(protein* _prot, UInt _chainIndex, UInt _resIndex)
{	
	UInt allowedRotsSize, randrot, restype;
	UIntVec allowedRots;
	restype = _prot->getTypeFromResNum(_chainIndex, _resIndex);
	allowedRots = _prot->getAllowedRotamers(_chainIndex, _resIndex, restype, 0);
	allowedRotsSize = allowedRots.size();
	if (allowedRotsSize > 2)
	{				
		randrot = rand() % allowedRotsSize;
		_prot->setRotamerWBC(_chainIndex, _resIndex, 0, allowedRots[randrot]);
	}
	return;
}

