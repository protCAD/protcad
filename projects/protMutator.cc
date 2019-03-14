//*******************************************************************************************************
//*******************************************************************************************************
//**************************************                       ******************************************
//**************************************    protMutator 1.1    ******************************************
//**************************************                       ******************************************
//*******************************************************************************************************
//*******************************************************************************************************


//--Program setup----------------------------------------------------------------------------------------
#include <iostream>
#include <string>
#include "ensemble.h"
#include "PDBInterface.h"
#include <sstream>
#include <iterator>
#include <vector>


void randomizeSideChain(protein* _prot, UInt _chainIndex, UInt _resIndex);

int main (int argc, char* argv[])
{
	//--Program setup
    if (argc !=2)
	{
    cout << "protMutator <inFile.pdb>" << endl;
	exit(1);
	}
	enum aminoAcid {A,R,N,D,Dh,C,Cx,Cf,Q,E,Eh,Hd,He,Hp,I,L,K,M,F,P,O,S,T,W,Y,V,G,dA,dR,dN,dD,dDh,dC,dCx,dCf,dQ,dE,dEh,dHd,dHe,dHp,dI,dL,dK,dM,dF,dP,dO,dS,dT,dW,dY,dV,Csf,Sf4,Hca,Eoc,Oec};
	string infile = argv[1];
    vector<vector<UInt> > resIDs;
    vector <UInt> v;

    UInt resID1[] = {A,dCf,A,dA,Cf,dA,A,dCf,A,dA,Cf,dA,A,dCf,A,dA,Cf,dA,A,dCf,A,dA,Cf,dA,A,dCf,A,dA,Cf,dA,A,dCf,A,dA,Cf,dA};// -design1a
    v.insert (v.begin(), resID1, resID1 + sizeof(resID1)/sizeof(resID1[0]));
    resIDs.push_back(v);
    v.clear();

   
    //--Mutate chains
    for (UInt h = 0; h < resIDs.size(); h++)
    {
		PDBInterface* thePDB = new PDBInterface(infile);
		ensemble* theEnsemble = thePDB->getEnsemblePointer();
		molecule* pMol = theEnsemble->getMoleculePointer(0);
		protein* bundle = static_cast<protein*>(pMol);
		UInt chainNum = bundle->getNumChains();
		for (UInt i = 0; i < chainNum; i ++)
		{
			UInt resNum = bundle->getNumResidues(i);
			for (UInt j = 0; j < resNum; j++)
			{
				if (j >= resIDs[i].size())
				{
					bundle->removeResidue(i,j);
				}
				else
				{
					bundle->activateForRepacking(i, j);
					bundle->mutateWBC(i, j, resID1[j]);
					//randomizeSideChain(bundle, i, j);
				}
			}
		}
		stringstream convert;
		string countStr, outFile;
		convert << h+1, countStr = convert.str();
		outFile = countStr + ".mut.pdb";
		pdbWriter(bundle, outFile);
		delete thePDB;
    }
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
