//*******************************************************************************************************
//*******************************************************************************************************
//**************************************                       ******************************************
//**************************************    mutantMaker 1.1    ******************************************
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
    cout << "mutantMaker <inFile.pdb>" << endl;
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
	//UInt resID_C[] = {dP,K,G,R,E,G,R,R,G,K,D,G,K,E,G,E,E,G,D,O,G,P,E,G,E,E,G,D,K,G};
	//UInt resID_A[] = {dP,O,G,D,R,G,D,O,G,E,E,G,P,O,G,D,D,G,R,E,G,K,D,G,R,R,G,P,K,G};
	//UInt resID_B[] = {dP,E,G,P,E,G,E,E,G,P,E,G,P,K,G,R,D,G,P,K,G,E,E,G,D,R,G,D,E,G};

	//UInt resID_A[] = {P,O,G,P,O,G,P,O,G,G,O,G,P,O,G,P,O,G,P,O,G,P,O,G};

	//UInt resID_A[] = {G,G,G,G,G,G,G,G,G,G,G,G,G,G,G,G,G,G,G,G,G,G,G,G,G,G,G,G,G,G,G,G,G,G,G,G,G,G,G,G,G,G,G,G,G,G,G,G,G,G,G,G,G,G,G};
	//UInt resID_A[] = {P,K,G,P,K,G,P,K,G,P,K,G,P,K,G,P,K,G,P,K,G,P,K,G,P,K,G,P,K,G};
	//UInt resID_B[] = {D,K,G,D,K,G,D,K,G,D,K,G,D,K,G,D,K,G,D,K,G,D,K,G,D,K,G,D,K,G};
	//UInt resID_C[] = {E,P,G,E,P,G,E,P,G,E,P,G,E,P,G,E,P,G,E,P,G,E,P,G,E,P,G,E,P,G};
	//UInt resID_A[] = {P,K,G,P,K,G,D,O,G,P,O,G,D,K,G,D,K,G,P,K,G,P,O,G,D,K,G,P,O,G};
	//UInt resID_B[] = {P,O,G,D,O,G,D,K,G,P,O,G,P,O,G,D,K,G,D,O,G,D,K,G,P,K,G,D,O,G};
	//UInt resID_C[] = {P,K,G,P,O,G,P,K,G,D,K,G,P,O,G,P,O,G,D,K,G,P,O,G,D,O,G,D,O,G};
	//UInt resID_A[] = {P,O,G,P,O,G,P,O,G,P,O,G,P,I,G,L,I,G,P,O,G,P,O,G,P,O,G,P,O,G};
	//UInt resID_A[] = {P,O,G,P,O,G,P,O,G,P,O,G,L,I,G,L,I,G,P,O,G,P,O,G,P,O,G,P,O,G};
	//UInt resID_A[] = {P,O,G,P,O,G,P,O,G,P,I,G,L,I,G,L,I,G,P,O,G,P,O,G,P,O,G,P,O,G};
	//UInt resID_A[] = {P,O,G,P,O,G,P,O,G,L,I,G,L,I,G,L,I,G,P,O,G,P,O,G,P,O,G,P,O,G};
    //UInt resID_A[] = {Q,R,L,R,L,R,L,E,N,V,G,S,N,K,G,A,R,L,R,L,R,L,G,G,V,V};
	//UInt resID_A[] = {K,K,G,K,K,G,K,K,G,P,O,G,P,O,G,P,O,G,P,O,G,P,O,G,P,O,G};
	//UInt resID_B[] = {D,D,G,D,D,G,D,D,G,P,O,G,P,O,G,P,O,G,P,O,G,P,O,G,P,O,G};

	//UInt resID_A[] = {P,I,G,P,P,G,P,R,G,N,R,G,E,R,G,S,E,G,S,P,G,He,P,G,M,P,G,P,P,G,P,P,G,A,P,G,P,C,C,G,G};
	//--triad
	//UInt resID_A[] = {S,M,E,S,L,S,K,T,He,He,Y,R,G,P,P,G,P,P,G,P,P,G,P,P,G,P,P,G,P,P,G,P,P,G,P,P,G,P,P,G,P,P};
	//UInt resID_A[] = {G,G,G,G,G,G,G,G,G,G,G,G,G,P,P,G,P,P,G,P,P,G,P,P,G,P,P,G,P,P,G,P,P,G,P,P,G,P,P,G,P,P};
	//UInt resID_A[] = {G,G,G,G,G,G,G,G,G,G,G,G,G,G,G,G,G,G,G,G};
	
	//--helix
	//UInt resID_A[] = {dA,dA,dA,dA,dA,dA,dA,dA,dA,dA,dA,dA,dA,dA,dA,G,G,A,A,A,A,A,A,A,A,A,A,A,A,A,A,A};
	//UInt resID_A[] = {V,V,V,V,V,V,V,V,V,V,V,V,V,V,V,V,V,V,V,V,V,V,V,V,V,V,V,V,V,V,V,V,V,V,V,V,V,V,V,V,V,V,V,V,V,V,V,V,V,V,V,V,V,V,V,V,V,V};
	//--ldhelix
	//UInt resID_A[] = {dA,G,G,G,G,G,dE,G,G,dQ,G,dHe,G,dE,G,D,G,D,G,G,W,G,G,G,G,G,G,G,G,G,E,D}; //selection2
	//UInt resID_A[] = {dA,dA,dA,dA,dA,dA,dA,dA,dA,dA,dA,dA,dA,dA,dA,G,G,A,A,A,A,A,A,A,A,A,A,A,A,A,A,A}; 

	//--CoiledCoil
	//UInt resID_A[] = {dE,V,Q,E,L,E,Q,K,V,Q,E,L,E,Q,K,V,Q,E,L,E,Q,K,V,Q,E,L,E,Q,K,V,dK};
    //UInt resID_A[] = {R,M,K,Q,L,E,D,K,N,E,E,L,L,S,K,N,Y,He,L,E,N,E,N,A,R,L,K,KEnergy.clear();,L,I,G};

	//UInt resID_A[] = {M,V,S,K,G,E,E,D};
	//UInt resID_A[] = {K,Y,L,E,D,M,G,G};

	//--tubeKnot
	//UInt resID_A[] = {S,dF,dF,dF,F,dV,V,dV,V,dV,V,dV,V,dV,V,dV,V,dV,V,dV,V,dV,V,dV,V,dV,V,dV,V,dV,V,dV,V,dV,V,dV,F,dS};
	//UInt resID_A[] = {A,dA,A,dA,A,dA,A,dA,A,dA,A,dA,A,dA,A};

	//--doubleTube
	//UInt resID_A[] = {A,dA,A,dA,A,dA,A,dA,A,dA,A,dA,A,dV,V,dV,V,dV,V,dV,V,dV,V,dV,V,dV,V,dV,V,dV,V,dV,V,dV,V,dV,V,dV,V,dV,V,dV

	//--Beta meander
	//UInt resID_A[] = {E,Q,Q,V,E,L,Q,dE,K,E,V,E,L,Q,V,K,L,E,E,dK,Q,E,Q,Q,V,E,L,Q,E,K};

	//--ns1rbd
	//UInt resID_A[] = {dD,A,A,A,A,A,A,A,A,A,A,dQ,G};

	//--effector domain
	//UInt resID_A[] = {N,Y,F,Y,S,L,F,Q,G};

	//--soybean
	//UInt resID_A[] = {N,E,C,Q,I,Q,K,L,N,A,L,K,P,D,N,R,I,E,S,E,G,G,F,I,E,T,W,N,P,N,N,K,P,F,Q,C,A,G,V,A,L,S,R,C,T,L,N,R,N,A,L,R,R,P,S,Y,T,N,G,P,Q,E,I,Y,I,Q,Q,G,N,G,I,F,G,M,I,F,P,G,C,P,S,T,Y,Q,E,P,Q,E,S,He,Q,K,V,He,R,F,R,E,G,D,L,I,A,V,P,T,G,V,A,W,W,M,Y,N,N,E,D,T,P,V,V,A,V,S,I,I,D,T,N,S,L,E,N,Q,L,D,Q,M,P,R,R,F,Y,L,A,G,N,Q,E,Q,E,Q,Q,E,E,E,N,E,G,S,N,I,L,S,G,F,A,P,E,F,L,K,E,A,F,G,V,N,M,Q,I,V,R,N,L,Q,G,E,N,E,E,E,D,S,I,V,T,V,K,G,G,L,R,V,T,A,P,A,M,G,I,D,E,T,I,C,T,M,R,L,R,Q,N,I,G,Q,N,S,S,P,D,I,Y,N,P,Q,A,G,S,I,T,T,A,T,S,L,D,F,P,A,L,W,L,L,K,L,S,A,Q,Y,G,S,L,R,K,N,A,M,F,V,P,He,Y,T,L,N,A,N,S,I,I,Y,A,L,N,G,R,A,L,V,Q,V,V,N,C,N,G,E,R,V,F,D,G,E,L,Q,E,G,G,V,L,I,V,P,Q,N,F,A,V,A,A,K,S,Q,S,D,N,F,E,Y,V,S,F,K,T,N,D,R,P,S,I,G,N,L,N,S,L,L,N,A,L,P,E,E,V,I,Q,He,T,F,N,L,K,S,Q,Q,A,R,Q,V,K,N,N,N,P,F,S,F,L,V,P,P};

	//--lenc1
	//UInt resID_A[] = {S,N,P,F,I,F,K,S,N,R,F,Q,T,I,Y,E,N,E,N,G,He,I,R,L,L,Q,R,F,D,K,R,S,K,I,F,E,N,L,Q,N,Y,R,L,L,E,Y,K,S,K,P,He,T,I,F,L,P,Q,F,T,D,A,D,F,I,L,V,V,L,S,G,K,A,I,L,T,V,L,N,S,N,D,R,N,S,F,N,L,E,R,G,D,T,I,K,L,P,A,G,T,I,A,Y,L,A,N,R,D,D,N,E,D,L,R,V,L,D,L,A,I,P,V,N,R,P,G,Q,L,Q,S,F,L,L,S,G,T,Q,N,Q,P,S,F,L,S,G,F,S,K,N,I,L,E,A,A,F,N,T,E,Y,E,E,I,E,K,V,L,L,E,E,Q,E,E,D,V,I,V,K,V,S,R,E,Q,E,P,F,N,L,R,S,R,N,P,I,Y,S,N,K,F,G,K,F,F,E,I,T,P,E,K,N,P,Q,L,Q,D,L,D,I,F,V,N,S,V,E,I,K,E,G,S,L,L,L,P,N,Y,N,S,R,A,I,V,I,V,T,V,N,E,G,K,G,D,F,E,L,V,G,Q,Q,V,Q,R,Y,R,A,R,L,S,P,G,D,V,L,V,I,P,A,G,He,P,V,A,I,N,A,S,S,D,L,N,L,I,G,F,G,I,N,A,K,N,N,Q,R,N,F,L,A,G,E,E,D,N,V,I,S,Q,I,Q,R,P,V,K,E,L,A,F,P,G,S,S,R,E,V,D,R,L,L,T,N,Q,K,Q,S,He,F,A,N,A,Q};
	//--5ht2a
	//UInt resID_A[] = {W,S,A,L,L,T,A,V,V,I,I,L,T,I,A,G,N,I,L,V,I,M,A,V,S,L,E,K,K,L,Q,N,A,T,N,Y,F,L,M,S,L,A,I,A,D,M,L,L,G,F,L,V,M,P,V,S,M,L,T,I,L,Y,G,Y,R,W,P,L,P,S,K,L,C,A,V,W,I,Y,L,D,V,L,F,S,T,A,S,I,M,He,L,C,A,I,S,L,D,R,Y,V,A,I,Q,N,P,I,He,He,S,R,F,N,S,R,T,K,A,F,L,K,I,I,A,V,W,T,I,S,V,G,I,S,M,P,I,P,V,F,G,L,Q,D,D,S,K,V,F,S,C,L,L,A,D,D,R,F,G,N,F,V,L,I,G,S,F,V,S,F,F,I,P,L,T,I,M,V,I,T,Y,F,L,T,I,K,S,L,Q,K,A,C,K,V,L,G,I,V,F,F,L,F,V,V,M,W,C,P,F,F,I,T,N,I,M,A,V,I,C,E,S,C,N,E,D,V,I,G,A,L,L,N,V,F,V,W,I,G,Y,L,S,S,A,V,N,P,L,V,Y,T,L,F,N,K,T,Y,R,S,A,F,S,R,Y,I,Q,C,Q,Y,K};	
	

	//--3bdc
	//,L,He,K,E,P,A,T,L,I,K,A,I,D,G,D,T,V,K,L,M,Y,K,G,Q,P,M,T,F,R,L,L,L,V,D,T,P,E,F,N,E,K,Y,G,P,E,A,S,A,F,T,K,K,M,V,E,N,A,K,K,I,E,V,E,F,D,K,G,Q,R,T,D,K,Y,G,R,G,L,A,Y,I,Y,A,D,G,K,M,V,N,E,A,L,V,R,Q,G,L,A,K,V,A,Y,V,Y,K,G,N,N,T,He,E,Q,L,L,R,K,A,E,A,Q,A,K,K,E,K,L,N,I,W,S
	//UInt resID_A[] = {L,He,K,E,P,A,T,L,I,K,A,I,D,G,D,T,V,K,L,M,Y,K,G,Q,P,M,T,F,R,L,E,L,V,D,T,P,E,F,N,E,K,Y,G,P,E,A,S,A,F,T,K,K,M,V,E,N,A,K,K,I,E,V,E,F,D,K,G,Q,R,T,D,K,Y,G,R,G,L,A,Y,I,Y,A,D,G,K,M,V,N,E,A,L,V,R,Q,G,L,A,K,V,A,Y,V,Y,K,G,N,N,T,He,E,Q,L,L,R,K,A,E,A,Q,A,K,K,E,K,L,N,I,W,S};

	//--trpcage triad
	//UInt resID_A[] = {N,D,He,I,Q,W,L,K,D,dQ,G,P,S,S,G,R,P,P,S,S};

     //UInt resID1[] = {Q,K,L,V,Cx,F,A,E,N,V,G,S,N,K,G,A,I,Cx,G,L,M,V,G,G,V,V};

    vector<vector<UInt> > resIDs;
    vector <UInt> v;

    UInt resID1[] = {G,Y,V,Q,E,R,T,I,F,F,K,D,D,G,N,Y,K,T,R,A,E,He,K,F,E,G,D,T,L,V,N,R,I,E,L,K,G,I,D,F,K,E,D,G,N,I,L,G,Q,K,L,E,Y,N,Y,N,S,V,N,V,Y,I,M,A,D,K,Q,K,N,G,I,K,V,N,F,K,I,R,V,N,I,E,D,G,S,V,Q,L,A,D,F,Y,Q,Q,N,T,P,I,G,D,G,P,V,L,L,P,D,N,Q,Y,L,S,T,Q,S,A,He,S,K,D,P,N,E,K,R,D,N,M,V,L,L,E,F,V,T,A,A,G,I,T,A,S,K,G,E,E,L,F,T,G,V,V,P,I,L,V,E,L,D,G,D,V,N,G,Q,K,F,S,V,S,G,E,G,E,G,D,A,T,Y,G,K,L,T,L,K,F,I,S,T,T};
    v.insert (v.begin(), resID1, resID1 + sizeof(resID1)/sizeof(resID1[0]));
    resIDs.push_back(v);
    v.clear();
    UInt resID2[] = {G,Y,V,Q,E,V,T,I,F,F,K,D,D,G,N,Y,K,T,R,A,E,He,K,F,E,G,D,T,L,V,N,R,I,E,L,K,G,I,D,F,K,E,D,G,N,I,L,G,Q,K,L,E,Y,N,Y,N,S,V,N,V,Y,I,M,A,D,K,Q,K,N,G,I,K,V,N,F,K,I,R,V,N,I,E,D,G,S,V,Q,L,A,D,F,Y,Q,Q,N,T,P,I,G,D,G,P,V,L,L,P,D,N,Q,Y,L,S,T,Q,S,A,He,S,K,D,P,N,E,K,R,D,N,M,V,L,L,A,F,V,T,A,A,G,I,T,A,S,K,G,E,E,L,F,T,G,V,V,P,I,L,V,E,L,D,G,D,V,N,G,Q,K,F,S,V,S,G,E,G,E,G,D,A,T,Y,G,K,L,T,L,K,F,I,S,T,T};
    v.insert (v.begin(), resID2, resID2 + sizeof(resID2)/sizeof(resID2[0]));
    resIDs.push_back(v);
    v.clear();
    UInt resID3[] = {G,Y,V,L,E,V,T,I,F,F,K,D,D,G,N,Y,K,I,R,A,E,He,K,F,E,G,D,T,L,V,A,R,I,E,L,K,G,I,D,F,K,E,D,G,N,I,L,G,Q,K,L,E,Y,N,Y,N,S,V,N,V,Y,I,M,A,D,K,Q,K,N,G,I,K,V,N,F,K,I,R,V,N,I,E,D,G,S,V,Q,L,A,D,F,Y,L,Q,A,T,P,I,G,D,G,P,V,L,L,P,D,N,Q,Y,L,S,I,Q,A,A,He,S,K,D,P,N,E,K,R,D,N,M,V,L,L,A,F,V,T,A,A,G,I,T,A,S,K,G,E,E,L,F,T,G,V,V,P,I,L,V,E,L,D,G,D,V,N,G,Q,K,F,S,V,S,G,E,G,E,G,D,A,T,Y,G,K,L,T,L,K,F,I,S,T,T};
    v.insert (v.begin(), resID3, resID3 + sizeof(resID3)/sizeof(resID3[0]));
    resIDs.push_back(v);
    v.clear();
    UInt resID4[] = {G,Y,V,L,E,V,T,I,F,F,K,D,D,G,N,Y,K,I,R,A,E,He,K,F,E,G,D,T,F,V,A,R,I,E,L,K,G,I,D,F,K,E,D,G,N,I,L,G,Q,K,L,E,Y,N,Y,N,S,V,N,L,Y,I,M,A,D,K,Q,K,N,G,I,K,V,N,F,K,I,R,V,N,I,E,D,G,S,V,Q,L,A,D,F,Y,L,Q,A,T,P,I,G,D,G,P,V,L,L,P,D,N,Q,Y,L,S,I,Q,A,A,He,S,K,D,P,N,E,K,R,D,N,M,V,L,L,A,F,L,T,A,A,G,I,T,A,S,K,G,E,E,L,F,T,G,V,V,P,I,L,V,E,L,D,G,D,V,N,G,Q,K,F,S,V,S,G,E,G,E,G,D,A,T,Y,G,K,L,T,L,K,F,I,Y,T,T};
    v.insert (v.begin(), resID4, resID4 + sizeof(resID4)/sizeof(resID4[0]));
    resIDs.push_back(v);
    v.clear();
    UInt resID5[] = {G,Y,V,L,E,V,T,I,F,F,K,D,D,G,N,Y,K,I,R,A,E,He,K,F,E,G,D,T,F,V,A,R,I,E,L,K,G,I,D,F,K,E,D,G,N,I,L,G,Q,K,L,E,Y,N,Y,N,S,V,N,L,Y,I,M,A,D,K,Q,K,N,G,I,K,V,N,F,K,I,R,V,N,I,E,D,G,S,V,Q,L,A,D,F,Y,L,Q,A,T,P,I,G,D,G,P,V,L,L,P,D,N,Q,Y,L,S,I,Q,A,A,He,S,K,D,P,N,E,K,R,D,N,M,V,L,L,A,F,L,T,A,A,G,I,T,A,S,K,G,E,E,L,F,T,G,V,V,P,I,L,V,E,L,D,G,D,V,N,G,Q,K,F,S,V,S,G,E,G,E,G,D,A,T,Y,G,K,L,T,L,K,F,I,Y,T,T};
    v.insert (v.begin(), resID5, resID5 + sizeof(resID5)/sizeof(resID5[0]));
    resIDs.push_back(v);
    v.clear();
    UInt resID6[] = {G,Y,V,L,E,V,T,I,F,F,K,D,D,G,N,F,K,I,R,A,E,He,K,F,E,G,D,T,F,V,A,R,I,E,L,K,G,I,D,F,K,E,D,G,N,I,L,G,Q,K,Y,E,Y,N,Y,N,S,V,N,L,Y,I,M,A,D,K,Q,K,N,G,I,K,V,N,F,K,I,R,V,N,I,E,D,G,S,V,Q,L,A,D,F,Y,L,Q,A,T,P,I,G,D,G,P,V,L,L,P,D,N,Q,Y,L,S,I,Q,A,A,He,S,K,D,P,N,E,K,R,D,N,M,V,L,L,A,F,L,T,A,A,G,I,T,A,S,K,G,E,E,L,F,T,G,V,V,P,I,L,V,E,L,D,G,D,V,N,G,Q,K,F,S,V,S,G,E,G,E,G,D,A,T,Y,G,K,L,T,L,K,F,I,Y,T,T};
    v.insert (v.begin(), resID6, resID6 + sizeof(resID6)/sizeof(resID6[0]));
    resIDs.push_back(v);
    v.clear();
    /*UInt resID7[] = {G,F,O,G,A,A,G,R,T,G,P,O,G,P,S,G,I,S,G,P,O,G,P,O,G,P,O,G,P,O,G,V};
    v.insert (v.begin(), resID7, resID7 + sizeof(resID7)/sizeof(resID7[0]));
    resIDs.push_back(v);
    v.clear();
    UInt resID8[] = {G,S,O,G,A,Q,G,L,Q,G,P,R,G,L,O,G,T,O,G,P,O,G,P,O,G,P,O,G,P,O,G,V};
    v.insert (v.begin(), resID8, resID8 + sizeof(resID8)/sizeof(resID8[0]));
    resIDs.push_back(v);
    v.clear();
    UInt resID9[] = {G,G,K,G,D,A,G,A,O,G,E,R,G,P,O,G,L,A,G,P,O,G,P,O,G,P,O,G,P,O,G,V};
    v.insert (v.begin(), resID9, resID9 + sizeof(resID9)/sizeof(resID9[0]));
    resIDs.push_back(v);
    v.clear();
    UInt resID10[] = {G,D,A,G,A,O,G,E,R,G,P,O,G,L,A,G,A,O,G,P,O,G,P,O,G,P,O,G,P,O,G,V};
    v.insert (v.begin(), resID10, resID10 + sizeof(resID10)/sizeof(resID10[0]));
    resIDs.push_back(v);
    v.clear();
    UInt resID11[] = {G,A,O,G,E,R,G,P,O,G,L,A,G,A,O,G,L,R,G,P,O,G,P,O,G,P,O,G,P,O,G,V};
    v.insert (v.begin(), resID11, resID11 + sizeof(resID11)/sizeof(resID11[0]));
    resIDs.push_back(v);
    v.clear();
    UInt resID12[] = {G,E,R,G,P,O,G,L,A,G,A,O,G,L,R,G,G,A,G,P,O,G,P,O,G,P,O,G,P,O,G,V};
    v.insert (v.begin(), resID12, resID12 + sizeof(resID12)/sizeof(resID12[0]));
    resIDs.push_back(v);
    v.clear();
    UInt resID13[] = {G,P,O,G,L,A,G,A,O,G,L,R,G,G,A,G,P,O,G,P,O,G,P,O,G,P,O,G,P,O,G,V};
    v.insert (v.begin(), resID13, resID13 + sizeof(resID13)/sizeof(resID13[0]));
    resIDs.push_back(v);
    v.clear();
    UInt resID14[] = {G,L,A,G,A,O,G,L,R,G,G,A,G,P,O,G,P,E,G,P,O,G,P,O,G,P,O,G,P,O,G,V};
    v.insert (v.begin(), resID14, resID14 + sizeof(resID14)/sizeof(resID14[0]));
    resIDs.push_back(v);
    v.clear();
    UInt resID36[] = {G,P,O,G,A,O,G,P,L,G,I,A,G,I,T,G,A,R,G,L,A,G,P,O,G,P,O,G,P,O,G,P,O,G,G};
    v.insert (v.begin(), resID36, resID36 + sizeof(resID36)/sizeof(resID36[0]));
    resIDs.push_back(v);
    v.clear();
    UInt resID15[] = {G,P,R,G,N,R,G,E,R,G,S,E,G,S,O,G,He,O,G,Q,O,G,P,O,G,P,O,G,A,O,G,V};
    v.insert (v.begin(), resID15, resID15 + sizeof(resID15)/sizeof(resID15[0]));
    resIDs.push_back(v);
    v.clear();
    UInt resID16[] = {G,L,A,G,E,O,G,K,O,G,I,O,G,L,O,G,R,A,G,P,O,G,P,O,G,P,O,G,P,O,G,V};
    v.insert (v.begin(), resID16, resID16 + sizeof(resID16)/sizeof(resID16[0]));
    resIDs.push_back(v);
    v.clear();
    UInt resID17[] = {G,E,R,G,E,R,G,E,K,G,E,R,G,E,Q,G,R,D,G,P,O,G,P,O,G,P,O,G,P,O,G,V};
    v.insert (v.begin(), resID17, resID17 + sizeof(resID17)/sizeof(resID17[0]));
    resIDs.push_back(v);
    v.clear();
    UInt resID18[] = {G,P,O,G,P,O,G,P,O,G,K,R,G,K,O,G,P,O,G,P,O,G,P,O,G,G};
    v.insert (v.begin(), resID18, resID18 + sizeof(resID18)/sizeof(resID18[0]));
    resIDs.push_back(v);
    v.clear();
    UInt resID19[] = {G,P,O,G,P,O,G,R,O,G,K,R,G,K,O,G,P,O,G,P,O,G,P,O,G,G};
    v.insert (v.begin(), resID19, resID19 + sizeof(resID19)/sizeof(resID19[0]));
    resIDs.push_back(v);
    v.clear();
    UInt resID20[] = {G,P,O,G,P,O,G,R,O,G,K,R,G,K,Q,G,P,O,G,P,O,G,P,O,G,G};
    v.insert (v.begin(), resID20, resID20 + sizeof(resID20)/sizeof(resID20[0]));
    resIDs.push_back(v);
    v.clear();
    UInt resID21[] = {G,P,O,G,P,O,G,R,Q,G,K,R,G,K,O,G,P,O,G,P,O,G,P,O,G,G};
    v.insert (v.begin(), resID21, resID21 + sizeof(resID21)/sizeof(resID21[0]));
    resIDs.push_back(v);
    v.clear();
    UInt resID22[] = {G,P,O,G,P,O,G,R,O,G,K,R,G,K,Q,G,Q,K,G,P,O,G,P,O,G,G};
    v.insert (v.begin(), resID22, resID22 + sizeof(resID22)/sizeof(resID22[0]));
    resIDs.push_back(v);
    v.clear();
    UInt resID23[] = {G,P,O,G,P,O,G,P,O,G,R,K,G,R,O,G,P,O,G,P,O,G,P,O,G,G};
    v.insert (v.begin(), resID23, resID23 + sizeof(resID23)/sizeof(resID23[0]));
    resIDs.push_back(v);
    v.clear();
    UInt resID24[] = {G,P,O,G,P,O,G,R,O,G,R,K,G,R,O,G,P,O,G,P,O,G,P,O,G,G};
    v.insert (v.begin(), resID24, resID24 + sizeof(resID24)/sizeof(resID24[0]));
    resIDs.push_back(v);
    v.clear();
    UInt resID25[] = {G,P,O,G,R,O,G,R,K,G,R,O,G,V,R,G,P,R,G,P,O,G,P,O,G,P,O,G,P,O,G};
    v.insert (v.begin(), resID25, resID25 + sizeof(resID25)/sizeof(resID25[0]));
    resIDs.push_back(v);
    v.clear();
    UInt resID26[] = {G,P,O,G,R,O,G,A,A,G,A,O,G,V,R,G,P,R,G,P,O,G,P,O,G,P,O,G,P,O,G};
    v.insert (v.begin(), resID26, resID26 + sizeof(resID26)/sizeof(resID26[0]));
    resIDs.push_back(v);
    v.clear();
    UInt resID27[] = {G,P,O,G,A,O,G,R,K,G,R,O,G,V,A,G,P,A,G,P,O,G,P,O,G,P,O,G,P,O,G};
    v.insert (v.begin(), resID27, resID27 + sizeof(resID27)/sizeof(resID27[0]));
    resIDs.push_back(v);
    v.clear();
    UInt resID28[] = {G,P,O,G,R,O,G,K,R,G,K,T,G,L,K,G,D,I,G,P,O,G,P,O,G,P,O,G,P,O,G};
    v.insert (v.begin(), resID28, resID28 + sizeof(resID28)/sizeof(resID28[0]));
    resIDs.push_back(v);
    v.clear();
    UInt resID29[] = {G,P,O,G,L,O,G,M,L,G,Q,K,G,E,M,G,P,K,G,P,O,G,P,O,G,P,O,G,P,O,G};
    v.insert (v.begin(), resID29, resID29 + sizeof(resID29)/sizeof(resID29[0]));
    resIDs.push_back(v);
    v.clear();
    UInt resID30[] = {G,P,O,G,R,O,G,K,R,G,K,O,G,V,R,G,P,R,G,P,O,G,P,O,G,P,O,G,P,O,G};
    v.insert (v.begin(), resID30, resID30 + sizeof(resID30)/sizeof(resID30[0]));
    resIDs.push_back(v);
    v.clear();
    UInt resID31[] = {G,P,O,G,P,O,G,R,O,G,R,K,G,R,O,G,P,O,G,P,O,G,P,O,G,P,O,G,P,O,G};
    v.insert (v.begin(), resID31, resID31 + sizeof(resID31)/sizeof(resID31[0]));
    resIDs.push_back(v);
    v.clear();
    UInt resID32[] = {G,P,O,G,P,O,G,R,O,G,K,R,G,K,Q,G,Q,K,G,P,O,G,P,O,G,P,O,G,P,O,G};
    v.insert (v.begin(), resID32, resID32 + sizeof(resID32)/sizeof(resID32[0]));
    resIDs.push_back(v);
    v.clear();
    UInt resID33[] = {G,I,N,G,F,O,G,K,D,G,R,D,G,T,K,G,E,K,G,E,O,G,P,O,G,P,O,G,P,O,G,P,O,G,G};
    v.insert (v.begin(), resID33, resID33 + sizeof(resID33)/sizeof(resID33[0]));
    resIDs.push_back(v);
    v.clear();
    UInt resID34[] = {G,P,O,G,R,O,G,R,R,G,R,O,G,L,K,G,E,Q,G,P,O,G,P,O,G,P,O,G,P,O,G,Y};
    v.insert (v.begin(), resID34, resID34 + sizeof(resID34)/sizeof(resID34[0]));
    resIDs.push_back(v);
    v.clear();
    UInt resID35[] = {G,P,O,G,I,R,G,P,K,G,Q,K,G,E,O,G,L,O,G,P,O,G,P,O,G,P,O,G,P,O,G,Y};
    v.insert (v.begin(), resID35, resID35 + sizeof(resID35)/sizeof(resID35[0]));
    resIDs.push_back(v);
    v.clear();*/

	//--Mutate chains
    UInt chainNum;
    delete thePDB;
    for (UInt h = 0; h < resIDs.size(); h++)
    {
        PDBInterface* thePDB = new PDBInterface(infile);
        ensemble* theEnsemble = thePDB->getEnsemblePointer();
        molecule* pMol = theEnsemble->getMoleculePointer(0);
        protein* bundle = static_cast<protein*>(pMol);
        chainNum = bundle->getNumChains();
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
                    bundle->mutateWBC(i, j, resIDs[i][j]);
                    randomizeSideChain(bundle, i, j);
                }
            }
        }
        stringstream convert;
        string countStr, outFile;
        convert << h+1, countStr = convert.str();
        outFile = countStr;
        pdbWriter(bundle, outFile);
        delete thePDB;
        PDBInterface* thePDB2 = new PDBInterface(outFile);
        ensemble* theEnsemble2 = thePDB2->getEnsemblePointer();
        molecule* pMol2 = theEnsemble2->getMoleculePointer(0);
        protein* bundle2 = static_cast<protein*>(pMol2);
        bundle2->protOptSolvent(500);
        double Energy = bundle2->intraSoluteEnergy(false);
        cout << Energy << endl;
        pdbWriter(bundle2, outFile);
        delete thePDB2;
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
