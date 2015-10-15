//*******************************************************************************************************
//*******************************************************************************************************
//************************************                     **********************************************
//************************************  dielectricFit 1.0  **********************************************
//************************************                     **********************************************
//*******************************************************************************************************
//**********************  -												- *******************************
//*******************************************************************************************************

/////// Just specify a infile and outfile, it will output modified structure.

//--Program setup----------------------------------------------------------------------------------------
#include <iostream>
#include <fstream>
#include <string>
#include <time.h>
#include <sstream>
#include "ensemble.h"
#include "PDBInterface.h"
#include "residueTemplate.h"

void randomizeSideChain(protein* _prot, UInt _chainIndex, UInt _resIndex);

int main (int argc, char* argv[])
{
	//--Running parameters
	if (argc !=1)
	{
		cout << "dielectricFit" << endl;
		exit(1);
	}
	enum aminoAcid {A,R,N,D,Dh,C,Cx,Q,E,Eh,Hd,He,Hn,Hp,I,L,K,M,F,P,O,S,T,W,Y,V,G,dA,dR,dN,dD,dDh,dC,dCx,dQ,dE,dEh,dHd,dHe,dHn,dHp,dI,dL,dK,dM,dF,dP,dO,dS,dT,dW,dY,dV};
	//string aminoAcidString[] = {"A","R","N","D","Dh","C","Cx","Q","E","Eh","Hd", "He","Hn","Hp","I","L","K","M","F","P","O","S","T","W","Y", "V","G","dA","dR","dN","dD","dDh","dC","dCx","dQ","dE","dEh","dHd","dHe","dHn","dHp","dI","dL","dK","dM","dF","dP","dO","dS","dT","dW","dY","dV"};
    residue::setCutoffDistance(9.0);
	rotamer::setScaleFactor(0.0);
	amberVDW::setScaleFactor(1.0);
	amberVDW::setRadiusScaleFactor(1.0);
	amberVDW::setLinearRepulsionDampeningOff();
	amberElec::setScaleFactor(1.0);
	srand (time(NULL));

	//--Initialize variables for loop
	string inFile, countstr, outFile;

	//training
	//UInt resID[] = {K,E,K,K,E,E,E,E,K,E,E,E,K,K,K,E,E,E,K,K,E};
	UInt mutants[] = {96,45,77,27,59,79,90,112,18,29,31,87,34,49,91,16,32,53,61,86,78};

	//test
	//UInt resID[] = {E,E,E,E,K,K,K,K,E,K,E,K,E,E,E,K,K,K,E,E,K};
	//UInt mutants[] = {96,119,77,21,59,79,90,112,12,23,24,87,28,49,91,10,26,53,61,86,78};

	//UInt resID[] = {K,E};
	//UInt mutants[] = {32,30};

	UInt mutantsSize = sizeof(mutants)/sizeof(mutants[0]);
	vector <double> chargeDensity(3);
	chargeDensity[0] = 0.0;
	chargeDensity[1] = 0.0;
	chargeDensity[2] = 0.0;


	//--loop
	for (UInt i = 0; i < mutantsSize; i++)
	{
		stringstream convert;
		convert << i, countstr = convert.str();
		inFile = countstr + ".opt.pdb";;
		outFile = countstr + ".opt.pdb";
		PDBInterface* theFramePDB = new PDBInterface(inFile);
		ensemble* theFrameEnsemble = theFramePDB->getEnsemblePointer();
		molecule* frameMol = theFrameEnsemble->getMoleculePointer(0);
		protein* frame = static_cast<protein*>(frameMol);
		frame->silenceMessages();
        //chargeDensity = frame->getChargeDensity(0, mutants[i], 8);
		//pdbWriter(frame, outFile);
        //cout << chargeDensity[0] << "   " << chargeDensity[1] << "   " << chargeDensity[2] << endl;
		delete theFramePDB;
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
