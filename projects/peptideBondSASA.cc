//*******************************************************************************************************
//*******************************************************************************************************
//********************************                       ************************************************
//********************************   amberAnalyzer 1.0   ************************************************
//********************************                       ************************************************
//*******************************************************************************************************
//***************   -Stability Selective Protein Evolution in Implicit Solvent-   ***********************
//*******************************************************************************************************

/////// Just specify infile structure and it will evolve in hetero-oligameric stability

//--Included files and functions-------------------------------------------------------------------------
#include <iostream>
#include <fstream>
#include <string>
#include <time.h>
#include <sstream>
#include "ensemble.h"
#include "PDBInterface.h"


//--Program setup----------------------------------------------------------------------------------------
int main (int argc, char* argv[])
{
	//--Running parameters
	if (argc !=1)
	{
		cout << "amberAnalyzer" << endl;
		exit(1);
	}
	enum aminoAcid {A,R,N,D,Dh,C,Cx,Q,E,Eh,Hd,He,Hn,Hp,I,L,K,M,F,P,O,S,T,W,Y,V,G,dA,dR,dN,dD,dDh,dC,dCx,dQ,dE,dEh,dHd,dHe,dHn,dHp,dI,dL,dK,dM,dF,dP,dO,dS,dT,dW,dY,dV};
	//string aminoAcidString[] = {"A","R","N","D","Dh","C","Cx","Q","E","Eh","Hd", "He","Hn","Hp","I","L","K","M","F","P","O","S","T","W","Y", "V","G","dA","dR","dN","dD","dDh","dC","dCx","dQ","dE","dEh","dHd","dHe","dHn","dHp","dI","dL","dK","dM","dF","dP","dO","dS","dT","dW","dY","dV"};
	residue::setCutoffDistance(10.0);
	rotamer::setScaleFactor(0.0);
	amberVDW::setScaleFactor(1.0);
	amberVDW::setRadiusScaleFactor(1.0);
	amberVDW::setLinearRepulsionDampeningOff();
	amberElec::setScaleFactor(1.0);
	srand (time(NULL));

	//--initialize variables for loop
	UInt chainNum, resNum, name;
	double SA;
	stringstream convert;
	string inFrame, outFile;
	//ofstream proteinFile("alaSASA");
	//ofstream residueFile("residue");
	//delete bundle;

	//--Run multiple independent evolutions
	for (UInt a = 4000; a < 4501; a++)
	{ 
		name = a;
		stringstream convert; 
		string countstr;
		convert << name, countstr = convert.str();
		inFrame = "pC_Frame." + countstr + ".pdb";
		PDBInterface* theFramePDB = new PDBInterface(inFrame);
		ensemble* theFrameEnsemble = theFramePDB->getEnsemblePointer();
		molecule* frameMol = theFrameEnsemble->getMoleculePointer(0);
		protein* frame = static_cast<protein*>(frameMol);
		frame->silenceMessages();
		chainNum = frame->getNumChains();
		for (UInt i = 0; i < chainNum; i++)
		{
			resNum = frame->getNumResidues(i);
			for (UInt j = 0; j < resNum; j++)
			{
				frame->activateForRepacking(i, j);
				frame->mutateWBC(i, j, A);
			}
		}
		
		cout << a << " " << (a*.02) << " ";
		frame->initializeSpherePoints();
		frame->removeSpherePoints();
		for (UInt i = 0; i < chainNum; i++)
		{
			resNum = frame->getNumResidues(i);
			for (UInt j = 0; j < resNum-1; j++)
			{
				SA = ((frame->tabulateSurfaceArea(i, j, 2))+(frame->tabulateSurfaceArea(i, j+1, 0)));
				cout << SA << " " ;
			}
		}	
		cout << endl;
		delete theFramePDB;
	}
	return 0;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////
//////// functions //////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////

	
