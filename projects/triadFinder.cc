//*******************************************************************************************************
//*******************************************************************************************************
//********************************                       ************************************************
//********************************     triadFiner 1.0    ************************************************
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
        cout << "triadFinder" << endl;
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
    UInt name;
    string inFrame;
    //ofstream donors("donors");
    //ofstream histidines("histidines");
    //ofstream acceptors("acceptors");
    //ofstream proteinFile("triads5all");
    double dist;
	//ofstream residueFile("residue");
	//delete bundle;

	//--Run multiple independent evolutions
    cout << "frame time distance" << endl;
	for (UInt a = 1; a < 19981; a++)
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


		dblVec Ocoords = frame->getCoords(0, 90, "HH");
        dblVec Ncoords = frame->getCoords(0, 14, "NE2");
        dist = CMath::distance(Ocoords, Ncoords);
        cout << a << " " << (a*.02) << " " << dist << endl;
		delete theFramePDB;
	}
	cout << "Complete" << endl << endl;
	return 0;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////
//////// functions //////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////

	
