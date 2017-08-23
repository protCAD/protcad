//*******************************************************************************************************
//*******************************************************************************************************
//********************************                       ************************************************
//********************************   database_phipsi 1.0   ************************************************
//********************************                       ************************************************
//*******************************************************************************************************
//***************   -Stability Selective Protein Evolution in Implicit Solvent-   ***********************
//*******************************************************************************************************

/////// Just specify infile structure and it will evolve in hetero-oligameric stability

//--Included files and functions-------------------------------------------------------------------------
#include <iostream>
#include <dirent.h>
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
	enum aminoAcid {A,R,N,D,Dh,C,Cx,Cf,Q,E,Eh,Hd,He,Hp,I,L,K,M,F,P,O,S,T,W,Y,V,G,dA,dR,dN,dD,dDh,dC,dCx,dCf,dQ,dE,dEh,dHd,dHe,dHp,dI,dL,dK,dM,dF,dP,dO,dS,dT,dW,dY,dV,Csf,Hca,Oec};
	string aminoAcidString[] = {"A","R","N","D","Dh","C","Cx","Cf","Q","E","Eh","Hd","He","Hn","Hp","I","L","K","M","F","P","O","S","T","W","Y","V","G","dA","dR","dN","dD","dDh","dC","dCx","dQ","dE","dEh","dHd","dHe","dHn","dHp","dI","dL","dK","dM","dF","dP","dO","dS","dT","dAT","dW","dY","dV","Hce","Pch","Csf","dCf"};
	residue::setCutoffDistance(8.0);
	rotamer::setScaleFactor(0.0);
	amberVDW::setScaleFactor(0.0);
	amberVDW::setRadiusScaleFactor(0.0);
	amberVDW::setLinearRepulsionDampeningOff();
	amberElec::setScaleFactor(0.0);
	srand (time(NULL));

	double phi, psi;
	string inFrame;
	DIR *pdir;
	struct dirent *pent;
	pdir=opendir(".");
	while ((pent=readdir(pdir)))
	{
		inFrame = pent->d_name;
		if (inFrame.find(".ent") != std::string::npos)
		{
			PDBInterface* theFramePDB = new PDBInterface(inFrame);
			ensemble* theFrameEnsemble = theFramePDB->getEnsemblePointer();
			molecule* frameMol = theFrameEnsemble->getMoleculePointer(0);
			protein* frame = static_cast<protein*>(frameMol);
			for (UInt i = 0; i < frame->getNumChains(); i++)
			{
				for (UInt j = 0; j < frame->getNumResidues(i); j++)
				{
					if (j != 0 && j != frame->getNumResidues(i)-1)
					{	phi = frame->getPhi(i,j), psi = frame->getPsi(i,j);
						cout << inFrame << " " << j << " " << phi << " " << psi << " " << phi-psi << endl; }
				}
				break;
			}
			delete theFramePDB;
		}
	}
	closedir(pdir);
	return 0;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////
//////// functions //////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////

	
