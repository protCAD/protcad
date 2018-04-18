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
    enum aminoAcid {A,R,N,D,Dh,C,Cx,Cf,Q,E,Eh,Hd,He,Hp,I,L,K,M,F,P,O,S,T,W,Y,V,G,dA,dR,dN,dD,dDh,dC,dCx,dCf,dQ,dE,dEh,dHd,dHe,dHp,dI,dL,dK,dM,dF,dP,dO,dS,dT,dW,dY,dV,Csf,Sf4,Hca,Eoc,Oec};
    string aminoAcidString[] = {"A","R","N","D","Dh","C","Cx","Cf","Q","E","Eh","Hd","He","Hp","I","L","K","M","F","P","O","S","T","W","Y","V","G","dA","dR","dN","dD","dDh","dC","dCx","dCf","dQ","dE","dEh","dHd","dHe","dHp","dI","dL","dK","dM","dF","dP","dO","dS","dT","dW","dY","dV","Csf","Sf4","Hca","Eoc","Oec"};
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
        if (inFrame.find(".trim") != std::string::npos)
		{
			PDBInterface* theFramePDB = new PDBInterface(inFrame);
			ensemble* theFrameEnsemble = theFramePDB->getEnsemblePointer();
			molecule* frameMol = theFrameEnsemble->getMoleculePointer(0);
			protein* frame = static_cast<protein*>(frameMol);
			for (UInt i = 0; i < frame->getNumChains(); i++)
			{
                for (UInt j = 1; j < frame->getNumResidues(i)-4; j++)
				{
                    if (frame->getNumResidues(i) > 4)
                    {
                        if (frame->getTypeFromResNum(i,j) == C && frame->getTypeFromResNum(i,j+3)== C)
                        {
                            cout << inFrame << " " << i << " " << aminoAcidString[frame->getTypeFromResNum(i,j+1)] << " " << aminoAcidString[frame->getTypeFromResNum(i,j+2)] << " " << frame->getPhi(i,j) << " " << frame->getPsi(i,j) << " " << frame->getPhi(i,j+1) << " " << frame->getPsi(i,j+1) << " " << frame->getPhi(i,j+2) << " " << frame->getPsi(i,j+2) << " " << frame->getPhi(i,j+3) << " " << frame->getPsi(i,j+3) << endl;
                        }
                    }
                }
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

	
