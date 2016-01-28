//*******************************************************************************************************
//*******************************************************************************************************
//***********************************                          ******************************************
//***********************************      protSorter 1.0      ******************************************
//***********************************                          ******************************************
//*******************************************************************************************************
//***************   -                                                         -   ***********************
//*******************************************************************************************************


//--Included files and functions-------------------------------------------------------------------------
#include <iostream>
#include <string>
#include <time.h>
#include <dirent.h>
#include <sstream>
#include <unistd.h>
#include "ensemble.h"
#include "PDBInterface.h"


//--Program setup----------------------------------------------------------------------------------------
int main (int argc, char* argv[])
{
    enum aminoAcid {A,R,N,D,Dh,C,Cx,Cf,Q,E,Eh,Hd,He,Hn,Hp,I,L,K,M,F,P,O,S,T,W,Y,V,G,dA,dR,dN,dD,dDh,dC,dCx,dQ,dE,dEh,dHd,dHe,dHn,dHp,dI,dL,dK,dM,dF,dP,dO,dS,dT,dAT,dW,dY,dV,Hce,Pch};
	//--Running parameters
        if (argc !=1)
	{
                cout << "protSorter" << endl;
		exit(1);
	}
    double Energy;
    string inFrame;
    DIR *pdir;
    struct dirent *pent;
    pdir=opendir(".");
    UInt count = 0;

    //--get sequence evolution results for position
    while ((pent=readdir(pdir)))
    {
        inFrame = pent->d_name;
        if (inFrame.find(".pdb") != std::string::npos)
        {
            count++;
            PDBInterface* theModelPDB = new PDBInterface(inFrame);
            ensemble* theModelEnsemble = theModelPDB->getEnsemblePointer();
            molecule* modelMol = theModelEnsemble->getMoleculePointer(0);
            protein* model = static_cast<protein*>(modelMol);
            model->silenceMessages();
            residue::setCutoffDistance(9.0);
            rotamer::setScaleFactor(0.0);
            amberVDW::setScaleFactor(0.95);
            amberVDW::setRadiusScaleFactor(1.0);
            amberVDW::setLinearRepulsionDampeningOff();
            amberElec::setScaleFactor(1.0);

            Energy = model->protEnergy();
            UInt hces = 0;
            UInt chainNum = model->getNumChains();
            for (UInt i = 0; i < chainNum; i ++)
            {
                UInt resNum = model->getNumResidues(i);
                for (UInt j = 0; j < resNum; j ++)
                {
                    UInt type = model->getTypeFromResNum(i,j);
                    if (type == Hce)
                    {
                        hces++;
                    }
                }
            }
            if (hces == 7)
            {
                stringstream convert;
                string outFile;
                UInt name = count;
                convert << name, outFile = convert.str();
                string tempModel = outFile + "_7.pdb";
                pdbWriter(model, tempModel);
                cout << inFrame << " " << Energy << endl;
            }
            delete theModelPDB;
        }
    }
    closedir(pdir);
    return 0;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////
//////// functions //////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////
