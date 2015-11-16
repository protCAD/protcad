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
        if (inFrame.find(".evo.pdb") != std::string::npos)
        {
            count++;
            PDBInterface* theModelPDB = new PDBInterface(inFrame);
            ensemble* theModelEnsemble = theModelPDB->getEnsemblePointer();
            molecule* modelMol = theModelEnsemble->getMoleculePointer(0);
            protein* model = static_cast<protein*>(modelMol);
            model->silenceMessages();
            residue::setCutoffDistance(9.0);
            rotamer::setScaleFactor(0.0);
            amberVDW::setScaleFactor(1.0);
            amberVDW::setRadiusScaleFactor(1.0);
            amberVDW::setLinearRepulsionDampeningOff();
            amberElec::setScaleFactor(1.0);

            Energy = model->protEnergy();
            if (Energy < -322)
            {
                stringstream convert;
                string outFile;
                UInt name = count;
                convert << name, outFile = convert.str();
                string tempModel = outFile + "_best.evo.pdb";
                pdbWriter(model, tempModel);
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
