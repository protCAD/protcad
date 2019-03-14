//*******************************************************************************************************
//*******************************************************************************************************
//********************************                       ************************************************
//********************************        protDB         ************************************************
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
	if (argc !=2)
	{
		cout << "protDB <infile.pdb>" << endl;
		exit(1);
	}
	
	residue::setElectroSolvationScaleFactor(1.0);
	residue::setHydroSolvationScaleFactor(0.0);
	amberElec::setScaleFactor(1.0);
	amberVDW::setScaleFactor(1.0);
	residue::setPolarizableElec(true);
	string inFrame, infile = argv[1];
	
	PDBInterface* thePDB = new PDBInterface(infile);
	ensemble* theEnsemble = thePDB->getEnsemblePointer();
	molecule* Mol = theEnsemble->getMoleculePointer(0);
	protein* _prot = static_cast<protein*>(Mol);
	
	DIR *pdir;
	struct dirent *pent;
	pdir=opendir(".");
	while ((pent=readdir(pdir)))
	{
		inFrame = pent->d_name;
        if (inFrame.find(".fold.pdb") != std::string::npos)
		{
			PDBInterface* theFramePDB = new PDBInterface(inFrame);
			ensemble* theFrameEnsemble = theFramePDB->getEnsemblePointer();
			molecule* frameMol = theFrameEnsemble->getMoleculePointer(0);
			protein* frame = static_cast<protein*>(frameMol);
			cout << inFrame << " " << frame->protEnergy() << " " << _prot->getRMSD(frame) << endl;
			delete theFramePDB;
		}
	}
	closedir(pdir);
	return 0;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////
//////// functions //////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////

	
