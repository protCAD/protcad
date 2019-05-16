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
	if (argc !=1)
	{
		cout << "protDB" << endl;
		exit(1);
	}
	
	DIR *pdir = opendir("."); struct dirent *pent; string inFrame;
	while ((pent=readdir(pdir)))
	{
		inFrame = pent->d_name;
		if (inFrame.find(".pdb") != std::string::npos)
		{
			PDBInterface* theFramePDB = new PDBInterface(inFrame);
			ensemble* theFrameEnsemble = theFramePDB->getEnsemblePointer();
			molecule* frameMol = theFrameEnsemble->getMoleculePointer(0);
			protein* frame = static_cast<protein*>(frameMol);
			for (UInt i = 0; i < frame->getNumResidues(0); i++)
			{
				cout << inFrame << " " << i << " " << frame->getBackboneSequenceType(0,i) << endl;
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

	
