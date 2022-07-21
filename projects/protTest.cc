//*******************************************************************************************************
//*******************************************************************************************************
//*****************************                         *************************************************
//*****************************        protTest         *************************************************
//*****************************                         *************************************************
//*******************************************************************************************************
//*******************************************************************************************************

#include <iostream>
#include <string>
#include <dirent.h>
#include "ensemble.h"
#include "PDBInterface.h"
#include <time.h>

//--Program setup-------------------------------------------------------------
int main (int argc, char* argv[])
{	
	if (argc !=1)
	{
		cout << "protTest" << endl;
		exit(1);
	}
	string aminoAcidString[] = {"A","R","N","D","D","C","C","C","Q","E","E","H","H","H","I","L","K","M","F","P","O","S","T","W","Y","V","G","A","R","N","D","D","C","C","C","Q","E","E","H","H","H","I","L","K","M","F","P","O","S","T","W","Y","V","A","R","N","D","D","C","C","C","Q","E","E","H","H","H","I","L","K","M","F","P","O","S","T","W","Y","V","G","A","R","N","D","D","C","C","C","Q","E","E","H","H","H","I","L","K","M","F","P","O","S","T","W","Y","V","A","R","N","D","D","C","C","C","Q","E","E","H","H","H","I","L","K","M","F","P","O","S","T","W","Y","V","G","A","R","N","D","D","C","C","C","Q","E","E","H","H","H","I","L","K","M","F","P","O","S","T","W","Y","V"};
	fstream aa;
	aa.open ("bests.faa", fstream::in | fstream::out | fstream::app);
	string inFrame;
	DIR *pdir;
	struct dirent *pent;
	pdir=opendir(".");
	while ((pent=readdir(pdir)))
	{
		inFrame = pent->d_name;
        if (inFrame.find(".evo.pdb") != std::string::npos)
		{
			PDBInterface* theFramePDB = new PDBInterface(inFrame);
			ensemble* theFrameEnsemble = theFramePDB->getEnsemblePointer();
			molecule* frameMol = theFrameEnsemble->getMoleculePointer(0);
			protein* frame = static_cast<protein*>(frameMol);
			double complexE = frame->protEnergy();
			double receptorE = frame->protEnergy(0);
			double ligandE = frame->protEnergy(1);
			double bindingEnergy = complexE-(ligandE+receptorE);
			if (bindingEnergy < -400 && complexE < -1500)
			{
				bool write = true;
				UInt numRes = frame->getNumResidues(0);
				for (UInt j = 0; j < numRes; j++)
				{	
					if (write){
						aa << ">" << inFrame << endl;
						write = false;
					}
					UInt restype = frame->getTypeFromResNum(0,j);
					aa << aminoAcidString[restype];
				}
				aa << endl;
			}
			delete theFramePDB;
		}
	}
	aa.close();
	closedir(pdir);
	return 0;
}
