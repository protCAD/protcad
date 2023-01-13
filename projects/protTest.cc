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
	if (argc !=2)
	{
		cout << "protTest <inFile.pdb>" << endl;
		exit(1);
	}
	string infile = argv[1];

	
#ifdef __CUDA__
	
	PDBInterface* thePDB = new PDBInterface(infile);
	ensemble* theEnsemble = thePDB->getEnsemblePointer();
	molecule* pMol = theEnsemble->getMoleculePointer(0);
	protein* prot = static_cast<protein*>(pMol);
	prot->loadDeviceMemClash();
	clock_t start = clock();
	double Energy = prot->getNumClashesCU();
	clock_t end = clock();
	double GPUtime = double(end - start)/CLOCKS_PER_SEC;
	cout << "GPU - Energy: " << Energy << " Time(s): " << GPUtime << endl;
#endif

	PDBInterface* thePDB2 = new PDBInterface(infile);
	ensemble* theEnsemble2 = thePDB2->getEnsemblePointer();
	molecule* pMol2 = theEnsemble2->getMoleculePointer(0);
	protein* prot2 = static_cast<protein*>(pMol2);

	clock_t start2 = clock();
	Energy = prot2->getNumHardClashes();
	clock_t end2 = clock();	

	double CPUtime = double(end2 - start2)/CLOCKS_PER_SEC;
	cout << "CPU - Energy: " << Energy << " Time(s): " << CPUtime << endl;
	cout << "GPU  Speedup: " << int(CPUtime/GPUtime) << "x" << endl;
	string outpdb = "testout.pdb";
	pdbWriter(prot2, outpdb);
	/*string aminoAcidString[] = {"A","R","N","D","D","C","C","C","Q","E","E","H","H","H","I","L","K","M","F","P","O","S","T","W","Y","V","G","A","R","N","D","D","C","C","C","Q","E","E","H","H","H","I","L","K","M","F","P","O","S","T","W","Y","V","A","R","N","D","D","C","C","C","Q","E","E","H","H","H","I","L","K","M","F","P","O","S","T","W","Y","V","G","A","R","N","D","D","C","C","C","Q","E","E","H","H","H","I","L","K","M","F","P","O","S","T","W","Y","V","A","R","N","D","D","C","C","C","Q","E","E","H","H","H","I","L","K","M","F","P","O","S","T","W","Y","V","G","A","R","N","D","D","C","C","C","Q","E","E","H","H","H","I","L","K","M","F","P","O","S","T","W","Y","V"};
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
	closedir(pdir);*/
	return 0;
}
