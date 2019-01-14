//*******************************************************************************************************
//*******************************************************************************************************
//*************************************                      ********************************************
//*************************************       protEnsemble   ********************************************
//*************************************                      ********************************************
//*******************************************************************************************************
//*********************** -generate multiple conformations of input fold- *******************************
//*******************************************************************************************************

/////// Just specify a infile and outfile, it will optimize to a generally effective minimum.

#include "ensemble.h"
#include "PDBInterface.h"
#include <sstream>
#include <time.h>
#include <unistd.h>

int main (int argc, char* argv[])
{
	if (argc !=2)
	{	cout << "protEnsemble <inFile.pdb>" << endl;
		exit(1); }
	srand (getpid());
	
	PDBInterface* thePDB = new PDBInterface(infile);
	ensemble* theEnsemble = thePDB->getEnsemblePointer();
	molecule* pMol = theEnsemble->getMoleculePointer(0);
	protein* _prot = static_cast<protein*>(pMol);
	
	
	return 0;
}

vector <double> getRandPhiPsifromRPTBin(UInt _RPT)
{
	vector <double> angles(2);
	int b, psi, phi;

	if (_RPT == 1){
		do{
			b = 153 + (rand() % 27);
			phi = (rand() % 360)-180;
			psi = -1 * phi + b;
		}while(psi > 180 || psi < -180);
		angles[0] = phi;
		angles[1] = psi;
		return angles;
	}
	if (_RPT == 2){
		do{
			b = 124 + (rand() % 29);
			phi = (rand() % 360)-180;
			psi = -1 * phi + b;
		}while(psi > 180 || psi < -180);
		angles[0] = phi;
		angles[1] = psi;
		return angles;
	}
	if (_RPT == 3){
		do{
			b = 95 + (rand() % 29);
			phi = (rand() % 360)-180;
			psi = -1 * phi + b;
		}while(psi > 180 || psi < -180);
		angles[0] = phi;
		angles[1] = psi;
		return angles;
	}
	if (_RPT == 4){
		do{
			b = 58 + (rand() % 37);
			phi = (rand() % 360)-180;
			psi = -1 * phi + b;
		}while(psi > 180 || psi < -180);
		angles[0] = phi;
		angles[1] = psi;
		return angles;
	}
	if (_RPT == 5){
		do{
			b = (rand() % 58);
			phi = (rand() % 360)-180;
			psi = -1 * phi + b;
		}while(psi > 180 || psi < -180);
		angles[0] = phi;
		angles[1] = psi;
		return angles;
	}
	if (_RPT == 6){
		do{
			b = -58 + (rand() % 58);
			phi = (rand() % 360)-180;
			psi = -1 * phi + b;
		}while(psi > 180 || psi < -180);
		angles[0] = phi;
		angles[1] = psi;
		return angles;
	}
	if (_RPT == 7){
		do{
			b = -95 + (rand() % 37);
			phi = (rand() % 360)-180;
			psi = -1 * phi + b;
		}while(psi > 180 || psi < -180);
		angles[0] = phi;
		angles[1] = psi;
		return angles;
	}
	if (_RPT == 8){
		do{
			b = -124 + (rand() % 29);
			phi = (rand() % 360)-180;
			psi = -1 * phi + b;
		}while(psi > 180 || psi < -180);
		angles[0] = phi;
		angles[1] = psi;
		return angles;
	}
	if (_RPT == 9){
		do{
			b = -153 + (rand() % 29);
			phi = (rand() % 360)-180;
			psi = -1 * phi + b;
		}while(psi > 180 || psi < -180);
		angles[0] = phi;
		angles[1] = psi;
		return angles;
	}
	if (_RPT == 10){
		do{
			b = -180 + (rand() % 27);
			phi = (rand() % 360)-180;
			psi = -1 * phi + b;
		}while(psi > 180 || psi < -180);
		angles[0] = phi;
		angles[1] = psi;
		return angles;
	}
}

