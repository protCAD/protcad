//*******************************************************************************************************
//*******************************************************************************************************
//*****************************                         *************************************************
//*****************************        protEnergy       *************************************************
//*****************************                         *************************************************
//*******************************************************************************************************
//*******************************************************************************************************

#include <iostream>
#include <string>
#include "ensemble.h"
#include "PDBInterface.h"
#include <time.h>

//--Program setup-------------------------------------------------------------
int main (int argc, char* argv[])
{	

	if (argc !=3)
	{
		cout << "protEnergy <inFile.pdb> <predictTM(t/f)>" << endl;
		exit(1);
	}
	string infile = argv[1];
	string predicttm = argv[2];

	bool predictTM = false;
	if (predicttm == "t"){predictTM = true;}

	PDBInterface* thePDB = new PDBInterface(infile);
	ensemble* theEnsemble = thePDB->getEnsemblePointer();
	molecule* pMol = theEnsemble->getMoleculePointer(0);
	protein* prot = static_cast<protein*>(pMol);
	
	amberElec::setScaleFactor(1.0);
	amberVDW::setScaleFactor(1.0);
	residue::setElectroSolvationScaleFactor(1.0);
	residue::setHydroSolvationScaleFactor(1.0);
	residue::setPolarizableElec(true);
	residue::setEntropyFactor(1.0);
	
	double Energy = 0.0;
	if (predictTM){
		double standardE;
		double tm;
		double temp;
		bool notm = true;
		residue::setTemperature(273);
		for (UInt i=0; i < 100; i++)
		{
			Energy = prot->protEnergy();
			temp = residue::getTemperature();
			if (Energy > 0 && notm) {tm = temp-273; notm = false;}
			if (temp == 300){standardE = Energy;}
			prot->setMoved(true,0);
			residue::setTemperature(274+i);
		}
		cout << infile << " " << standardE << " kcal/mol at 300K Predicted TM: " << tm << "C" << endl;
	}
	else{
		residue::setTemperature(300);
		cout << infile << " " << prot->protEnergy() << " kcal/mol at 300K" << endl;
	}
	return 0;
}
