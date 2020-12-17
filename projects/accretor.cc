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

	if (argc !=2)
	{
		cout << "protEnergy <inFile.pdb>" << endl;
		exit(1);
	}
	clock_t start, end;
	double cpu_time_used;
	string infile = argv[1];
	start = clock();
	PDBInterface* thePDB = new PDBInterface(infile);
	ensemble* theEnsemble = thePDB->getEnsemblePointer();
	molecule* pMol = theEnsemble->getMoleculePointer(0);
	protein* bundle = static_cast<protein*>(pMol);
	end = clock();
	
	
	residue::setElectroSolvationScaleFactor(1.0);
	residue::setHydroSolvationScaleFactor(1.0);
	residue::setPolarizableElec(true);
	amberElec::setScaleFactor(1.0);
	amberVDW::setScaleFactor(1.0);
	residue::setTemperature(300);
	
	
	double Energy = bundle->protEnergy();
	
	cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
	cout << infile << " " << Energy << " kcal/mol time: " << cpu_time_used << endl;
	string outFile = infile;
	pdbWriter(bundle, outFile);
	
	start = clock();
	UInt clashes = bundle->getNumHardClashes();
	end = clock();
	cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
	cout << "Clashes: " << clashes << " clashes time: " << cpu_time_used << endl;
	
	start = clock();
	clashes = bundle->getNumHardBackboneClashes();
	end = clock();
	cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
	cout << "Backbone Clashes: " << clashes << " clashes time: " << cpu_time_used << endl;
	return 0;
}
