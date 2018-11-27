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
	string infile = argv[1];
	PDBInterface* thePDB = new PDBInterface(infile);
	ensemble* theEnsemble = thePDB->getEnsemblePointer();
	molecule* pMol = theEnsemble->getMoleculePointer(0);
	protein* bundle = static_cast<protein*>(pMol);
	clock_t start, end;
	double cpu_time_used;
	
	residue::setElectroSolvationScaleFactor(1.0);
	residue::setHydroSolvationScaleFactor(1.0);
	amberElec::setScaleFactor(1.0);
	amberVDW::setScaleFactor(1.0);
	
	start = clock();
	double Energy = bundle->protEnergy();
	end = clock();
	cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
	cout << "protEnergy: " << Energy << " kcal/mol time: " << cpu_time_used << endl;
	string outFile = "protEnergy_out.pdb";
	pdbWriter(bundle, outFile);
	
	start = clock();
	Energy = bundle->protEnergy();
	end = clock();
	cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
	cout << "protEnergy: " << Energy << " kcal/mol time: " << cpu_time_used << endl;
	outFile = "protEnergy_out2.pdb";
	pdbWriter(bundle, outFile);
	
	bundle->setMoved(true);
	start = clock();
	Energy = bundle->protEnergy();
	end = clock();
	cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
	cout << "protEnergy: " << Energy << " kcal/mol time: " << cpu_time_used << endl;
	outFile = "protEnergy_out3.pdb";
	pdbWriter(bundle, outFile);
	
	
	
	/*start = clock();
	UInt clashes = bundle->getNumHardClashes();
	end = clock();
	cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
	cout << clashes << " clashes time: " << cpu_time_used << endl;
	string outFile = "protRelax_out.pdb";
	pdbWriter(bundle, outFile);
	
	start = clock();
	clashes = bundle->getNumHardClashes();
	end = clock();
	cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
	cout << clashes << " clashes time: " << cpu_time_used << endl;
	outFile = "protRelax_out.pdb";
	pdbWriter(bundle, outFile);*/
	return 0;
}
