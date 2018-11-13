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
	
	
	start = clock();
	double Energy = bundle->intraSoluteEnergy();
	end = clock();
	cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
	cout << "intraSoluteEnergy: " << Energy << " kcal/mol time: " << cpu_time_used << endl;
	
	start = clock();
	Energy = bundle->protEnergy();
	end = clock();
	cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
	cout << "protEnergy: " << Energy << " kcal/mol time: " << cpu_time_used << endl;
	
	start = clock();
	Energy = bundle->intraSoluteEnergy();
	end = clock();
	cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
	cout << "intraSoluteEnergy: " << Energy << " kcal/mol time: " << cpu_time_used << endl;
	
	start = clock();
	UInt clashes = bundle->getNumHardClashes();
	end = clock();
	cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
	cout << clashes << " clashes time: " << cpu_time_used << endl;
	
	start = clock();
	clashes = bundle->getNumHardClashes();
	end = clock();
	cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
	cout << clashes << " clashes time: " << cpu_time_used << endl;
	
	start = clock();
	Energy = bundle->protEnergy();
	end = clock();
	cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
	cout << "protEnergy: " << Energy << " kcal/mol time: " << cpu_time_used << endl;

	return 0;
}
