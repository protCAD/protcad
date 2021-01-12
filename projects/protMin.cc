//*******************************************************************************************************
//*******************************************************************************************************
//*************************************                      ********************************************
//*************************************       protMin        ********************************************
//*************************************                      ********************************************
//*******************************************************************************************************
//****************************** sidechain and backbone optimization ************************************
//*******************************************************************************************************

/////// Just specify a infile and outfile, it will optimize to a generally effective minimum.

#include "ensemble.h"
#include "PDBInterface.h"
#include <sstream>
#include <time.h>
#include <unistd.h>
int main (int argc, char* argv[])
{
	if (argc !=4)
	{   cout << "protMin <backboneRelax(t/f)> <inFile.pdb> <outFile.pdb>" << endl;
		exit(1); }

	string relax = argv[1];
	string infile = argv[2];
	string outFile = argv[3];
	PDBInterface* thePDB = new PDBInterface(infile);
	ensemble* theEnsemble = thePDB->getEnsemblePointer();
	molecule* pMol = theEnsemble->getMoleculePointer(0);
	protein* _prot = static_cast<protein*>(pMol);
	bool backbone = false;
	if (relax == "t"){backbone = true;}
	clock_t start, end;
	double cpu_time_used;
	int seed = (int)getpid()*(int)gethostid(); srand (seed);

	residue::setElectroSolvationScaleFactor(1.0);
	residue::setHydroSolvationScaleFactor(1.0);
	residue::setPolarizableElec(true);
	amberElec::setScaleFactor(1.0);
	amberVDW::setScaleFactor(1.0);
	residue::setTemperature(300);

	double startE = _prot->protEnergy();
	cout << "Starting Energy: " << startE << " kcal/mol" << endl;
	UInt startclashes = _prot->getNumHardClashes();
	start = clock();
	_prot->cofactorRelax(1000);
	end = clock();
	UInt endclashes = _prot->getNumHardClashes();
	double endE = _prot->protEnergy();
	cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
	cout << "Ending Energy: "  << endE << " kcal/mol" << endl;
	cout << "delta Energy: " << endE-startE << " kcal/mol" << endl;
	cout << "Clashes cleared: " << (int)startclashes-(int)endclashes << endl;
	cout << "time: " << cpu_time_used << endl;
	pdbWriter(_prot, outFile);

	return 0;
}

