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
	if (argc !=6)
	{   cout << "protMin <clearclashes(t/f)> <minimize(t/f)> <minbackboneRelax(t/f)> <inFile.pdb> <outFile.pdb>" << endl;
		exit(1); }

	string clearclashes = argv[1];
	string minimize = argv[2];
	string bbrelax = argv[3];
	string infile = argv[4];
	string outFile = argv[5];

	PDBInterface* thePDB = new PDBInterface(infile);
	ensemble* theEnsemble = thePDB->getEnsemblePointer();
	molecule* pMol = theEnsemble->getMoleculePointer(0);
	protein* _prot = static_cast<protein*>(pMol);
	bool min = true;
	if (minimize == "f"){min = false;}
	bool backbone = false;
	if (bbrelax == "t"){backbone = true;}
	bool clash = false;
	if (clearclashes == "t"){clash = true;}
	clock_t start, end;
	double cpu_time_used;

	residue::setElectroSolvationScaleFactor(1.0);
	residue::setHydroSolvationScaleFactor(1.0);
	residue::setPolarizableElec(true);
	amberElec::setScaleFactor(1.0);
	amberVDW::setScaleFactor(1.0);
	residue::setTemperature(300);
  
	if (clash){
		UInt startclashes = _prot->getNumHardClashes();
		cout << "Starting Clashes: " << startclashes << endl;
		_prot->protRelax(1000, backbone);
		UInt endclashes = _prot->getNumHardClashes();
		cout << "Clashes cleared: " << (int)startclashes-(int)endclashes << endl;
	}
	if (min){
		double startE = _prot->protEnergy();
		cout << "Starting Energy: " << startE << " kcal/mol" << endl;
		start = clock();
		_prot->protMin(backbone);
		end = clock();
		double endE = _prot->protEnergy();
		cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
		cout << "Ending Energy: "  << endE << " kcal/mol" << endl;
		cout << "delta Energy: " << endE-startE << " kcal/mol" << endl;
		cout << "time: " << cpu_time_used << endl;
	}
	pdbWriter(_prot, outFile);

	return 0;
}

