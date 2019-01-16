//*******************************************************************************************************
//*******************************************************************************************************
//*************************************                      ********************************************
//*************************************       protMin        ********************************************
//*************************************                      ********************************************
//*******************************************************************************************************
//******** -sidechain and backbone minimization with a burial-based scaling of electrostatics- **********
//*******************************************************************************************************

/////// Just specify a infile and outfile, it will optimize to a generally effective minimum.

#include "ensemble.h"
#include "PDBInterface.h"
#include <sstream>
#include <time.h>
int main (int argc, char* argv[])
{
	if (argc !=3)
	{   cout << "protMin <infile.pdb> <outfile.pdb>" << endl;
		exit(1); }

	string infile = argv[1];
	string outfile = argv[2];
	double meanEnergy, sumEnergy = 0.0, bestEnergy = 1E10;
	UInt size = 10;
	vector <double> Energies(size);
	
	//#pragma omp parallel for
	for (UInt i = 0; i < size; i++)
	{
		PDBInterface* thePDB = new PDBInterface(infile);
		ensemble* theEnsemble = thePDB->getEnsemblePointer();
		molecule* pMol = theEnsemble->getMoleculePointer(0);
		protein* _prot = static_cast<protein*>(pMol);
		_prot->protMin();
		double Energy = _prot->protEnergy();
		sumEnergy += Energy;
		Energies[i] = Energy;
		if (Energy < bestEnergy){
			bestEnergy = Energy;
			pdbWriter(_prot, outfile);
		}
		delete thePDB;
	}
	meanEnergy = sumEnergy/size;
	sumEnergy = 0.0;
	for (UInt i = 0; i < Energies.size(); i++)
	{
		sumEnergy += pow(Energies[i]-meanEnergy,2);
	}
	sumEnergy /= size;
	double stdev = sqrt(sumEnergy);
	cout << "Energy: " << bestEnergy << "Mean: " << meanEnergy << " +/- " << stdev << endl;

	return 0;
}

