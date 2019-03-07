//*******************************************************************************************************
//*******************************************************************************************************
//*************************************                      ********************************************
//*************************************       protOpt        ********************************************
//*************************************                      ********************************************
//*******************************************************************************************************
//******** -sidechain and backbone optimization with a burial-based scaling of electrostatics- **********
//*******************************************************************************************************

/////// Just specify a infile and outfile, it will optimize to a generally effective minimum.

#include "ensemble.h"
#include "PDBInterface.h"
#include <sstream>
#include <time.h>
int main (int argc, char* argv[])
{
	if (argc !=3)
	{   cout << "protMin <infile.pdb> <temp(K)>" << endl;
		exit(1); }

	string infile = argv[1];
	string temp = argv[2];
	double temperature = atof(temp.c_str());
	string iterate;
	double meanEnergy, sumEnergy = 0.0;
	UInt size = 10;
	vector <double> Energies(size);
	residue::setTemperature(temperature);

	for (UInt i = 0; i < size; i++)
	{
		PDBInterface* thePDB = new PDBInterface(infile);
		ensemble* theEnsemble = thePDB->getEnsemblePointer();
		molecule* pMol = theEnsemble->getMoleculePointer(0);
		protein* _prot = static_cast<protein*>(pMol);
		_prot->protMin(true);
		double Energy = _prot->protEnergy();
		sumEnergy += Energy;
		Energies[i] = Energy;
		stringstream convert;
		convert << i+1, iterate = convert.str();
		string minModel = iterate + "_min400.pdb";
		pdbWriter(_prot, minModel);
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
	cout << "Energy: " << meanEnergy << " +/- " << stdev << endl;

	return 0;
}

