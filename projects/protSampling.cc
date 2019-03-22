//*******************************************************************************************************
//*******************************************************************************************************
//*************************************                      ********************************************
//*************************************       protMin        ********************************************
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
	{   cout << "protSampling <infile.pdb> <temp(K)>" << endl;
		exit(1); }

	string infile = argv[1], temp = argv[2], iterate, samplingModel;
	double temperature = atof(temp.c_str());
	double meanEnergy, sumEnergy = 0.0;
	UInt size = 1000;
	vector <double> Energies(size);
	residue::setTemperature(temperature);

	PDBInterface* thePDB = new PDBInterface(infile);
	ensemble* theEnsemble = thePDB->getEnsemblePointer();
	molecule* pMol = theEnsemble->getMoleculePointer(0);
	protein* _prot = static_cast<protein*>(pMol);
	for (UInt i = 0; i < size; i++)
	{
		_prot->protSampling(1);
		double Energy = _prot->protEnergy();
		cout << Energy << endl;
		sumEnergy += Energy;
		Energies[i] = Energy;
		stringstream convert;
		convert << i+1, iterate = convert.str();
		samplingModel = iterate + "_" + temp + ".pdb";
		pdbWriter(_prot, samplingModel);
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

