//*******************************************************************************************************
//*******************************************************************************************************
//*************************************                      ********************************************
//*************************************       protFolder     ********************************************
//*************************************                      ********************************************
//*******************************************************************************************************
//******** -sidechain and backbone optimization with a burial-based scaling of electrostatics- **********
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
	{	cout << "protFolder <inFile.pdb>" << endl;
		exit(1); }

	//--Assign model nomenclature
	string infile = argv[1];
	stringstream convert;
	string startstr;
	srand (getpid());
	UInt name = rand() % 100000000;
	convert << name, startstr = convert.str();
	string foldModel = startstr + "_fold.pdb";
	
	//--Build protein object and calculate starting energy
	PDBInterface* thePDB = new PDBInterface(infile);
	ensemble* theEnsemble = thePDB->getEnsemblePointer();
	molecule* pMol = theEnsemble->getMoleculePointer(0);
	protein* _prot = static_cast<protein*>(pMol);
	double Energy = 1E10;
	double newEnergy = Energy;
	
	//--sample folds and save energetic minima
	while (true)
	{
		_prot->protSampling();
		newEnergy = _prot->protEnergy();
		cout << startstr << " " << newEnergy;
		if (newEnergy < Energy){
			pdbWriter(_prot, foldModel);
			cout << " minima!" << endl;
			Energy = newEnergy;
		}
		else{cout << endl;}
	}
	return 0;
}
