//*******************************************************************************************************
//*******************************************************************************************************
//*****************************                         *************************************************
//*****************************        protBindingEnergy       *************************************************
//*****************************                         *************************************************
//*******************************************************************************************************
//*******************************************************************************************************

#include <iostream>
#include <string>
#include "ensemble.h"
#include "PDBInterface.h"

//--Program setup-------------------------------------------------------------
int main (int argc, char* argv[])
{	

	if (argc !=3)
	{
		cout << "protBindingEnergy <inFile.pdb> <ligand chain ID>" << endl;
		exit(1);
	}
	string infile = argv[1];
	UInt ligandChainID;
	sscanf(argv[2], "%d", &ligandChainID);
	PDBInterface* thePDB = new PDBInterface(infile);
	ensemble* theEnsemble = thePDB->getEnsemblePointer();
	molecule* pMol = theEnsemble->getMoleculePointer(0);
	protein* complex = static_cast<protein*>(pMol);
	protein* ligand = new protein(*complex);
	protein* apo = new protein(*complex);

	residue::setElectroSolvationScaleFactor(1.0);
	residue::setHydroSolvationScaleFactor(1.0);
	residue::setPolarizableElec(true);
	amberElec::setScaleFactor(1.0);
	amberVDW::setScaleFactor(1.0);
	residue::setTemperature(300);

	// delete every chain but ligand chain
	for (UInt i = 0; i < complex->getNumChains(); i++)
	{
		if (i != ligandChainID)
		{
			ligand->removeChain(i);
		}
	}
	cout << "ligand " << endl;

	// delete ligand chain
	//apo->removeChain(ligandChainID);
	cout << "apo" << endl;

	double complexE = complex->protEnergy();
	cout << "CE" << endl;
	double ligandE = ligand->protEnergy();
	cout << "LE" << endl;
	double apoE = apo->protEnergy();
	cout << "AE" << endl;
	double bindingEnergy = complexE-(ligandE+apoE);
	 cout << infile << " " << bindingEnergy << endl;
	return 0;
}
