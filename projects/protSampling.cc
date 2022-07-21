//*******************************************************************************************************
//*******************************************************************************************************
//*****************************                         *************************************************
//*****************************       protSampling      *************************************************
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
		cout << "protSampling <infile.pdb> <outfile.pdb>" << endl;
		exit(1);
	}
	random_device rd; srand((int)rd());
	string infile = argv[1];
	string outfile = argv[2];

	amberElec::setScaleFactor(1.0);
	amberVDW::setScaleFactor(1.0);
	residue::setElectroSolvationScaleFactor(1.0);
	residue::setHydroSolvationScaleFactor(1.0);
	residue::setPolarizableElec(true);
	residue::setEntropyFactor(1.0);
	residue::setTemperature(300);

	PDBInterface* thePDB = new PDBInterface(infile);
	ensemble* theEnsemble = thePDB->getEnsemblePointer();
	molecule* pMol = theEnsemble->getMoleculePointer(0);
	protein* prot = static_cast<protein*>(pMol);
	pdbWriter(prot, outfile);
	double oldEnergy = prot->protEnergy(), newEnergy = oldEnergy;

	for (UInt i = 0; i < 1000; i++)
	{
		delete thePDB;
		thePDB = new PDBInterface(outfile);
		theEnsemble = thePDB->getEnsemblePointer();
		pMol = theEnsemble->getMoleculePointer(0);
		prot = static_cast<protein*>(pMol);
		oldEnergy = prot->protEnergy();

		prot->protSampling(50);
		newEnergy = prot->protEnergy();
		
		PDBInterface* theOldPDB = new PDBInterface(outfile);
		ensemble* theOldEnsemble = theOldPDB->getEnsemblePointer();
		molecule* pMolOld = theOldEnsemble->getMoleculePointer(0);
		protein* protOld = static_cast<protein*>(pMolOld);
		oldEnergy = protOld->protEnergy();
		delete theOldPDB;

		if (newEnergy < oldEnergy)
		{
			pdbWriter(prot, outfile);
		}
	} 
	return 0;
}
