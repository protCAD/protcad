#include "typedef.h"
#include "ensemble.h"
#include "PDBInterface.h"
#include "ran.h"
#include <sstream>
#include<cstdlib>
#include<ctime>



void printTable(protein* _prot);

int main (int argc, char* argv[])
{

        enum aminoAcid {A, R, N, D, C, Q, E, G, H, I, L, K, M, F, P, S, T, W, Y, V, X};

	string inputFileName = argv[1];

       	// read in prot structure
        PDBInterface* thePDB = new PDBInterface(inputFileName);
        ensemble* theEnsemble = thePDB->getEnsemblePointer();
        molecule* theMol = theEnsemble->getMoleculePointer(0);
        protein* prot = static_cast<protein*>(theMol);

	// set energy function parameters - key for this program is radius scale, set to 0.9
        residue::setCutoffDistance(8.0);
        pmf::setScaleFactor(0.0);
        rotamer::setScaleFactor(0.0);
        microEnvironment::setScaleFactor(0.0);
        amberVDW::setScaleFactor(1.0);
        amberVDW::setRadiusScaleFactor(1.0);
        amberVDW::setLinearRepulsionDampeningOff();
        amberElec::setScaleFactor(0.0);
        solvation::setItsScaleFactor(0.0);

	prot->activateAllForRepacking(0);

	

	UInt clashes = prot->getNumHardClashes();
	double energy = prot->intraEnergy();

	for (UInt n = 0; n < 3; n ++)
	{

		for (UInt i = 0; i < prot->getNumResidues(0); i ++)
		{
			prot->setPhi(0,i, 180.0);
			prot->setPsi(0,i, 180.0);
			prot->mutate(0,i, 0);
			prot->setPhi(0,i, -65.0);
			prot->setPsi(0,i, -42.0);
	
		}
	
		clashes = prot->getNumHardClashes();
		energy = prot->intraEnergy();

		if (n == 0) pdbWriter(prot, "right1.pdb");

		//printTable(prot);	
	
		cout << "HELIX RIGHT " << n << " - C " << clashes << " E " << energy << endl;

		for (UInt i = 0; i < prot->getNumResidues(0); i ++)
		{
			prot->setPhi(0,i, 65.0);
			prot->setPsi(0,i, 42.0);
	
		}

		clashes = prot->getNumHardClashes();
		energy = prot->intraEnergy();

		if ( n == 0) pdbWriter(prot, "left1.pdb");

		//printTable(prot);		

		cout << "HELIX LEFT " << n << " - C " << clashes << " E " << energy << endl;
	}


	return 0;
}

void printTable(protein* _prot)
{
	for (UInt i = 0; i < _prot->getNumResidues(0); i ++)
	{
		cout << i << " PHI " << _prot->getPhi(0,i) << " PSI " << _prot->getPsi(0,i) << endl;
	}
	return;
}

