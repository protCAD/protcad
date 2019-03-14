#include "typedef.h"
#include "ensemble.h"
#include "PDBInterface.h"
#include "ran.h"
#include <sstream>
int main (int argc, char* argv[])
{

        enum aminoAcid {A, R, N, D, C, Q, E, G, H, I, L, K, M, F, P, S, T, W, Y, V};


	string inputFileName = argv[1];

       // read in prot structure
        PDBInterface* thePDB = new PDBInterface(inputFileName);
        ensemble* theEnsemble = thePDB->getEnsemblePointer();
        molecule* theMol = theEnsemble->getMoleculePointer(0);
        protein* prot = static_cast<protein*>(theMol);


       residue::setCutoffDistance(8.0);
        pmf::setScaleFactor(0.0);
        rotamer::setScaleFactor(0.0);
        microEnvironment::setScaleFactor(0.0);
        amberVDW::setScaleFactor(1.0);
        amberVDW::setRadiusScaleFactor(0.9);
        amberVDW::setLinearRepulsionDampeningOff();
        amberElec::setScaleFactor(0.0);
        solvation::setItsScaleFactor(0.0);

	UInt res;
	sscanf(argv[2], "%u", &res);

	prot->activateAllForRepacking(0);

	for (UInt i = 0; i < prot->getNumResidues(0); i ++)
	{
		prot->setPhi(0,i, -180.0);
		prot->setPsi(0,i, -180.0);
		prot->mutate(0,i, res);
	}
	for (double psi = -180.0; psi < 180.0; psi = psi + 5.0)
	{
		cout << "output: ";
		for (double phi = -180.0; phi < 180.0; phi = phi +  5.0)
		{
			prot->setPhi(0,1,phi);
			prot->setPsi(0,1,psi);
		
			cout << prot->getNumHardClashes() << " ";
		}
		cout << endl;
	}
	pdbWriter (prot, "output.pdb");

	return 0;
}
