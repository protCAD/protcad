#include "typedef.h"
#include "ensemble.h"
#include "PDBInterface.h"
#include "ran.h"
#include <sstream>

int main (int argc, char* argv[])
{
	enum aminoAcid {A, R, N, D, C, Q, E, G, H, I, L, K, M, F, P, S, T, W, Y, V};

	string inputFileName = argv[1];

    PDBInterface* thePDB = new PDBInterface(inputFileName);
    ensemble* theEnsemble = thePDB->getEnsemblePointer();
    molecule* theMol = theEnsemble->getMoleculePointer(0);
    protein* prot = static_cast<protein*>(theMol);

    residue::setCutoffDistance(4.0);
    pmf::setScaleFactor(0.0);
    rotamer::setScaleFactor(0.0);
    microEnvironment::setScaleFactor(0.0);
    amberVDW::setScaleFactor(1.0);
    amberVDW::setRadiusScaleFactor(1.0);
    amberVDW::setLinearRepulsionDampeningOff();
    amberElec::setScaleFactor(0.0);
    solvation::setItsScaleFactor(0.0);

	prot->activateForRepacking(0,10);
	for (double phi = -180.0; phi < 180.0; phi += 5.0)
	{
		for (double psi = -180.0; psi < 180.0; psi += 5.0)
		{
			prot->setPhi(0,19, phi);
			prot->setPsi(0,19, psi);
			dblVec start = prot->getCoords(0,0,"CA");
			dblVec end = prot->getCoords(0,39, "CA");
			dblVec diff = start - end;

			double distance = sqrt(CMath::dotProduct(diff,diff));
		
			double energy = prot->intraEnergy();

			cout << "map " << phi << " " << psi << " " << energy << " " << distance <<  endl;
		}
	}
	return 0;
}

