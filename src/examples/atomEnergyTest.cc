#include <iostream>
#include <string>
#include <vector>
#include "ensemble.h"
#include "PDBInterface.h"
#include "ran.h"
#include "generalio.h"

int main (int argc, char* argv[])
{
	string inputFileName = argv[1];

	PDBInterface* thePDB = new PDBInterface(inputFileName);
	ensemble* theEnsemble = thePDB->getEnsemblePointer();
	molecule* theMol = theEnsemble->getMoleculePointer(0);
	protein* prot = static_cast<protein*>(theMol);

	residue::setCutoffDistance(12.0);
	pmf::setScaleFactor(0.0);
	rotamer::setScaleFactor(0.0);
	microEnvironment::setScaleFactor(0.0);
	amberVDW::setScaleFactor(1.0);
	amberVDW::setRadiusScaleFactor(1.0);
	amberVDW::setLinearRepulsionDampeningOff();
	amberElec::setScaleFactor(1.0);
	amberElec::setDielectricConstant(4.0);
	amberElec::distanceDependance = false;
	amberElec::setHighEnergyCutOff(false);
	solvation::setItsScaleFactor(0.0);

	UInt chain1, res1, atom1, chain2, res2, atom2;

	sscanf(argv[2], "%d", &chain1);
	sscanf(argv[3], "%d", &res1);
	sscanf(argv[4], "%d", &atom1);
	sscanf(argv[5], "%d", &chain2);
	sscanf(argv[6], "%d", &res2);
	sscanf(argv[7], "%d", &atom2);

	cout << "intraenergy:  " << prot->getIntraEnergy(chain1, res1, atom1, chain2, res2, atom2) << endl;

	return 0;
}
