#include <iostream>
#include <string>
#include <vector>
#include "ensemble.h"
#include "PDBInterface.h"
#include "ran.h"
#include "generalio.h"

int main (int argc, char* argv[])
{
    residue::setCutoffDistance(6.0);
    pmf::setScaleFactor(0.0);
    rotamer::setScaleFactor(0.0);
    microEnvironment::setScaleFactor(0.0);
    amberVDW::setScaleFactor(1.0);
    amberVDW::setRadiusScaleFactor(1.0);
    amberVDW::setLinearRepulsionDampeningOff();
    amberElec::setScaleFactor(0.0);
    solvation::setItsScaleFactor(0.0);

    string protFile = argv[1];
    PDBInterface* thePDB = new PDBInterface(protFile);
    ensemble* theEnsemble = thePDB->getEnsemblePointer();
    molecule* theMol = theEnsemble->getMoleculePointer(0);
    protein* prot = static_cast<protein*>(theMol);
	prot->silenceMessages();

	UInt thrpos;
	sscanf(argv[2], "%u", &thrpos);

	prot->makeAtomSilent(0,thrpos,6);
	for (double chi = -180.0; chi < 180.0 ; chi = chi + 5.0) 
	{
		prot->setChi(0,thrpos,0,0,chi);
		double energy = prot->getPositionEnergy(0, thrpos);
		cout << "OG " << chi << " " << energy << endl;
	}
	
	return 0;
}

