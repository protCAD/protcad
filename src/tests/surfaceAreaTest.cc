#include <iostream>
#include <string>
#include "ensemble.h"
#include "PDBInterface.h"
int main(int argc, char* argv[])
{
	string fileName = argv[1];
	cout << "reading in file " << fileName << endl;

	PDBInterface* thePDB = new PDBInterface(fileName);
	cout << "PDBinterface invoked!" << endl;
	ensemble* theEnsemble = thePDB->getEnsemble();
	cout << "ensemble created!" << endl;
	molecule* pTheMolecule = theEnsemble->getMoleculePointer(0);
	cout << "molecule created!" << endl;
	protein* pProt = static_cast<protein*>(pTheMolecule);
	cout << "protein created." << endl;
	if (pTheMolecule == 0)
	{
		cout << "failed!" << endl;
		exit(1);
	}

		pProt->initializeSpherePoints();
		pProt->removeSpherePoints();
		cout << "SURFACE AREA " << pProt->tabulateSurfaceArea() << endl;
		cout << "SOLVATION ENERGY " <<pProt->tabulateSolvationEnergy(0) << endl;

/*
	atom* atom1 = new atom;
	atom* atom2 = new atom;
	atom* atom3 = new atom;

	double radius;
	sscanf(argv[2], "%lf", &radius);

	atom1->setAtomicRadius(radius);
	atom2->setAtomicRadius(radius);
	atom3->setAtomicRadius(radius);
	atom1->setItsProbeRadius(0);
	atom2->setItsProbeRadius(0);
	atom3->setItsProbeRadius(0);
	atom3->setCoords(0,0,0);
//	atom1->generateNewSpherePoints();
//	atom2->generateNewSpherePoints();
//	sscanf(argv[3], "%lf", &coord);

	for (double coord = 5; coord >= 0; coord -= 0.05)
	{
		atom1->setCoords(0,0,-1* coord);
		atom2->setCoords(0,0, 1* coord);

		atom1->clearSpherePoints();
		atom2->clearSpherePoints();	
		atom3->clearSpherePoints();

		atom3->removeBlockedPoints(atom2);
		atom3->removeBlockedPoints(atom1);
		cout << coord << " " << atom3->calculateExposedSASA() << endl;
		
	}
*/	return 0;
}
