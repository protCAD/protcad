#include <iostream>
#include <string>
#include "PDBInterface.h"
#include "ruler.h"


int main(int argc, char* argv[])
{	
	if (argc != 2)
	{	cout << "Usage: PDBInterfacetest pdbname" << endl;
		exit (1);
	}	

	string fname = argv[1];
	PDBInterface* thePDB = new PDBInterface(fname);
	ensemble* theEnsemble = thePDB->getEnsemble();

	molecule* pTheMolecule = theEnsemble->getMoleculePointer(0);
	if (pTheMolecule == 0)
	{       cout << "Failed" << endl;
		return 1;
	}
	protein* pTheProtein = static_cast<protein*>(pTheMolecule);
	protein* pNewProtein = new protein(*pTheProtein);

	axis myAxis = X_axis;
	pNewProtein->rotate(myAxis, 45.0);
	pNewProtein->translate(0.0,10.0,1.0);

	ruler* rmsd = new ruler();

	rmsd->setStationaryMolecule(pTheProtein);
	rmsd->setMobileMolecule(pNewProtein);

	rmsd->appendToList(0,20,40,0,"CA");
	rmsd->appendToList(1,20,40,0,"CA");

	rmsd->superimposeProteins();
	cout << "rmsd: " << rmsd->getrmsd() << endl;
	dblVec rotaxis = rmsd->getAxisOfRotation();
	cout << "axis of rotation: "<< rotaxis << endl;

	return 0;
}
