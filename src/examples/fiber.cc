#include <iostream>
#include <string>
#include "ensemble.h"
#include "PDBInterface.h"

void createFiber(protein* _prot, double _rad1, double _rad2, double _rot1);

int main (int argc, char* argv[])
{
	string infile = argv[1];
	PDBInterface* thePDB = new PDBInterface(infile);
	ensemble* theEnsemble = thePDB->getEnsemblePointer();
	molecule* pMol = theEnsemble->getMoleculePointer(0);
	protein* prot = static_cast<protein*>(pMol);
	
	double rad1, rad2, rot1;

	string rad1s = argv[2];
	string rad2s = argv[3];
	string rot1s = argv[4];
	
	sscanf(rad1s.c_str(), "%lf", &rad1);
	sscanf(rad2s.c_str(), "%lf", &rad2);
	sscanf(rot1s.c_str(), "%lf", &rot1);
	
	createFiber(prot, rad1, rad2, rot1);

	string filename = "fiberout_" + rad1s;  
	filename = filename + "_";
	filename = filename + rad2s;
	filename = filename + "_";
	filename = filename + rot1s;
	filename = filename + ".pdb";	
	pdbWriter(prot, filename);
}

void createFiber(protein* _prot, double _rad1, double _rad2, double _rot1)
{
	_prot->rotate(Y_axis, _rot1);
	_prot->translate(0.5 * _rad2, 0, 0);  // move a-f by half radius2
	_prot->translate(0, 0, _rad1, 0); // move a,d up
	_prot->translate(3, 0, _rad1, 0);
	_prot->translate(2, 0, -1.0 * _rad1, 0); // move c,f down
	_prot->translate(5, 0, -1.0 * _rad1, 0); 
	_prot->rotate(3, Z_axis, 180.0);
	_prot->rotate(4, Z_axis, 180.0);
	_prot->rotate(5, Z_axis, 180.0);

	return;
}
