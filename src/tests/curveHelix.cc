#include <iostream>
#include <string>
#include "ensemble.h"
#include "PDBInterface.h"

#define PI 3.14159

int main (int argc, char* argv[])
{

	string infile = argv[1];
	string outfile = argv[2];
	
	// read in protein
	PDBInterface* thePDB = new PDBInterface(infile);
	ensemble* theEnsemble = thePDB->getEnsemblePointer();
	molecule* pMol = theEnsemble->getMoleculePointer(0);
	protein* prot = static_cast<protein*>(pMol);

	double radius;
	sscanf(argv[3], "%lf", &radius);

	cout << "radius " << radius << endl;
	for (UInt res = 0; res < prot->getNumResidues(0); res ++)
	{
		for (UInt atom = 0; atom < prot->getNumAtoms(0,res); atom ++)
		{
			dblVec coord = prot->getCoords(0, res, atom);
			double z = coord[2];
			double r = radius + coord[0]; // bend on x
			
			double theta = z/radius;
			double y = r * sin(theta);
			double x = r * cos(theta);

			cout << "old coords " << coord[0] << " " << coord[1] << " " << coord[2] << endl;

			dblVec newCoord(3);
			newCoord[0] = x;
			newCoord[1] = coord[1];
			newCoord[2] = y;

			cout << "new coords " << newCoord[0] << " " << newCoord[1] << " " << newCoord[2] << endl;
			
			prot->setCoords(0,res,atom,newCoord);
			cout << "atom " << atom << " res " << res << " set" << endl;
		}
	}
	
	pdbWriter(prot, outfile);

}



