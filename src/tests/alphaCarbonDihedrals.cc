#include "typedef.h"
#include "ensemble.h"
#include "PDBInterface.h"

int main (int argc, char* argv[])
{
	enum aminoAcid {A,R,N,D,C,Q,E,G,H,I,L,K,M,F,P,S,T,W,Y,V};

	string inFile = argv[1];

	double step;
	sscanf(argv[2], "%lf", &step);

	//Input file should be four amino acids

	PDBInterface* thePDB = new PDBInterface(inFile);
	ensemble* theEnsemble = thePDB->getEnsemblePointer();
	molecule* theMol = theEnsemble->getMoleculePointer(0);
	protein* prot = static_cast<protein*>(theMol);


			/*dblVec CA1 = prot->getCoords(0,0,"CA");
			dblVec CA2 = prot->getCoords(0,1,"CA");
			dblVec CA3 = prot->getCoords(0,2,"CA");
			dblVec CA4 = prot->getCoords(0,3,"CA");

			cout << CMath::dihedral(CA1,CA2,CA3,CA4) << endl;*/


	double startTorsion = -180.0;
	double psi0,psi1,phi2,phi1;

	for (psi0 = startTorsion; psi0 < 180.0; psi0 += step)
	{ 
		for (phi1 = startTorsion; phi1 < 180.0; phi1 += step)
		{ 
			prot->setPsi(0,0,psi0);
			prot->setPhi(0,1,phi1*-1); 
			prot->setPsi(0,1,psi0*-1);
			prot->setPhi(0,2,phi1); 
			prot->setPsi(0,2,psi0);
			prot->setPhi(0,3,phi1*-1); 

			dblVec CA1 = prot->getCoords(0,0,"CA");
			dblVec CA2 = prot->getCoords(0,1,"CA");
			dblVec CA3 = prot->getCoords(0,2,"CA");
			dblVec CA4 = prot->getCoords(0,3,"CA");

			double angle = CMath::dihedral(CA1,CA2,CA3,CA4);

			cout <<  phi1 << " " << psi0 << " " << angle << endl;
 		} 
	}
 

	return 0;
}

