#include "typedef.h"
#include "ensemble.h"
#include "PDBInterface.h"



int main(int argc, char* argv[])
{
	string inputFileName = argv[1];
	PDBInterface* thePDB = new PDBInterface(inputFileName);
	ensemble* theEnsemble = thePDB->getEnsemblePointer();
	molecule* theMol = theEnsemble->getMoleculePointer(0);
	protein* prot = static_cast<protein*>(theMol);

	UInt chainA, chainB;

	sscanf(argv[2], "%u", &chainA);
	sscanf(argv[3], "%u", &chainB);

	vector < vector <double> > CAA;
	vector < vector <double> > CAB;

	for (UInt i = 0; i < prot->getNumResidues(chainA); i ++)
	{
		dblVec coords = prot->getCoords(chainA, i, "CA");
		vector < double > CAcoords(0);
		CAcoords.push_back(coords[0]);
		CAcoords.push_back(coords[1]);
		CAcoords.push_back(coords[2]);

		CAA.push_back(CAcoords);
	}

	for (UInt i = 0; i < prot->getNumResidues(chainB); i ++)
	{
		dblVec coords = prot->getCoords(chainB, i, "CA");
		vector < double > CAcoords(0);
		CAcoords.push_back(coords[0]);
		CAcoords.push_back(coords[1]);
		CAcoords.push_back(coords[2]);

		CAB.push_back(CAcoords);
	}

	vector < vector <double> > Aaxis = HelixAxis(CAA);
	vector < vector <double> > Baxis = HelixAxis(CAB);

	vector < vector <double> > closestPoints(0);

        closestPoints = closest_points(Aaxis[0], Aaxis[1], Baxis[0], Baxis[1]);

cout << closestPoints[0][0] << " " << closestPoints[0][1] << " "<< closestPoints[0][2] << " --- "
                                << closestPoints[1][0] << " "<< closestPoints[1][1] << " "<< closestPoints[1][2] << endl;

        vector <double> BC = closestPoints[1] - closestPoints[0];

        double angle = (180.0/3.14159) * GetAngle(Aaxis[1], BC, Baxis[1]);

        double distance = sqrt(BC[0]*BC[0] + BC[1]*BC[1] + BC[2]*BC[2]);

        cout << argv[1] << " cross " << angle << " dist " << distance << endl;

	return 0;
}
