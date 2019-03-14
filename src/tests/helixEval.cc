#include "typedef.h"
#include "ensemble.h"
#include "PDBInterface.h"
#include "ran.h"
#include <sstream>
#include "/home/vikas/math.h"

int main (int argc, char* argv[])
{
	string inputFileName = argv[1];
	PDBInterface* thePDB = new PDBInterface(inputFileName);
	ensemble* theEnsemble = thePDB->getEnsemblePointer();
	molecule* theMol = theEnsemble->getMoleculePointer(0);
	protein* prot = static_cast<protein*>(theMol);

	UInt start = 5;
	UInt end = 5;

	double midPhi = prot->getPhi(0, start);
	double midPsi = prot->getPsi(0, start);

	int handedness = 0;  // -1 left, +1 right
	if (midPhi >= -120 && midPhi <= 0 && midPsi >= -120 && midPsi <= 0) handedness = 1;
	if (midPhi <= 120 && midPhi >= 0 && midPsi <= 120 && midPsi >= 0) handedness = -1;

	if (handedness == 0)
	{
		cout << "middle is not helical - aborting" << endl;
		exit(1);
	}
	
	for (UInt i = start - 1; i > 0; i --)
	{
		double phi = prot->getPhi(0, i);
		double psi = prot->getPsi(0, i);
		int itsHandedness = 0;
		if (phi >= -120 && phi <= 0 && psi >= -120 && psi <= 0) itsHandedness = 1;
		if (phi <= 120 && phi >= 0 && psi <= 120 && psi >= 0) itsHandedness = -1;
		if (itsHandedness == handedness && start == i + 1) start = i;
	}

	for (UInt i = end + 1; i < prot->getNumResidues(0); i ++)
	{
		double phi = prot->getPhi(0, i);
		double psi = prot->getPsi(0, i);
		int itsHandedness = 0;
		if (phi >= -120 && phi <= 0 && psi >= -120 && psi <= 0) itsHandedness = 1;
		if (phi <= 120 && phi >= 0 && psi <= 120 && psi >= 0) itsHandedness = -1;
		if (itsHandedness == handedness && end == i -1 ) end = i;
	}

	vector < vector < double > > CA(0);
	for (UInt i = 0; i < (end - start + 1); i ++)
	{
		dblVec coords =  prot->getCoords(0, start + i, "CA");
		vector < double > CAcoords(0);
		CAcoords.push_back(coords[0]);
		CAcoords.push_back(coords[1]);
		CAcoords.push_back(coords[2]);
		
		CA.push_back(CAcoords);
	}
	cout << "coords read in " << endl;
	
	vector < vector <double> > helixAxis = HelixAxis(CA); // first arg is point on line, second is direction

	cout << helixAxis[0][0] << " " << helixAxis[0][1] << " "<< helixAxis[0][2] << " --- " 
		<< helixAxis[1][0] << " "<< helixAxis[1][1] << " "<< helixAxis[1][2] << endl;

	CA.resize(0);	
	for (UInt i = 3; i < 15; i ++)
	{
		dblVec coords =  prot->getCoords(1, i, "CA");
		vector < double > CAcoords(0);
		CAcoords.push_back(coords[0]);
		CAcoords.push_back(coords[1]);
		CAcoords.push_back(coords[2]);
		
		CA.push_back(CAcoords);
	}

	vector < vector < double > > zAxis = HelixAxis(CA);
	
	vector < vector <double> > closestPoints(0);

	closestPoints = closest_points(zAxis[0], zAxis[1], helixAxis[0], helixAxis[1]);
	cout << closestPoints[0][0] << " " << closestPoints[0][1] << " "<< closestPoints[0][2] << " --- "
		                << closestPoints[1][0] << " "<< closestPoints[1][1] << " "<< closestPoints[1][2] << endl;
	
	
	vector <double> BC = closestPoints[1] - closestPoints[0];

	double angle = GetAngle(helixAxis[1], BC, zAxis[1]);

	double distance = sqrt(BC[0]*BC[0] + BC[1]*BC[1] + BC[2]*BC[2]);

	cout << argv[1] << " cross " << angle << " dist " << distance << endl;



	/*
	cout << handedness << " Helix starts at " << start << " and ends at " << end << endl;

	dblVec helixStart(3);
	dblVec helixEnd(3);

	dblVec n1 = prot->getCoords(0,start,"CA");
	dblVec n2 = prot->getCoords(0,start+1,"CA");
	dblVec n3 = prot->getCoords(0,start+2,"CA");
	dblVec n4 = prot->getCoords(0,start+3,"CA");

	dblVec c1 = prot->getCoords(0,end,"CA");
	dblVec c2 = prot->getCoords(0,end-1,"CA");
	dblVec c3 = prot->getCoords(0,end-2,"CA");
	dblVec c4 = prot->getCoords(0,end-3,"CA");

	helixStart = (n1 + n2 + n3 + n4) / 4.0;
	helixEnd = (c1 + c2 + c3 + c4) / 4.0;

	dblVec helixMid = (helixStart + helixEnd) / 2.0;

	double zmid = helixMid[2];

	helixMid[2] = 0.0;
	helixStart[2] = helixStart[2] - zmid;

	dblVec zax(3);
	zax[0] = 0.0;
	zax[1] = 0.0;
	zax[2] = 1.0;

	dblVec org(3);
	org[0] = 0.0;
	org[1] = 0.0;
	org[2] = 0.0;

	dblVec BA = org - zax;
	dblVec BC = org - helixMid;

	dblVec CB = helixMid - org;
	dblVec CD = helixMid - helixStart;

	dblVec P = CMath::cross(BA,BC);
	dblVec Q = CMath::cross(CB,CD);

	double sign = CMath::dotProduct(BA,Q) / (sqrt(CMath::dotProduct(BA,BA))*sqrt(CMath::dotProduct(Q,Q)));
	if (sign >= 0.0) sign = -1.0;
	else sign = 1.0;
	
	double angle = sign*(180.0 / 3.14159) * acos(CMath::dotProduct(P,Q) / (sqrt(CMath::dotProduct(P,P))*sqrt(CMath::dotProduct(Q,Q))));
	cout << "cross " << argv[1] << " "  << angle << endl;
	*/
	return 0;
}

