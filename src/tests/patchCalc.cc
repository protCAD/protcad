#include <iostream>
#include <string>
#include "ensemble.h"
#include "PDBInterface.h"

#define THRESHOLD 10
#define PI 3.14159
vector <UIntVec> populatePatch (protein* _prot, UIntVec _center, vector <UIntVec> _surfacePos);
int main (int argc, char* argv[])
{
 //   	if (argc != 4)
 //   	{
 //       	cout << "dimerBuilder  <infilename.pdb> <target.pdb>  <outfilename.pdb> " << endl;
 //   	exit(1);
 //   	}
	enum aminoAcid {A, R, N, D, C, Q, E, G, H, I, L, K, M, F, P, S, T, W, Y, V};

    	string infile = argv[1];
    	PDBInterface* thePDB = new PDBInterface(infile);
    	ensemble* theEnsemble = thePDB->getEnsemblePointer();
    	molecule* pMol = theEnsemble->getMoleculePointer(0);
    	protein* prot = static_cast<protein*>(pMol);

	string sasaBaselinefile = argv[2];

	ifstream inFile;
	inFile.open(sasaBaselinefile.c_str());

	vector < double > baseLineSASA;
	
	string currentLine;	
	vector <string> parsedStrings;

	while (getline(inFile, currentLine, '\n'))
	{
		double sasatemp;
		parsedStrings = Parse::parse(currentLine);
		sscanf(parsedStrings[1].c_str(), "%lf", &sasatemp);
		baseLineSASA.push_back(sasatemp);
	}

    	residue::setCutoffDistance(6.0);
    	pmf::setScaleFactor(0.0);
    	rotamer::setScaleFactor(1.0);
    	microEnvironment::setScaleFactor(0.0);
    	amberVDW::setScaleFactor(1.0);
    	amberVDW::setRadiusScaleFactor(1.0);
    	amberVDW::setLinearRepulsionDampeningOn();
    	amberElec::setScaleFactor(0.0);
    	solvation::setItsScaleFactor(0.0);

	prot->initializeSpherePoints();
	prot->removeSpherePoints();

	vector <UIntVec> surfacePos;
	vector <double> surfacePosSASA;

	for (UInt i = 0; i < prot->getNumChains(); i ++)
	{
		for (UInt j = 0; j < prot->getNumResidues(i); j ++)
		{
			UInt type = prot->getTypeFromResNum(i,j);
			double sasaratio = prot->tabulateSurfaceArea(i,j) / baseLineSASA[type];
			cout << "chain " << i << " res " << j << " " << sasaratio << endl;
			if (sasaratio > 0.2)
			{
				UIntVec tempPos;
				tempPos.push_back(i);
				tempPos.push_back(j);
				surfacePos.push_back(tempPos);
				surfacePosSASA.push_back(sasaratio);
			}
		}
	}

	vector < vector <UIntVec> >  patches;

	for (UInt i = 0; i < surfacePos.size(); i ++)
	{
		if (prot->getTypeFromResNum(surfacePos[i][0],surfacePos[i][1]) == K)
		{
			vector <UIntVec> tempPatch = populatePatch(prot, surfacePos[i], surfacePos);
			patches.push_back(tempPatch);
		}
	}

	for (UInt i = 0; i < patches.size(); i ++)
	{
		cout << "patch " << i << " central residue is " << prot->getTypeStringFromResNum(patches[i][0][0], patches[i][0][1]) << " " <<  patches[i][0][1] << endl;
		for (UInt j = 1; j < patches[i].size(); j ++)
		{
			cout << "\t" << prot->getTypeStringFromResNum(patches[i][j][0], patches[i][j][1]) << " " << patches[i][j][1] << endl;
		}
	}

	return 0;
}




vector <UIntVec> populatePatch (protein* _prot, UIntVec _center, vector <UIntVec> _surfacePos)
{
	dblVec centralResCoords = _prot->getCoords(_center[0], _center[1], "CA");

	double thresholdsq = THRESHOLD*THRESHOLD;

	vector <UIntVec> nearNeighbors;

	for (UInt i = 0; i < _surfacePos.size(); i ++)
	{
		if (!(_center == _surfacePos[i]))
		{
			dblVec tempPosCoords = _prot->getCoords(_surfacePos[i][0], _surfacePos[i][1], "CA");
			dblVec diff = centralResCoords - tempPosCoords;
			double distsq = CMath::dotProduct(diff,diff);
			if (distsq <= thresholdsq)
			{
				nearNeighbors.push_back(_surfacePos[i]);
			}
		}
	}

	cout << "nn size " << nearNeighbors.size() << endl;

	dblVec centerOfMassCoords = centralResCoords;
	for (UInt i = 0; i < nearNeighbors.size(); i ++)
	{
		dblVec tempPosCoords = _prot->getCoords(nearNeighbors[i][0], nearNeighbors[i][1], "CA");
		centerOfMassCoords = centerOfMassCoords + tempPosCoords;
	}
	
	centerOfMassCoords = centerOfMassCoords / (double)(nearNeighbors.size() + 1);

	dblVec centerVec = centralResCoords - centerOfMassCoords;

	vector <UIntVec> patch;

	patch.push_back(_center);  // the central amino acid is the first in the patch

	for (UInt i = 0; i < nearNeighbors.size(); i ++)
	{
		dblVec tempPosCoords = _prot->getCoords(nearNeighbors[i][0], nearNeighbors[i][1], "CA");
		dblVec localVec = tempPosCoords - centerOfMassCoords;

		double angle = acos( CMath::dotProduct(centerVec,localVec) / sqrt(CMath::dotProduct(centerVec,centerVec) * CMath::dotProduct(localVec,localVec)));

		if (angle * 180/PI <= 110)
		patch.push_back(nearNeighbors[i]);
	}

	return patch;	

}

