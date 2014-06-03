#include <iostream>
#include <string>
#include "PDBInterface.h"
#include "ruler.h"
#include "lineSegment.h"

void getHelicalAxis(protein* _pTheProtein, lineSegment& _theAxis, dblVec& _centroid,
		dblVec& _axisVector, vector<UInt> helixlimits)
{
	protein* pNewProtein = new protein(*_pTheProtein);

	ruler* rmsd = new ruler();

	rmsd->setStationaryMolecule(_pTheProtein);
	rmsd->setMobileMolecule(pNewProtein);

	rmsd->appendToList(0,helixlimits[0],helixlimits[1]-1,0,"CA");
	rmsd->appendToList(1,helixlimits[0]+1,helixlimits[1],0,"CA");

	rmsd->superimposeProteins();
	_centroid = rmsd->getCentroidOfStationary();
	_axisVector = rmsd->getAxisOfRotation();
	return;
}

void calculateHelixAxes(protein* _pTheProtein, vector<lineSegment> & _theAxes, vector<dblVec> & _theCentroids, vector<dblVec> theAxisVectors)
{
	vector<vector<UInt> > helixLimits;

	vector<UInt> tempVec;
	tempVec.resize(2);
	tempVec[0] = 34;
	tempVec[1] = 65;
	helixLimits.push_back(tempVec);
	tempVec[0] = 71;
	tempVec[1] = 100;
	helixLimits.push_back(tempVec);
	tempVec[0] = 106;
	tempVec[1] = 141;
	helixLimits.push_back(tempVec);
	tempVec[0] = 151;
	tempVec[1] = 173;
	helixLimits.push_back(tempVec);
	tempVec[0] = 200;
	tempVec[1] = 226;
	helixLimits.push_back(tempVec);
	tempVec[0] = 246;
	tempVec[1] = 276;
	helixLimits.push_back(tempVec);
	tempVec[0] = 288;
	tempVec[1] = 310;
	helixLimits.push_back(tempVec);
/*
	Helix 1	:	34,65
	Helix 2 : 	71,100
	Helix 3 :	106,141
	Helix 4 : 	151,173
	Helix 5 :	200,226
	Helix 6 : 	246,276
	Helix 7 :	288,310
*/

	double null = 0.0;
	for (UInt i=0; i<helixLimits.size(); i++)
	{ 	lineSegment tempLine();
		dblVec tempCentroid(3,null);
		dblVec tempAxisVector(3,null);
		getHelicalAxis(_pTheProtein,tempAxis,tempCentroid,tempAxisVector,helixLimits[i]);
		_theAxes.push_back(tempAxisVector);
		_theCentroids.push_back(tempCentroid);
	}
	return;
}

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
	vector<lineSegment> theAxes;
	vector<dblVec> theAxisVectors;
	vector<dblVec> theCentroids;
	calculateHelixAxes(pTheProtein, theAxes, theCentroids, theAxisVectors);

	UInt null = 0;
	Matrix<UInt> contactMatrix(theCentroids.size(),theCentroids.size(),null);

/*
	for (UInt i=0; i<theCentroids.size(); i++)
	{
		for (UInt j=0; j<theCentroids.size(); j++)
		{
			if (CMath::distance(theCentroids[i],theCentroids[j])

	*/
	cout << "theCentroids.size() = " << theCentroids.size() << endl;

	return 0;
}
