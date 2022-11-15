#include "ruler.h"

ruler::ruler()
{
	double null = 0.0;
	itsCentroid1.newsize(3);
	itsCentroid2.newsize(3);
	itsRotationMatrix.newsize(3,3);
	for (UInt i=0; i<3; i++)
	{	itsCentroid1[i] = null;
		itsCentroid2[i] = null;
		for (UInt j=0; j<3;j++)
		{	itsRotationMatrix[i][j] = null;
		}
	}
	pItsStationaryMolecule = 0;
	pItsMobileMolecule = 0;
	itsrmsd = null;
}

ruler::~ruler()
{
}

void ruler::superimposeProteins ()
{
	vector<dblVec> coord1;
	vector<dblVec> coord2;

	protein* pProt1 = static_cast<protein*>(pItsStationaryMolecule);
	protein* pProt2 = static_cast<protein*>(pItsMobileMolecule);

	//cout << "pItsStationaryMolecule = " << pItsStationaryMolecule << endl;
	//cout << "pItsMobileMolecule = " << pItsMobileMolecule << endl;
	cout << endl;
	cout << "Atom equivalence list:" << endl;
	for (UInt i=0; i<atomList1.size(); i++)
	{
		cout << atomList1[i][0] <<" "<< atomList1[i][1]<<" "<< atomList1[1][2] << "     ";
		cout << atomList2[i][0] <<" "<< atomList2[i][1]<<" "<< atomList2[1][2];
		cout << endl;
		}

	atom* pAtom;
	dblVec pDV;

	if (atomList1.size() != atomList2.size())
	{
		cout << "Error! in superimpose(): different number";
		cout << " of atoms in two lists!" << endl;
		cout << "list1: " << atomList1.size() << "  ";
		cout << "list2: " << atomList2.size() << endl;
		return;
	}

	for (UInt i=0;i<atomList1.size();i++)
	{
		pAtom =((pProt1->itsChains[atomList1[i][0]])->itsResidues[atomList1[i][1]])->itsAtoms[atomList1[i][2]]; 
		if (pAtom)
			coord1.push_back(pAtom->getCoords());
	}

	for (UInt i=0;i<atomList2.size();i++)
	{
		pAtom =((pProt2->itsChains[atomList2[i][0]])->itsResidues[atomList2[i][1]])->itsAtoms[atomList2[i][2]]; 
		if (pAtom)
			coord2.push_back(pAtom->getCoords());
	}
	
	//cout << "Reached waypoint 1" << endl;

	if (coord1.size() > 4000)
	{
		cout << "Too many atoms to do superposition" << endl;
		cout << "you need to edit the source code" << endl;
	}

	// Now unroll the coordinates into a single long vector
	// of size 3*Natm1
	int numatm1 = coord1.size();
	int numatm2 = coord2.size();
	int list1[2000];
	int list2[2000];
	int nat = numatm1;
	double newCoord1[6000];
	double newCoord2[6000];
	double newCoord3[6000];
	//double newWeights[2000];
	double rotmat[9];
	double centroid1[3];
	double centroid2[3];
	double rmsdat[2000];
	double rmsd = 0;
	int ierr = 0;
	//cout << "Reached waypoint 2" << endl;
	for (int i=0; i<numatm1; i++)
	{	for (UInt c1=0; c1<3;c1++)
		{	newCoord1[ (i*3) + c1] = coord1[i][c1];
			newCoord2[ (i*3) + c1] = coord2[i][c1];
		}
		list1[i] = i+1;
		list2[i] = i+1;
		//newWeights[i] = itsWeights[i];
	}

	//cout << "Reached waypoint 4" << endl;
	bestfit_(newCoord1, &numatm1, newCoord2,
		&numatm2, &nat, newCoord3, list1, list2,
		&rmsd, &ierr, rotmat, centroid1, centroid2, rmsdat);

	//cout << "Reached waypoint 5" << endl;

	// now roll everything back up into a format which
	// is less like fortran and more like c++


	for (UInt i=0; i<3; i++)
	{	itsCentroid1[i] = centroid1[i];
	}
	for (UInt i=0; i<3; i++)
	{	itsCentroid2[i] =  centroid2[i];
	}
	for (UInt i=0; i<3; i++)
	{	for (UInt j=0; j<3; j++)	
			itsRotationMatrix[i][j] = rotmat[i*3 + j];
	}
	itsrmsd = rmsd;
	//cout << "Reached waypoint 6" << endl;
	return;
}

dblVec ruler::getAxisOfRotation()
{
	dblVec axis(3);
	axis[0] = itsRotationMatrix[2][1] - itsRotationMatrix[1][2];
	axis[1] = itsRotationMatrix[0][2] - itsRotationMatrix[2][0];
	axis[2] = itsRotationMatrix[1][0] - itsRotationMatrix[0][1];
	double norm = CMath::dotProduct(axis,axis);
	norm = sqrt(norm);
	axis = axis / norm;
	return axis;
}

void ruler::appendToList(UInt _listnum, UInt _startRes, UInt _endRes,
	UInt _chain, string _atomNames)
{
	if (_listnum != 0 && _listnum != 1)
	{
		cout << "Error from ruler::generateList" << endl;
		cout << "Cannot generate list - invalid listum: ";
		cout << _listnum << endl;
		return;
	}
	protein* pProt;
	vector<string> atomnames = parseString(_atomNames);
	vector<vector<UInt> > tempList;
	if (_listnum ==0)
	{
		pProt = static_cast<protein*>(pItsStationaryMolecule);
	}
	else
	{
		pProt = static_cast<protein*>(pItsMobileMolecule);
	}

	UInt startIndex = pProt->getIndexFromResNum(_chain,_startRes);
	UInt stopIndex = pProt->getIndexFromResNum(_chain,_endRes);
	UInt numnames = atomnames.size();
	if (numnames == 0)
	{
		cout << "Error from ruler::generateList" << endl;
		cout << "Cannot generate list - no atomnames in string: ";
		cout << _atomNames << endl;
		return;
	}
	for (UInt i=startIndex; i<=stopIndex; i++)
	{
		UInt numatm = pProt->itsChains[_chain]->itsResidues[i]->getNumAtoms();
		for (UInt j=0; j<numatm; j++)
		{	
			string tempName = pProt->itsChains[_chain]->itsResidues[i]->itsAtoms[j]->getName();
			for (UInt k=0; k<numnames; k++)
			{
				if (tempName == atomnames[k])
				{
					vector<UInt> tempVec;
					tempVec.push_back(_chain);
					tempVec.push_back(i);
					tempVec.push_back(j);
					tempList.push_back(tempVec);
					itsWeights.push_back(1.0);
				}
			}
		}
	}
	for (UInt i=0; i<tempList.size(); i++)
	{
		if (_listnum ==0)
		{
			atomList1.push_back(tempList[i]);
		}
		else
		{
			atomList2.push_back(tempList[i]);;
		}
	}
}
