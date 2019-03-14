#include "pdbReader.h"
#include "pdbWriter.h"
#include "atomIterator.h"
#include "CMath.h"
//#include "svmt.h"
#include "typedef.h"
#include <time.h>

int main()
{
		
	string fname1 = "1";
	string fname2 = "2";

	molecule* pTheMolecule1 = pdbReader(fname1);
	if (pTheMolecule1 == 0)
	{	return 1;
	}
	protein* pTheProtein1 = static_cast<protein*>(pTheMolecule1);


	molecule* pTheMolecule2 = pdbReader(fname2);
	if (pTheMolecule2 == 0)
	{	return 1;
	}
	protein* pTheProtein2 = static_cast<protein*>(pTheMolecule2);

	vector<double> theWeights;
	vector<dblVec> coord1;
	vector<dblVec> coord2;

	atom* pAtom;	
	atomIterator ai1(pTheProtein1);
	for ( ; !(ai1.last()); ai1++)
	{	
		pAtom = ai1.getAtomPointer();
		if (pAtom->getName() == "CA")
		{	coord1.push_back(pAtom->getCoords());
			theWeights.push_back(1.0);
		}
	}

	atomIterator ai2(pTheProtein2);
	for ( ; !(ai2.last()); ai2++)
	{	
		pAtom = ai2.getAtomPointer();
		if (pAtom->getName() == "CA")
			coord2.push_back(pAtom->getCoords());
	}
/*
	for (UInt i=0; i<coord2.size(); i++)
	{
		for (UInt j=0; j<3; j++)
		{	cout << coord2[i][j] << " ";
		}
		cout << endl;
	}
*/

	if (coord1.size() != coord2.size())
	{	cout << "error - not one to one correspondence " << endl;
	}
	cout << "size of coord1 = " << coord1.size() << endl;
	cout << "size of coord2 = " << coord2.size() << endl;
	cout << "size of weights  = " << theWeights.size() << endl;

	vector<dblVec>* pCoord1 = &coord1;
	vector<dblVec>* pCoord2 = &coord2;
	vector<double>* pTheWeights = &theWeights;

/*
	dblVec centroid1 = CMath::centroid(pCoord1,pTheWeights);
	cout << "returned from first centroid calculation" << endl;
	dblVec centroid2 = CMath::centroid(pCoord2,pTheWeights);
*/	
	dblVec centroid1(3);
	dblVec centroid2(3);
	dblMat rotmat(3,3);
	dblVec trnvec(3);

	CMath::fndmat(pTheWeights,pCoord1,pCoord2,centroid1,centroid2,
		rotmat,trnvec);

	cout << "Back in main" << endl;
	cout << "Centroid1:" << endl;
	cout << centroid1 << endl;
	cout << endl;
	cout << "Centroid2:" << endl;
	cout << centroid2 << endl;
	cout << endl;
	cout << "rotmat:" << endl;
	cout << rotmat << endl;
	cout << endl;
	cout << "trnvec:" << endl;
	cout << trnvec << endl;


	delete pTheProtein1;
	delete pTheProtein2;

	return 0;
}
