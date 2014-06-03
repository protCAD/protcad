#include "pdbReader.h"
#include "pdbWriter.h"
#include "atomIterator.h"
#include "typedef.h"
#include <time.h>

int main()
{	time_t startTime,endTime,runTime;
	time(&startTime);

	string filename = "1";
	molecule* pTheProtein1 = pdbReader(filename);
	if (pTheProtein1 == 0)
	{	return 1;
	}

	string filename2 = "2";
	molecule* pTheProtein2 = pdbReader(filename2);
	if (pTheProtein2 == 0)
	{	return 1;
	}

/*
    dblVec theVec(3,4.6);
    point origin(0.0,0.0,0.0);
    pTheProtein2->rotate(0,origin,theVec,60.0);
    pTheProtein2->translate(0,10.0,5.0,-3.5);

    string outfilename = "3";
    pdbWriter(static_cast<protein*>(pTheProtein2),outfilename);
*/

    double tempWeight = 1.0;
    vector<dblVec> coord1;
    vector<dblVec> coord2;
    vector<double> weights;

    atomIterator theIter1(static_cast<protein*>(pTheProtein1));
    atomIterator theIter2(static_cast<protein*>(pTheProtein2));

    atom* pAtom;
    dblVec pDV;

    for (;!(theIter1.last());theIter1++)
    {
       pAtom = theIter1.getAtomPointer(); 
       if ( pAtom->getName() == "CA" )
       {
            coord1.push_back(pAtom->getCoords());
            weights.push_back(tempWeight);
       }
    }

    for (;!(theIter2.last());theIter2++)
    {
       pAtom = theIter2.getAtomPointer(); 
       if ( pAtom->getName() == "CA" )
       {
            coord2.push_back(pAtom->getCoords());
       }
    }

	// Now unroll the coordinates into a single long vector
	// of size 3*Natm1
	int numatm1 = coord1.size();
	int numatm2 = coord2.size();
	int list1[numatm1];
	int list2[numatm1];
	int nat = numatm1;
	double newCoord1[numatm1*3];
	double newCoord2[numatm1*3];
	double newCoord3[numatm1*3];
	double newWeights[numatm1];
	double rotmat[9];
	double centroid1[3];
	double centroid2[3];
	double trnvec[3];
	double rmsd = 0;
	int ierr = 0;
	for (UInt i=0; i<numatm1; i++)
	{	for (UInt c1=0; c1<3;c1++)
		{	newCoord1[ (i*3) + c1] = coord1[i][c1];
			newCoord2[ (i*3) + c1] = coord2[i][c1];
		}
		list1[i] = i+1;
		list2[i] = i+1;
		newWeights[i] = weights[i];
	}

	bestfit_(newCoord1, &numatm1, newCoord2,
		&numatm2, &nat, newCoord3, list1, list2,
		&rmsd, &ierr, rotmat, centroid1, centroid2);


/*
    cout << "size of pCoord1 = " << pCoord1->size() << endl;
    cout << "size of pCoord2 = " << pCoord2->size() << endl;
    cout << "size of pWeights = " << pWeights->size() << endl;


    dblVec centroid1(3,0.0);
    dblVec centroid2(3,0.0);
    dblMat rotmat(3,3,0.0);
    dblVec trnvec(3,0.0);
    CMath::fndmat(pWeights,pCoord1,pCoord2,centroid1,centroid2,
                  rotmat,trnvec);
    cout << " Centroid1 " << centroid1 << endl;
    cout << " Centroid2 " << centroid2 << endl;
    cout << " rotation matrix " << rotmat << endl;
    cout << " translation vector " << trnvec << endl;
    cout << "rmsd = " << CMath::rmsd(pWeights,pCoord1,pCoord2) << endl;
*/
    cout << " Centroid1 ";
    for (UInt i=0; i<3; i++)
    {	cout  << centroid1[i] << "  ";
    }
    cout << endl;
    cout << " Centroid2 ";
    for (UInt i=0; i<3; i++)
    {	cout  << centroid2[i] << "  ";
    }
    cout << endl;
    cout << " rotation matrix " << endl;
    for (UInt i=0; i<3; i++)
    {	for (UInt j=0; j<3; j++)	
		cout  << rotmat[i*3 + j] << "  ";
	cout << endl;
    }
    cout << endl;


    cout << endl;
    cout << "rmsd = " << rmsd << endl;


/*
	time(&endTime);
	string outfilename = "time.dat";
	ofstream oFile;
	oFile.open(outfilename.c_str());
	oFile << "startTime = " << startTime << "\n";
	oFile << "endTime = " << endTime << "\n";
	runTime = endTime - startTime;
	oFile << "runTime = " << runTime << "\n";
	oFile.close();
*/

	return 0;
}
