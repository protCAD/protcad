//*******************************************************************************************************
//*******************************************************************************************************
//***********************************                          ******************************************
//***********************************         protAlign        ******************************************
//***********************************                          ******************************************
//*******************************************************************************************************
//******  -Calculates best fit and RMSD of two pdbs and aligns the second to the first-  ****************
//*******************************************************************************************************


#include "PDBInterface.h"
#include "ensemble.h"

int main (int argc, char* argv[])
{
	if (argc !=3)
	{   cout << "protAlign <inFile1.pdb> <inFile2.pdb>" << endl;
		exit(1); }

	string infile1 = argv[1];
	PDBInterface* thePDB1 = new PDBInterface(infile1);
	ensemble* theEnsemble1 = thePDB1->getEnsemblePointer();
	molecule* pMol1 = theEnsemble1->getMoleculePointer(0);
	protein* _prot1 = static_cast<protein*>(pMol1);
	
	string infile2 = argv[2];
	PDBInterface* thePDB2 = new PDBInterface(infile2);
	ensemble* theEnsemble2 = thePDB2->getEnsemblePointer();
	molecule* pMol2 = theEnsemble2->getMoleculePointer(0);
	protein* _prot2 = static_cast<protein*>(pMol2);

    vector<dblVec> coord1;
    vector<dblVec> coord2;
    atomIterator theIter1(static_cast<protein*>(_prot1));
    atomIterator theIter2(static_cast<protein*>(_prot2));
    atom* pAtom;
    
    // Load backbone atoms into vector for fit and alignment
    for (;!(theIter1.last());theIter1++)
    {
       pAtom = theIter1.getAtomPointer(); 
       if (  pAtom->getName() == "N" || pAtom->getName() == "CA" || pAtom->getName() == "C" || pAtom->getName() == "O" ){
            coord1.push_back(pAtom->getCoords());
       }
    }
    for (;!(theIter2.last());theIter2++)
    {
       pAtom = theIter2.getAtomPointer(); 
       if (  pAtom->getName() == "N" || pAtom->getName() == "CA" || pAtom->getName() == "C" || pAtom->getName() == "O" ){
            coord2.push_back(pAtom->getCoords());
       }
    }

	// Now unroll the coordinates into a single long vector
	int numatm1 = coord1.size(); int numatm2 = coord2.size();
	int list1[numatm1]; int list2[numatm1];
	int nat = numatm1;
	double newCoord1[numatm1*3]; double newCoord2[numatm1*3]; double newCoord3[numatm1*3];
	double rotmat[9]; double centroid1[3]; double centroid2[3]; double rmsd = 0; int ierr = 0;
	for (int i=0; i<numatm1; i++)
	{	for (UInt c1=0; c1<3;c1++)
		{	newCoord1[ (i*3) + c1] = coord1[i][c1];
			newCoord2[ (i*3) + c1] = coord2[i][c1];
		}
		list1[i] = i+1;
		list2[i] = i+1;
	}
	
	// Calculate best fit of backbone atoms, rotation matrix and rmsd using fortran algorithm based on Machlachlan
	bestfit_(newCoord1, &numatm1, newCoord2, &numatm2, &nat, newCoord3, list1, list2, &rmsd, &ierr, rotmat, centroid1, centroid2);
	
	// Load rotation vector into rotation matrix
	dblMat rotMat(3,3,3);
    for (UInt i=0; i<3; i++)
    {	for (UInt j=0; j<3; j++)
		{
			rotMat[i][j] = rotmat[i*3 + j];
		}
    }
    
    // Transform second pdb onto first using rotation matrix
    for (UInt i = 0; i < _prot2->getNumChains(); i++)
    {
		_prot2->transform(i,rotMat);
    }
	pdbWriter(_prot2,infile2);
    cout << "RMSD = " << rmsd << endl;
	return 0;
}
