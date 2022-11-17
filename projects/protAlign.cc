//*******************************************************************************************************
//*******************************************************************************************************
//***********************************                          ******************************************
//***********************************         protAlign        ******************************************
//***********************************                          ******************************************
//*******************************************************************************************************
//******  -Calculates best fit and RMSD of two pdbs and aligns the second protein-  *********************
//*******************************************************************************************************


#include "PDBInterface.h"
#include "ensemble.h"

int main (int argc, char* argv[])
{
	if (argc > 3 || argc < 2)
	{   cout << "protAlign <inFile1.pdb> <inFile2.pdb> - aligns two pdbs" << endl;
		cout << "protAlign <inFile1.pdb> - aligns protein to z-axis" << endl;
		exit(1); }
	string infile1 = argv[1];
	PDBInterface* thePDB1 = new PDBInterface(infile1);
	ensemble* theEnsemble1 = thePDB1->getEnsemblePointer();
	molecule* pMol1 = theEnsemble1->getMoleculePointer(0);
	protein* _prot1 = static_cast<protein*>(pMol1);
	
	if (argc == 2){
		_prot1->alignToAxis(Z_axis);
		pdbWriter(_prot1,infile1);
		return 0;
	}
	else{
		string infile2 = argv[2];
		PDBInterface* thePDB2 = new PDBInterface(infile2);
		ensemble* theEnsemble2 = thePDB2->getEnsemblePointer();
		molecule* pMol2 = theEnsemble2->getMoleculePointer(0);
		protein* _prot2 = static_cast<protein*>(pMol2);

		vector<dblVec> coord1;
		vector<dblVec> coord2;
		bool first = true;
		
		// Load backbone atoms into vector for fit and alignment
		for (UInt h = 0; h < _prot1->getNumChains(); h++)
		{
			for (UInt i = 0; i < _prot1->getNumResidues(h); i++)
			{
				coord1.push_back(_prot1->getCoords(h,i,0));
			}
		}
		for (UInt h = 0; h < _prot2->getNumChains(); h++)
		{
			for (UInt i = 0; i < _prot2->getNumResidues(h); i++)
			{
				coord2.push_back(_prot2->getCoords(h,i,0));
			}
		}
		int diff = 0;
		if(coord1.size() != coord2.size()){
			if (coord2.size() < coord1.size()){ diff = coord1.size()-coord2.size(); first = true;}
			else{diff = coord2.size()-coord1.size(); first = false;}
		}
		int maxsize, size;
		if (first){maxsize = coord2.size();}else{maxsize = coord1.size();}
		double rotmat[9]; double centroid1[3]; double centroid2[3]; double rmsd; double rmsdat[maxsize]; int ierr = 0;
		int list1[maxsize]; int list2[maxsize]; double bestRotMat[9]; double bestRmsd = 1E10;
		double newCoord1[maxsize*3]; double newCoord2[maxsize*3]; double newCoord3[maxsize*3];
		double bestcent1[3]; double bestcent2[3];
		for (int j = 0; j < maxsize; j++){rmsdat[j]=0;}
		double cutoff; int slide = 1;
		
		if (diff > 0){slide += diff;}
		for (int g = 0; g < slide; g++)
		{
			for (int h = 0; h < 5; h++)
			{
				size = 0;
				for (int i=0; i<maxsize; i++)
				{	
					if(rmsdat[i] < cutoff){
						for (int j=0; j<3; j++)
						{
							if (first) {
								newCoord1[ (size*3) + j] = coord1[i+g][j];
								newCoord2[ (size*3) + j] = coord2[i][j];
							}
							else{
								newCoord1[ (size*3) + j] = coord1[i][j];
								newCoord2[ (size*3) + j] = coord2[i+g][j];
							}
						}
						list1[size] = size+1;
						list2[size] = size+1;
						size++;
					}
				}

				// Calculate best fit of backbone atoms, rotation matrix and rmsd using fortran algorithm based on Machlachlan
				bestfit_(newCoord1, &size, newCoord2, &size, &size, newCoord3, list1, list2, &rmsd, &ierr, rotmat, centroid1, centroid2, rmsdat);
				if (rmsd < bestRmsd){
					bestRmsd = rmsd;
					for (int j = 0; j < 9; j++){bestRotMat[j]=rotmat[j];}
					for (int j = 0; j < 3; j++){bestcent1[j]=centroid1[j]; bestcent2[j]=centroid2[j];}
				}

				// Calculate two stdevs away from deviation mean for atom removal in next alignment cycle (pymol super method)
				double mean, sum=0.0, variance=0.0, stdev;
				for (int i = 0; i < size; i++){sum += rmsdat[i];}
				mean = sum/size;
				for(int i = 0; i < size; ++i)
				{
					variance += pow(rmsdat[i] - mean, 2);
					variance=variance/size;
					stdev = sqrt(variance);
				}
				cutoff = mean+(stdev*2);
			}
		}
		
		
		// Load rotation vector into rotation matrix
		dblMat rotMat(3,3,3);
		for (UInt i=0; i<3; i++)
		{	for (UInt j=0; j<3; j++)
			{
				rotMat[i][j] = bestRotMat[(j*3) + i];
			}
		}
		
		// Transform second or smaller pdb onto first using rotation matrix and centroid translation
		for (UInt i = 0; i < _prot2->getNumChains(); i++)
		{
			_prot2->translateChain(i, -bestcent2[0], -bestcent2[1], -bestcent2[2]);
			_prot2->transform(i,rotMat);
			_prot2->translateChain(i,bestcent1[0],bestcent1[1],bestcent1[2]);
		}
		pdbWriter(_prot2,infile2);
		cout << "RMSD: " << bestRmsd << endl;
	}
	return 0;
}
