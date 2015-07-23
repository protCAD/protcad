//*******************************************************************************************************
//*******************************************************************************************************
//*********************************                     *************************************************
//*********************************   mergeComplex 1.0  *************************************************
//*********************************                     *************************************************
//*******************************************************************************************************
//*******************************    -merge centroids-    ***********************************************
//*******************************************************************************************************

/////// Just specify a infile and outfile, it will diffuse to a generally effective minimum (totalsize).

//--Program setup----------------------------------------------------------------------------------------
#include <iostream>
#include <string>
#include <time.h>
#include <sstream>
#include "ensemble.h"
#include "PDBInterface.h"

dblVec carbonCentroid(protein* _prot, UInt _chain);

int main (int argc, char* argv[])
{
	//--Running parameters
	if (argc !=3)
	{
		cout << "mergeComplex <inFile.pdb> <outFile.pdb>" << endl;
		exit(1);
	}
	enum aminoAcid {A, R, N, D, Dh, C, Q, E, Eh, G, H, O, I, L, K, M, F, P, S, T, W, Y, V, dA, dR, dN, dD, dDh, dC, dQ, dE, dEh, dH, dO, dI, dL, dK, dM, dF, dP, dS, dT, dW, dY, dV};
	string infile = argv[1];
	PDBInterface* thePDB = new PDBInterface(infile);
	ensemble* theEnsemble = thePDB->getEnsemblePointer();
	molecule* pMol = theEnsemble->getMoleculePointer(0);
	protein* bundle = static_cast<protein*>(pMol);
	bundle->silenceMessages();
	residue::setCutoffDistance(10.0);
	rotamer::setScaleFactor(0.0);
	amberVDW::setScaleFactor(1.0);
	amberVDW::setRadiusScaleFactor(1.0);
	amberVDW::setLinearRepulsionDampeningOff();
	amberElec::setScaleFactor(1.0);
	srand (time(NULL));

	//--Initialize variables for loop
	double rotx, roty, rotz, transx, transy, transz, totaldistOld, totaldistNew;
	UInt test = 0, nobetter = 0;
	dblVec substrateCentroid = carbonCentroid(bundle, 0), ligandCentroid;
	string outFile;

	cout << endl << "\t*mergeComplex*" << endl;
	cout << endl << "Simulation running..." << endl;


//--Run optimizaiton loop till grand minimum-------------------------------------------------------------
	do
	{  
		//--Get centroids, distances and move values
		nobetter++, test++;
		ligandCentroid = carbonCentroid(bundle, 1);
		totaldistOld = CMath::distance(substrateCentroid, ligandCentroid);
		rotx = (rand() % 60)-30, roty = (rand() % 60)-30, rotz = (rand() % 60)-30;
		transx = (rand() % 3)-1, transy = (rand() % 3)-1, transz = (rand() % 3)-1;

		//--Make random motion
		bundle->rotateChain(1, X_axis, rotx);
		bundle->rotateChain(1, Y_axis, roty);
		bundle->rotateChain(1, Z_axis, rotz);
		bundle->translateChain(1, transx, 0, 0);
		bundle->translateChain(1, 0, transy, 0);
		bundle->translateChain(1, 0, 0, transz);
		
		//--Get new distance and test for merge
		ligandCentroid = carbonCentroid(bundle, 1);
		totaldistNew = CMath::distance(substrateCentroid, ligandCentroid);		
		if (totaldistNew < totaldistOld)
		{	
			test = 0, totaldistOld = totaldistNew;
			cout << totaldistNew << endl;
		}

		//--Revert if no improvement
		if (test != 0)
		{
			bundle->translateChain(1, 0, 0, (transz * -1));
			bundle->translateChain(1, 0, (transy * -1), 0);
			bundle->translateChain(1, (transx * -1), 0, 0);
			bundle->rotateChain(1, Z_axis, (rotz * -1));
			bundle->rotateChain(1, Y_axis, (roty * -1));
			bundle->rotateChain(1, X_axis, (rotx * -1));
		}	
	} while (totaldistNew > .5);	      


//--Print final energy and write a pdb file--------------------------------------------------------------
	cout << endl << "Simulation is complete!" << endl << endl;
	outFile = argv[2];
	pdbWriter(bundle, outFile);
	return 0;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////
//////// functions //////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////

///// get average coordinates of all carbons in chain ///////////////////////////////////////////////////
dblVec carbonCentroid(protein* _prot, UInt _chain)
{
	//--initialize and clear variables
	double number = 0;
	UInt numRes, numAtoms;
	string atomType;
	dblVec coords, coordsSum(3), coordsAve(3);
	coordsSum[0] = 0.0, coordsSum[1] = 0.0, coordsSum[2] = 0.0;
	coordsAve[0] = 0.0, coordsAve[1] = 0.0, coordsAve[2] = 0.0;
	
	//--loop through all atoms of all residues in search of carbon
	numRes = _prot->getNumResidues(_chain);
	for (UInt i = 0; i < numRes; i++)
	{
		numAtoms = _prot->getNumAtoms(_chain, i);
		for (UInt j = 0; j < numAtoms; j++)
		{
			atomType = _prot->getTypeStringFromAtomNum(_chain, i, j);
			if (atomType == "C")
			{
				number++;
				coords = _prot->getCoords(_chain, i, j);
				coordsSum[0] += coords[0];
				coordsSum[1] += coords[1];
				coordsSum[2] += coords[2];
			}
		}
	}

	//--get average of all carbon coordinates
	coordsAve[0] = ((coordsSum[0])/number);
	coordsAve[1] = ((coordsSum[1])/number);
	coordsAve[2] = ((coordsSum[2])/number);

	return coordsAve;
}



