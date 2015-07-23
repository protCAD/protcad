//*******************************************************************************************************
//*******************************************************************************************************
//*************************************                 *************************************************
//*************************************  foldingC 1.0   *************************************************
//*************************************                 *************************************************
//*******************************************************************************************************
//*********************************   -simulated folding-   *********************************************
//*******************************************************************************************************

/////// Just specify a infile and outfile, it will optimize to a generally effective minimum (totalsize).

//--Program setup----------------------------------------------------------------------------------------
#include <iostream>
#include <fstream>
#include <string>
#include <time.h>
#include <sstream>
#include "ensemble.h"
#include "PDBInterface.h"
int main (int argc, char* argv[])
{
	//--Running parameters
	if (argc !=2)
	{
		cout << "foldingC <inFile.pdb>" << endl;
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
	rotamer::setScaleFactor(1.0);
	amberVDW::setScaleFactor(1.0);
	amberVDW::setRadiusScaleFactor(1.0);
	amberVDW::setLinearRepulsionDampeningOff();
	amberElec::setScaleFactor(1.0);
	srand (time(NULL));

	//--Initialize variables for loop
	double deltaTheta, currentangle, sum;
	double goodEnergy, Energy, dist1, dist2, dist3, preArea, postArea;
	int direction, distance;
	int residues[] = {0,1,2,3,4,5,6,7};
	int resSize = sizeof(residues)/sizeof(residues[0]);
	dblVec coords1, coords2, coords3;
	UInt randchain, randres;
	UInt resNum, randtype, chainNum = bundle->getNumChains();
	string outPDB, countstr;
	stringstream convert;
	ofstream energyFile("Energy"); 
	ofstream areaFile("Area");
	delete thePDB;

//--Run optimizaiton loop till grand minimum-------------------------------------------------------------
	for (UInt a = 0; a<100; a++)
	{
		PDBInterface* thePDB = new PDBInterface(infile);
		ensemble* theEnsemble = thePDB->getEnsemblePointer();
		molecule* pMol = theEnsemble->getMoleculePointer(0);
		protein* bundle = static_cast<protein*>(pMol);
		UInt test = 0, open = 0, goodEnergy = 1.0;
		do
		{  
			//--Generate random residue & local Energy----------------
			randchain = rand() % chainNum;
			resNum = bundle->getNumResidues(randchain);
			randres = residues[rand() % resSize];

			//--backbone optimization----------------
			if (randres < (resNum-1) && randres > 0)
			{
				//--folding directions
				coords1 = bundle->getCoords(0, 8, "CA");
				coords2 = bundle->getCoords(1, 8, "CA");
				coords3 = bundle->getCoords(2, 8, "CA");
				dist1 = CMath::distance(coords1, coords2);
				dist2 = CMath::distance(coords1, coords3);
				dist3 = CMath::distance(coords2, coords3);
				sum = (dist1+dist2+dist3)/2;
				preArea = sqrt(sum * (sum-dist1) * (sum-dist2) * (sum-dist3));
						
				//--set phi or psi with NtoC or CtoN, local or global transformation
				direction = 0, randtype = rand() % 2, distance = rand() % 2;
				do 
				{
					deltaTheta = ((rand() % 3) -1);
				} while (deltaTheta == 0);

				currentangle = bundle->getAngle(randchain, randres, randtype);	
				bundle->setAngleLocal(randchain, randres, currentangle, deltaTheta, randtype, distance, direction);
				Energy = bundle->intraSoluteEnergy(true), test++;
			
				//--folding check				
				coords1 = bundle->getCoords(0, 8, "CA");
				coords2 = bundle->getCoords(1, 8, "CA");
				coords3 = bundle->getCoords(2, 8, "CA");
				dist1 = CMath::distance(coords1, coords2);
				dist2 = CMath::distance(coords1, coords3);
				dist3 = CMath::distance(coords2, coords3);
				sum = (dist1+dist2+dist3)/2;
				postArea = sqrt(sum * (sum-dist1) * (sum-dist2) * (sum-dist3));

				//--Energy test
				if (postArea > preArea && Energy < 0.0)	
				{					
					open++;			
					goodEnergy = Energy;
					energyFile << Energy << " ";
					areaFile << postArea << " ";
					test = 0;
				}

				//--Revert if not better										
				if (test != 0)
				{
					bundle->setAngleLocal(randchain, randres, currentangle, (deltaTheta*-1), randtype, distance, direction);
				}
			}
		} while (open < 200);
		energyFile << endl;
		areaFile << endl;
		stringstream convert; 
		string countstr;
		convert << a, countstr = convert.str();
		outPDB = countstr + ".pdb";
		pdbWriter(bundle, outPDB);
		delete thePDB;
	}
	return 0;
}

