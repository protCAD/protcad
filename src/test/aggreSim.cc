//*******************************************************************************************************
//*******************************************************************************************************
//**************************************                 ************************************************
//**************************************   aggreSim 1.0  ************************************************
//**************************************                 ************************************************
//*******************************************************************************************************
//**************************    -six-axis aggregation simulation-    ************************************
//*******************************************************************************************************

/////// Just specify a 2 chain infile with merged centroids and it will result in aggregates.

//--Included files and functions-------------------------------------------------------------------------
#include <iostream>
#include <string>
#include <time.h>
#include <sstream>
#include "ensemble.h"
#include "PDBInterface.h"

dblVec carbonCentroid(protein* _prot, UInt _chain);
void diffuseComplex(protein* _prot);

//--Program setup----------------------------------------------------------------------------------------
int main (int argc, char* argv[])
{
	//--Running parameters
	if (argc !=2)
	{
		cout << "aggreSim <inFile.pdb>" << endl;
		exit(1);
	}
	string infile = argv[1];
	enum aminoAcid {A,R,N,D,Dh,C,Cx,Q,E,Eh,Hd,He,Hn,Hp,I,L,K,M,F,P,O,S,T,W,Y,V,G,dA,dR,dN,dD,dDh,dC,dCx,dQ,dE,dEh,dHd,dHe,dHn,dHp,dI,dL,dK,dM,dF,dP,dO,dS,dT,dW,dY,dV};
	PDBInterface* thePDB = new PDBInterface(infile);
	ensemble* theEnsemble = thePDB->getEnsemblePointer();
	molecule* pMol = theEnsemble->getMoleculePointer(0);
	protein* bundle = static_cast<protein*>(pMol);
	bundle->silenceMessages();
	residue::setCutoffDistance(10.0);
	pmf::setScaleFactor(0.0);
	rotamer::setScaleFactor(0.0);
	microEnvironment::setScaleFactor(0.0);
	amberVDW::setScaleFactor(1.0);
	amberVDW::setRadiusScaleFactor(1.0);
	amberVDW::setLinearRepulsionDampeningOff();
	amberElec::setScaleFactor(1.0);
	solvation::setItsScaleFactor(0.0);
	srand (time(NULL));

	//--initialize variables for loop
	double Energy, pastEnergy, rotx, roty, rotz, transx, transy, transz;
	double totaldistOld, totaldistNew, finaldist, finaldistave;
	UInt name, totalsize = 50, nobetter, test;
	dblVec chain0Centroid, chain1Centroid;
	string outFile;
	delete thePDB;
	//cout << "sim# " << "distance " << "totalEnergy " << "totalBindingEnergy " << "187_bindingEnergy" << endl;
	
	//--Loop for multiple simulations
	for (int a = 1; a < 500; a++)
	{
		PDBInterface* thePDB = new PDBInterface(infile);
		ensemble* theEnsemble = thePDB->getEnsemblePointer();
		molecule* pMol = theEnsemble->getMoleculePointer(0);
		protein* bundle = static_cast<protein*>(pMol);
		diffuseComplex(bundle);
		pastEnergy = bundle->intraSoluteEnergy(true);
		chain0Centroid = carbonCentroid(bundle, 0);
		chain1Centroid = carbonCentroid(bundle, 1);
		nobetter = 0, test = 1;
		//--Run optimizaiton loop till grand minimum---------------------------------------------------
		do
		{  
			//--Move randomly constrained by energy and distance condition
			do
			{
				nobetter++;
				chain1Centroid = carbonCentroid(bundle, 1);
				totaldistOld = CMath::distance(chain1Centroid, chain0Centroid);
				rotx = (rand() % 3)-1, roty = (rand() % 3)-1, rotz = (rand() % 3)-1;
				transx = (rand() % 3)-1, transy = (rand() % 3)-1, transz = (rand() % 3)-1;
				do
				{
					//--Make random motion until Energy rises-----------
					bundle->rotateChain(1, X_axis, rotx), bundle->rotateChain(1, Y_axis, roty), bundle->rotateChain(1, Z_axis, rotz);
					bundle->translateChain(1, transx, 0, 0), bundle->translateChain(1, 0, transy, 0), bundle->translateChain(1, 0, 0, transz);
					test++;

					//--Get new distance and Energy
					chain1Centroid = carbonCentroid(bundle, 1);
					totaldistNew = CMath::distance(chain1Centroid, chain0Centroid);		
					Energy = bundle->intraSoluteEnergy(true);		
					if (totaldistNew < totaldistOld && Energy < (pastEnergy + (fabs(pastEnergy/9))))
					{	
						//cout << Energy << endl;
						nobetter = 0, test = 0, pastEnergy = Energy, totaldistOld = totaldistNew;
					}

					//--Revert if no improvement
					if (test != 0)
					{
						bundle->translateChain(1, 0, 0, (transz * -1)), bundle->translateChain(1, 0, (transy * -1), 0), bundle->translateChain(1, (transx * -1), 0, 0);
						bundle->rotateChain(1, Z_axis, (rotz * -1)), bundle->rotateChain(1, Y_axis, (roty * -1)), bundle->rotateChain(1, X_axis, (rotx * -1));
					}	
				}while (test == 0);
			}while (nobetter < (totalsize * .75));
			bundle->protOptSolvent(50);
		}while (nobetter < (totalsize * 1.5));
		bundle->protOptSolvent(150);
		
		//--Print final energy and write a pdb file----------------------------------------------------
		finaldist = CMath::distance(chain1Centroid, chain0Centroid);
		finaldistave = finaldist;
		name = rand() % 1000000;
		cout << name << " " << finaldistave << " " << bundle->intraSoluteEnergy(true) << " " << bundle->bindingSoluteEnergy(0,1) << endl;
		stringstream convert; 
		string countstr;
		convert << name, countstr = convert.str();
		outFile = countstr + ".pdb";
		pdbWriter(bundle, outFile);
		delete thePDB;
	}
	cout << "Complete" << endl << endl;
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


///// diffuse complex randomly until there is a null interaction energy /////////////////////////////////
void diffuseComplex(protein* _prot)
{
	//--Initialize variables for loop
	double Energy, rotx, roty, rotz, transx, transy, transz;
	rotx = (rand() % 720)-360, roty = (rand() % 720)-360, rotz = (rand() % 720)-360;
	do 
	{
		transx = (rand() % 15)-7;
		transy = (rand() % 15)-7;
		transz = (rand() % 15)-7;
	} while (transx == 0 && transy == 0 && transz == 0);
	
	//--Make random motion and test for null interaction
	do 
	{	
		_prot->rotateChain(1, X_axis, rotx);
		_prot->rotateChain(1, Y_axis, roty);
		_prot->rotateChain(1, Z_axis, rotz);
		_prot->translateChain(1, transx, 0, 0);
		_prot->translateChain(1, 0, transy, 0);
		_prot->translateChain(1, 0, 0, transz);
		Energy = _prot->bindingSoluteEnergy(0,1);
	} while (Energy != 0.0);
	return;
}  

		
