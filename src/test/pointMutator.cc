//*******************************************************************************************************
//*******************************************************************************************************
//***********************************                    ************************************************
//***********************************  pointMutator 1.5  ************************************************
//***********************************                    ************************************************
//*******************************************************************************************************
//**************   -point mutations, then backbone and sidechain optimization-   ************************
//*******************************************************************************************************

/////// Just specify a infile and preferred outfile name.

//--Program setup----------------------------------------------------------------------------------------
#include <iostream>
#include <string>
#include <time.h>
#include <sstream>
#include "ensemble.h"
#include "PDBInterface.h"

void protOpt (protein* _prot);

int main (int argc, char* argv[])
{
	//--Running parameters
	if (argc !=2)
	{
		cout << "pointMutator <inFile.pdb>" << endl;
		exit(1);
	}
	enum aminoAcid {A, R, N, D, Dh, C, Q, E, Eh, G, H, O, I, L, K, M, F, P, S, T, W, Y, V, dA, dR, dN, dD, dDh, dC, dQ, dE, dEh, dH, dO, dI, dL, dK, dM, dF, dP, dS, dT, dW, dY, dV};
	string residuenames[] = {"A","R","N","D","Dh","C","Q","E", "Eh","G","H","O","I","L","K","M","F","P","S","T","W","Y","V","dA","dR","dN","dD","dDh","dC","dQ","dE","dEh","dH","dO","dI","dL","dK","dM","dF","dP","dS","dT","dW","dY","dV"};
	string infile = argv[1];
	PDBInterface* thePDB = new PDBInterface(infile);
	ensemble* theEnsemble = thePDB->getEnsemblePointer();
	molecule* pMol = theEnsemble->getMoleculePointer(0);
	protein* bundle = static_cast<protein*>(pMol);
	bundle->silenceMessages();
	residue::setCutoffDistance(12.0);
	pmf::setScaleFactor(0.0);
	rotamer::setScaleFactor(0.0);
	microEnvironment::setScaleFactor(0.0);
	amberVDW::setScaleFactor(1.0);
	amberVDW::setRadiusScaleFactor(1.0);
	amberVDW::setLinearRepulsionDampeningOff();
	amberElec::setScaleFactor(1.0);
	solvation::setItsScaleFactor(0.0);
	srand (time(NULL));
	
//--Select and mutate sites------------------------------------------------------------------------------

	//--Selections
	int chains[] = {1,3,5,7,9,11};
	int chainsSize = sizeof(chains)/sizeof(chains[0]);
	int residues[] = {9};
	int resnums[] =  {10};
	int residuesSize = sizeof(residues)/sizeof(residues[0]);
	int resID[] = {E};
	int resIDsize = sizeof(resID)/sizeof(resID[0]);
	int replicates = 10;
	int count = 0;
	UInt mutant;
	cout << endl << "pdb " << "residue " << "site " << "mutant " << "Energy " << "interEnergy" << endl;

	//--Mutations
	for (int i = 0; i < residuesSize; i++)
	{
		for (int j = 0; j < resIDsize; j++)
		{
			for (int k = 0; k < replicates; k++)
			{
				PDBInterface* thePDB = new PDBInterface(infile);
				ensemble* theEnsemble = thePDB->getEnsemblePointer();
				molecule* pMol = theEnsemble->getMoleculePointer(0);
				protein* bundle = static_cast<protein*>(pMol);
				UInt restype = bundle->getTypeFromResNum(chains[0], (UInt)residues[i]);
				mutant = resID[j];
				for (int l = 0; l < chainsSize; l++)
				{
					bundle->activateForRepacking(chains[l], residues[i]);
					bundle->mutate(chains[l], residues[i], mutant);
				}
				protOpt(bundle);
				count++;
				double intra = bundle->intraEnergy(), sumIntra = 0, chainIntra, inter;
				for (int m = 0; m < chainsSize; m++)
				{
					chainIntra = bundle->intraEnergy(m);
					sumIntra = (chainIntra + sumIntra);
				}
				inter = intra - sumIntra;
				stringstream convert; 
				string countstr;
				convert << count, countstr = convert.str();
				string outFile = countstr + ".pdb";
				pdbWriter(bundle, outFile);
				cout << count << " " << residuenames[restype] << " " << resnums[i] << " " << residuenames[mutant] << " " << intra << " " << inter << endl;
				delete thePDB;				
			}
		}
	}
	return 0;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////
//////// functions //////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////

void protOpt (protein* _prot)
{
		//--Initialize variables for loop
	double deltaTheta = 0, currentangle, totalpreposE = 0, avepreposE = -1E10;
	double Energy, pastEnergy, preposE, currentposE;
	int direction, thisone, distance;
	UInt randchain, randres, randrestype, allowedRotsize, randrot, number = 0, nobetter = 0, check = 0;
	UInt resNum, randtype, chainNum = _prot->getNumChains(), totalsize = 100, rotbetter = 0, test = 0;
	vector < vector <double> > currentRot;
	UIntVec allowedRots;
	string outFile;
	
	//--Get starting energies
	pastEnergy = _prot->intraEnergy();

//--Run optimizaiton loop till grand minimum-------------------------------------------------------------
	do
	{  
		//--Generate random residue & local Energy----------------
		randchain = rand() % chainNum;
		resNum = _prot->getNumResidues(randchain);
		randres = rand() % resNum;
		randrestype = _prot->getTypeFromResNum(randchain, randres);
		nobetter++;

		//--backbone optimization----------------
		if (rotbetter >= totalsize && preposE > avepreposE)
		{
			//--set phi or psi with NtoC or CtoN, local or global transformation
			test = 0;
			for (UInt i = 0; i < 2; i++)
			{	
				direction = rand() % 2, randtype = rand() % 2, distance = rand() % 2;
				do 
				{
					deltaTheta = ((rand() % 7) -3);
				} while (deltaTheta == 0);
				do
				{
					preposE = _prot->getPositionEnergy(randchain, randres);
					currentangle = _prot->getAngle(randchain, randres, randtype);	
					_prot->setAngleLocal(randchain, randres, currentangle, deltaTheta, randtype, distance, direction);
					currentposE = _prot->getPositionEnergy(randchain, randres), thisone = 0, test++;

					//--Energy test
					if (currentposE < (preposE - .05))	
					{				
						//Energy = _prot->intraEnergy();					
						if (Energy < pastEnergy)
						{	
							//cout << Energy << endl;
							nobetter = 0, test--, thisone = 1, pastEnergy = Energy;
						}
					}
				} while (thisone == 1);

				//--Revert if not better										
				if (thisone == 0)
				{
					_prot->setAngleLocal(randchain, randres, currentangle, (deltaTheta*-1), randtype, distance, direction);
				}
				if (test < 0)
				{
					break;
				}
			}
		}

		//--Rotamer optimization------------
		preposE = _prot->getPositionEnergy(randchain, randres);
		if (preposE > avepreposE)
		{
			//--Get current rotamer and allowed	
			currentRot = _prot->getSidechainDihedrals(randchain, randres);
			allowedRots = _prot->getAllowedRotamers(randchain, randres, randrestype, 0);
			allowedRotsize = (allowedRots.size()/3), rotbetter++;			
			
			//--Try 1/3 of allowed rotamers keep first improvement or revert to previous			
			for (UInt j = 0; j < allowedRotsize; j ++)
			{
				randrot = rand() % allowedRots.size();
				_prot->setRotamerWBC(randchain, randres, 0, allowedRots[randrot]);
				currentposE = _prot->getPositionEnergy(randchain, randres), check++;
				if (currentposE < (preposE - .05))
				{	
					Energy = _prot->intraEnergy();
					if (Energy < pastEnergy)	
					{					
						//cout << Energy << endl;
						rotbetter --, nobetter = 0, check = 0, pastEnergy = Energy;
						break;
					}								
				}
				_prot->setSidechainDihedralAngles(randchain, randres, currentRot);		
			}
		}

		//--check status of optimization-----------
		if (number == totalsize)
		{
			number = 0, totalpreposE = 0;
		}
		number++, totalpreposE = (totalpreposE + preposE), avepreposE = (totalpreposE/number);
	} while (nobetter < totalsize * 1.2);	      
	return;
}

