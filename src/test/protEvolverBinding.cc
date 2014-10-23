//*******************************************************************************************************
//*******************************************************************************************************
//***********************************                          ******************************************
//***********************************     protEvolver 2.0      ******************************************
//***********************************                          ******************************************
//*******************************************************************************************************
//***************   -Folding Selective Protein Evolution in Implicit Solvent-   ***********************
//*******************************************************************************************************

/////// Just specify infile structure, active chains and residues indexes, and it will evolve a sequence favorable for folding

//--Included files and functions-------------------------------------------------------------------------
#include <iostream>
#include <string>
#include <time.h>
#include <sstream>
#include "ensemble.h"
#include "PDBInterface.h"

void randomizeSideChains(protein* _prot, UInt _chainIndex);
vector < UInt > getChainSequence(protein* _prot, UInt _chainIndex);

//--Program setup----------------------------------------------------------------------------------------
int main (int argc, char* argv[])
{
	//--Running parameters
	if (argc !=2)
	{
		cout << "protEvolver <inFile.pdb>" << endl;
		exit(1);
	}
	string infile = argv[1];
	enum aminoAcid {A,R,N,D,Dh,C,Cx,Q,E,Eh,Hd,He,Hn,Hp,I,L,K,M,F,P,O,S,T,W,Y,V,G,dA,dR,dN,dD,dDh,dC,dCx,dQ,dE,dEh,dHd,dHe,dHn,dHp,dI,dL,dK,dM,dF,dP,dO,dS,dT,dW,dY,dV};
	string aminoAcidString[] = {"A","R","N","D","Dh","C","Cx","Q","E","Eh","Hd", "He","Hn","Hp","I","L","K","M","F","P","O","S","T","W","Y", "V","G","dA","dR","dN","dD","dDh","dC","dCx","dQ","dE","dEh","dHd","dHe","dHn","dHp","dI","dL","dK","dM","dF","dP","dO","dS","dT","dW","dY","dV"};
	PDBInterface* thePDB = new PDBInterface(infile);
	ensemble* theEnsemble = thePDB->getEnsemblePointer();
	molecule* pMol = theEnsemble->getMoleculePointer(0);
	protein* bundle = static_cast<protein*>(pMol);
	bundle->silenceMessages();
	residue::setCutoffDistance(9.0);
	pmf::setScaleFactor(0.0);
	rotamer::setScaleFactor(0.0);
	microEnvironment::setScaleFactor(0.0);
	amberVDW::setScaleFactor(1.0);
	amberVDW::setRadiusScaleFactor(1.0);
	amberVDW::setLinearRepulsionDampeningOff();
	amberElec::setScaleFactor(1.0);
	solvation::setItsScaleFactor(0.0);
	srand (time(NULL));

	//--inputs for mutation
	UInt activeChains[] = {1};
	//UInt allowedLResidues[] = {A,R,N,D,Q,E,He,I,L,K,M,F,P,S,T,W,Y,V,G};
	//UInt activeResidues[] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19};
	UInt allowedLResidues[] = {A,N,D,Q,E,He,L,K,M,F,S,W,Y};
	UInt activeResidues[] = {1,4,5,6};
	//UInt allowedLResidues[] = {A,R,N,D,Q,E,He,I,L,K,M,F,P,S,T,W,Y,V,G};
	//UInt activeResidues[] = {0,1,2,3,5,6,7,9,10,12,13,14,15,16,17,18,19,21,22,24,25,26,28,29,31};
	//UInt allowedLResidues[] = {E,R,K,D,P,P};
	//UInt activeResidues[] = {0,1,3,4,6,7,9,10,12,13,15,16,18,19,21,22,24,25,27,28};
	UInt allowedDResidues[] = {dA,dR,dN,dD,dQ,dE,dHe,dI,dL,dK,dM,dF,dP,dS,dT,dW,dY,dV,G};

	double phi, bestEnergy, foldingEnergy, posE, totalposE, aveposE;
	UInt nobetter = 0, activeResiduesSize = sizeof(activeResidues)/sizeof(activeResidues[0]), activeChainsSize = sizeof(activeChains)/sizeof(activeChains[0]);
	UInt dResidues = sizeof(allowedDResidues)/sizeof(allowedDResidues[0]), lResidues = sizeof(allowedLResidues)/sizeof(allowedLResidues[0]);
	UInt name, mutant, numResidues, plateau = (lResidues*activeResiduesSize*activeChainsSize*1.5), fib = 0;
	vector < UInt > mutantPosition, chainSequence, sequencePosition;
	vector < vector < UInt > > proteinSequence, finalSequence;
	vector < double > Energy(2);
	UInt randres, randchain, rescount = 0;
	stringstream convert;
	string startstr, outFile;
	name = rand() % 1000000;
	convert << name, startstr = convert.str();
	string tempModel = startstr + "_temp.pdb";
	delete thePDB;
	cout << "name complexEnergy complexFoldingBindingEnergy complexBindingEnergy ligandEnergy ligandFoldingEnergy sequence" << endl;

	//--Run multiple independent evolution cycles-----------------------------------------------------
	for (UInt a = 1; a < 20; a++)
	{
		PDBInterface* thePDB = new PDBInterface(infile);
		ensemble* theEnsemble = thePDB->getEnsemblePointer();
		molecule* pMol = theEnsemble->getMoleculePointer(0);
		protein* bundle = static_cast<protein*>(pMol);

		//--load in initial pdb and mutate in random starting sequence on active chains and residues
		proteinSequence.clear(), chainSequence.clear(), Energy.clear(), fib = 0, nobetter = 0;
		for (UInt i = 0; i < activeChainsSize; i++)
		{
			for (UInt j = 0; j < activeResiduesSize; j++)
			{
				bundle->activateForRepacking(activeChains[i], activeResidues[j]);
				phi = bundle->getPhi(activeChains[i], activeResidues[j]);
				if (phi > 0 && phi < 180)
				{
					mutant = allowedDResidues[(rand() % dResidues)];
					bundle->mutateWBC(activeChains[i], activeResidues[j], mutant);
				}
				if (phi < 0 && phi > -180)
				{
					mutant = allowedLResidues[(rand() % lResidues)];
					bundle->mutateWBC(activeChains[i], activeResidues[j], mutant);
				}
			}
			randomizeSideChains(bundle, activeChains[i]);
			chainSequence = getChainSequence(bundle, activeChains[i]);
			proteinSequence.push_back(chainSequence);
		}
		bundle->protOptSolvent(300);

		//--Determine next mutation position
		rescount = 0, totalposE = 0, mutantPosition.clear();
		for (UInt i = 0; i < activeChainsSize; i++)
		{
			for (UInt j = 0; j < activeResiduesSize; j++)
			{
				rescount++;
				posE = bundle->getPositionSoluteEnergy(activeChains[i], activeResidues[j], true);
				totalposE = totalposE + posE;
			}
		}
		aveposE = totalposE/rescount;
		posE = -1E10;
		do
		{
			randres = activeResidues[rand() % activeResiduesSize];
			if (activeChainsSize > 1)
			{
				randchain = activeChains[rand() % activeChainsSize];
			}
			else
			{
				randchain = activeChains[0];
			}
			posE = bundle->getPositionSoluteEnergy(randchain, randres, true);
		}while (posE <= aveposE);
		mutantPosition.push_back(randchain);
		mutantPosition.push_back(randres);
		pdbWriter(bundle, tempModel);

		//--set Energy startpoint
		Energy.clear();
		Energy = bundle->chainFoldingBindingEnergy(false);
		foldingEnergy = Energy[1];
		bestEnergy = foldingEnergy;
		delete thePDB;
		
		//--Run single evolution loop till hitting plateau
		do
		{  
			PDBInterface* thePDB = new PDBInterface(infile);
			ensemble* theEnsemble = thePDB->getEnsemblePointer();
			molecule* pMol = theEnsemble->getMoleculePointer(0);
			protein* bundle = static_cast<protein*>(pMol);

			//--Mutate current sequence, new mutant and optimize system
			fib++, nobetter++;
			for (UInt i = 0; i < activeChainsSize; i++)
			{
				numResidues = bundle->getNumResidues(activeChains[i]);
				for (UInt j = 0; j < numResidues; j++)
				{
					bundle->activateForRepacking(activeChains[i],j);
					if (activeChains[i] == mutantPosition[0] && j == mutantPosition[1])
					{
						//--new mutant
						sequencePosition.push_back(i);
						sequencePosition.push_back(j);
						phi = bundle->getPhi(mutantPosition[0], mutantPosition[1]);
						if (phi > 0 && phi < 180)
						{
							mutant = allowedDResidues[(rand() % dResidues)];
							bundle->mutateWBC(mutantPosition[0],mutantPosition[1], mutant);
						}
						if (phi < 0 && phi > -180)
						{
							mutant = allowedLResidues[(rand() % lResidues)];
							bundle->mutateWBC(mutantPosition[0],mutantPosition[1], mutant);
						}
					}
					else
					{
						bundle->mutateWBC(activeChains[i],j, proteinSequence[i][j]);
					}
				}
				randomizeSideChains(bundle, activeChains[i]);
			}
			bundle->protOptSolvent(300);
			protein* tempBundle = new protein(*bundle);
			
			//--Determine next mutation position
			rescount = 0, totalposE = 0, mutantPosition.clear();
			for (UInt i = 0; i < activeChainsSize; i++)
			{
				for (UInt j = 0; j < activeResiduesSize; j++)
				{
					rescount++;
					posE = bundle->getPositionSoluteEnergy(activeChains[i], activeResidues[j], true);
					totalposE = totalposE + posE;
				}
			}
			aveposE = totalposE/rescount;
			posE = -1E10;
			do
			{
				randres = activeResidues[rand() % activeResiduesSize];
				if (activeChainsSize > 1)
				{
					randchain = activeChains[rand() % activeChainsSize];
				}
				else
				{
					randchain = activeChains[0];
				}
				posE = bundle->getPositionSoluteEnergy(randchain, randres, true);
			}while (posE <= aveposE);
			mutantPosition.push_back(randchain);
			mutantPosition.push_back(randres);

			//--Energy test and determination of next mutant position
			Energy.clear();
			Energy = bundle->chainFoldingBindingEnergy(false);
			if ((Energy[1] > Energy[0]) && (Energy[1] < (foldingEnergy+fib)))
			{
				if (Energy[1] < bestEnergy)
				{
					bestEnergy = Energy[1];
					pdbWriter(tempBundle, tempModel);
				}
				fib = 0, nobetter--, proteinSequence[sequencePosition[0]][sequencePosition[1]] = mutant, foldingEnergy = Energy[1];
			}
			sequencePosition.clear();
			delete thePDB;
			delete tempBundle;
		}while (nobetter < plateau);

		//--Print final energy and write a pdb file----------------------------------------------------
		PDBInterface* theModelPDB = new PDBInterface(tempModel);
		ensemble* theModelEnsemble = theModelPDB->getEnsemblePointer();
		molecule* modelMol = theModelEnsemble->getMoleculePointer(0);
		protein* model = static_cast<protein*>(modelMol);
		name = rand() % 1000000;
		stringstream convert;
		string countstr;
		convert << name, countstr = convert.str();
		outFile = countstr + ".evo.pdb";
		pdbWriter(model, outFile);
		finalSequence.clear(), chainSequence.clear(), Energy.clear();
		for (UInt i = 0; i < activeChainsSize; i++)
		{
			chainSequence = getChainSequence(model, activeChains[i]);
			finalSequence.push_back(chainSequence);
		}
		Energy = model->chainFoldingBindingEnergy(activeChains[0]);
		cout << name << " " << Energy[0] << " " << Energy[1] << " " << Energy[2] << " " << Energy[3] << " " << Energy[4] << " ";
		for (UInt i = 0; i < activeChainsSize; i++)
		{
			for (UInt j = 0; j < finalSequence[i].size(); j++)
			{
				cout << aminoAcidString[finalSequence[i][j]] << " ";
			}
		}
		cout << endl;
		delete theModelPDB;
	}
	return 0;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////
//////// functions //////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////

void randomizeSideChains(protein* _prot, UInt _chainIndex)
{	
	UInt allowedRotsSize, randrot, restype, numResidues;
	UIntVec allowedRots;
	numResidues = _prot->getNumResidues(_chainIndex);
	for (UInt j = 0; j < numResidues; j++)
	{
		restype = _prot->getTypeFromResNum(_chainIndex, j);
		allowedRots = _prot->getAllowedRotamers(_chainIndex, j, restype, 0);
		allowedRotsSize = allowedRots.size();
		if (allowedRotsSize > 2)
		{
			randrot = rand() % allowedRotsSize;
			_prot->setRotamerWBC(_chainIndex, j, 0, allowedRots[randrot]);
		}
	}
	return;
}

vector < UInt > getChainSequence(protein* _prot, UInt _chainIndex)
{	
	UInt restype, numResidues;
	vector < UInt > sequence;
	numResidues = _prot->getNumResidues(_chainIndex);
	for (UInt j = 0; j < numResidues; j++)
	{
		restype = _prot->getTypeFromResNum(_chainIndex, j);
		sequence.push_back(restype);
	}
	return sequence;
}
