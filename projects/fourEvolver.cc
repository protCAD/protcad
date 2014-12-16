//*******************************************************************************************************
//*******************************************************************************************************
//***********************************                   *************************************************
//***********************************  protEvolver 1.5  *************************************************
//***********************************                   *************************************************
//*******************************************************************************************************
//***************   -Stability Selective Protein Evolution in Implicit Solvent-   ***********************
//*******************************************************************************************************

/////// Just specify infile structure and it will evolve in hetero-oligameric stability

//--Included files and functions-------------------------------------------------------------------------
#include <iostream>
#include <string>
#include <time.h>
#include <sstream>
#include "ensemble.h"
#include "PDBInterface.h"

void randomizeSideChains(protein* _prot, UInt _chainIndex);
vector < UInt > getSequence(protein* _prot, UInt _chainIndex);
vector < UInt > getMutationPosition(protein* _prot, UInt _chainIndex);

//--Program setup----------------------------------------------------------------------------------------
int main (int argc, char* argv[])
{

cout << "Program Started" << endl;

	//--Running parameters
	if (argc !=2)
	{
		cout << "fourEvolver <inFile.pdb>" << endl;
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


	//--mutation input
	UInt allowedLResidues[] = {A,F,I,L,M,T,V,W,Y};
	UInt allowedDResidues[] = {dA};
    //UInt activeChain[] = {0,1,2,3};
    //UInt randomChain;

cout << "Loaded Mutation Candidates" << endl;

    double pastEnergy, bestEnergy, Energy, origEnergy = bundle->intraSoluteEnergy(true);
	UInt resCount, name, mutant, nobetter, numResidues, sequencePosition, totalRes = 0, plateau, fib = 0;
	UInt dResidues = sizeof(allowedDResidues)/sizeof(allowedDResidues[0]), lResidues = 9; //sizeof(allowedLResidues)/sizeof(allowedLResidues[0]);
	vector < UInt > mutantPosition, position, finalSequence, chainSequence;
	vector < vector < UInt > > sequence, originalSequence;

	for (UInt i = 0; i < bundle->getNumChains(); i++)
	{	chainSequence = getSequence(bundle, i);
		sequence.push_back(chainSequence);
		chainSequence.clear();
	}

cout << "Created Mutation Vector" << endl;

	stringstream convert;
	string startstr, outFile;
	name = rand() % 1000000;
	convert << name, startstr = convert.str();
	string goodModel = startstr + "_temp.pdb";

cout << "Temporary PDB Named" << endl;

	pdbWriter(bundle, goodModel);
cout << "Temporary PDB Generated" << endl;

	
	//--set energy plateau for evolution cycle end
	totalRes = 30;
	if (lResidues > 1)
	{
        plateau = (50);
        //plateau = 1;
	}
	else
	{
		plateau = (totalRes * dResidues);
	}
	 
	//--Run multiple independent evolution cycles-----------------------------------------------------
	for (UInt a = 1; a < 100; a++)
	{
		//--load in initial pdb
		sequence.clear();
		sequence = originalSequence;

        mutantPosition = getMutationPosition(bundle, 0);

		pastEnergy = origEnergy;
		bestEnergy = origEnergy;
		nobetter = 0, fib = 0;

		
		//--Run single evolution loop till hitting plateau
		do
		{  
			PDBInterface* thePDB = new PDBInterface(infile);
			ensemble* theEnsemble = thePDB->getEnsemblePointer();
			molecule* pMol = theEnsemble->getMoleculePointer(0);
			protein* bundle = static_cast<protein*>(pMol);
			
cout << "Started single evolution" << endl;

			//--Mutate current sequence, new mutant and optimize system
            mutantPosition = getMutationPosition(bundle, 0);
			resCount = 0, nobetter++, fib++;
			numResidues = bundle->getNumResidues(mutantPosition[0]);
            bundle->activateAllForRepacking(mutantPosition[0]);

			for (UInt j = 0; j < numResidues; j++)
			{
cout << "Mutant Position Chosen" << endl;
				if (j == mutantPosition[1])
				{
					//--new mutant
					sequencePosition = mutantPosition[1];
cout << "Sequence Position for mutation: " << sequencePosition <<endl;
                    mutant = allowedLResidues[(rand() % lResidues)];
cout << "Residue chosen" << endl;
                    bundle->mutateWBC(mutantPosition[0],j, mutant);

				}
                /*else
				{	
cout << "Reverting back to original sequence" << endl;
                    //bundle->mutateWBC(mutantPosition[0], j, sequence[mutantPosition[0]][j]);
                }*/
				resCount++;
			}
            //randomizeSideChains(bundle, mutantPosition[0]);
            bundle->protOptSolvent(25);

			//--Energy test and determination of next mutant position
			Energy = bundle->intraSoluteEnergy(true);
			cout << Energy << endl;
			if (Energy < (pastEnergy + fib) && (Energy < (pastEnergy-.5) || Energy > (pastEnergy+.5)))
			{	
				if (Energy < bestEnergy)
				{
					bestEnergy = Energy;
					pdbWriter(bundle, goodModel);
cout << "Temporary PDB Generated" << endl;
				}
				fib = 0, nobetter--, sequence[mutantPosition[0]][sequencePosition] = mutant, pastEnergy = Energy;
				mutantPosition.clear();
                mutantPosition = getMutationPosition(bundle, 0);
			}
			delete thePDB;
		}while (nobetter < plateau);

		//--Print final energy and write a pdb file----------------------------------------------------
		PDBInterface* theModelPDB = new PDBInterface(goodModel);
		ensemble* theModelEnsemble = theModelPDB->getEnsemblePointer();
		molecule* modelMol = theModelEnsemble->getMoleculePointer(0);
		protein* model = static_cast<protein*>(modelMol);
		finalSequence.clear();
        finalSequence = getSequence(model, 0);
		name = rand() % 1000000;
		cout << name << " " << model->intraSoluteEnergy(true) << " " << model->intraSoluteEnergy(true) << " " << endl;
		for (UInt i = 0; i < finalSequence.size(); i++)
		{
			cout << aminoAcidString[finalSequence[i]] << " ";
		}
		cout << endl;
		stringstream convert; 
		string countstr;
		convert << name, countstr = convert.str();
		outFile = countstr + ".evo.pdb";
		pdbWriter(model, outFile);
		delete theModelPDB;
	}
	cout << "Complete" << endl << endl;
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

vector < UInt > getSequence(protein* _prot, UInt _chainIndex)
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

vector < UInt > getMutationPosition(protein* _prot, UInt _chainIndex)
{	
	//--get average position energy

	UInt randomChain;
	UInt residues[] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,18,19,20,21,22,23,24,25,26,27,28,29,30};

cout << "Allowed Residues chosen" << endl;

		for (UInt j = 0; j < sizeof(residues)/sizeof(residues[j]); j++)
		{
			cout << residues[j] << " ";
		}
		cout << endl;


	randomChain = rand() % 4;
	UInt residuesSize = sizeof(residues)/sizeof(residues[0]);
	UInt randres, rescount = 0;
	double posE, totalposE = 0, aveposE;
	vector < UInt > position;
	
	for (UInt j = 0; j < residuesSize; j++)
	{
		rescount++;
		posE = _prot->bindingPositionSoluteEnergy(randomChain, residues[j], 0);
		totalposE = totalposE + posE;
	}
	aveposE = totalposE/rescount;

cout << "Created Position Solute Energies Vector" << endl;

	//--find position with worse than average energy
	posE = -1E10;

cout << "Set Baseline Starting Energy" << endl;

	do
	{
cout << "do start" << endl;
		randres = residues[rand() % residuesSize];
cout << "randres set" << endl;
		posE = _prot->bindingPositionSoluteEnergy(randomChain, randres, 0);
cout << "posE set" << endl;

	}while (posE <= aveposE);
	position.push_back(randomChain);
cout << "Chose Random Chain" << endl;
	position.push_back(randres);
cout << "Chose Random Residue" << endl;
	return position;
}
	




