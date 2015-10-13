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
#include <unistd.h>
#include "ensemble.h"
#include "PDBInterface.h"

void randomizeSideChains(protein* _prot, UInt _chainIndex);
vector <UInt> getChainSequence(protein* _prot, UInt _chainIndex);
vector <UInt> getMutationPosition(protein* _prot, UInt* _activeChains, UInt* _activeResidues, UInt _nobetter);

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
    enum aminoAcid {A,R,N,D,Dh,C,Cx,Q,E,Eh,Hd,He,Hn,Hp,I,L,K,M,F,P,O,S,T,W,Y,V,G,dA,dR,dN,dD,dDh,dC,dCx,dQ,dE,dEh,dHd,dHe,dHn,dHp,dI,dL,dK,dM,dF,dP,dO,dS,dT,dAT,dW,dY,dV,Hce,Pch};
    string aminoAcidString[] = {"A","R","N","D","Dh","C","Cx","Q","E","Eh","Hd","He","Hn","Hp","I","L","K","M","F","P","O","S","T","W","Y", "V","G","dA","dR","dN","dD","dDh","dC","dCx","dQ","dE","dEh","dHd","dHe","dHn","dHp","dI","dL","dK","dM","dF","dP","dO","dS","dT","dAT","dW","dY","dV","Hce","Pch"};
    PDBInterface* thePDB = new PDBInterface(infile);
    ensemble* theEnsemble = thePDB->getEnsemblePointer();
    molecule* pMol = theEnsemble->getMoleculePointer(0);
    protein* bundle = static_cast<protein*>(pMol);
    bundle->silenceMessages();
    residue::setCutoffDistance(9.0);
    rotamer::setScaleFactor(0.0);
    amberVDW::setScaleFactor(1.0);
    amberVDW::setRadiusScaleFactor(0.95);
    amberVDW::setLinearRepulsionDampeningOff();
    amberElec::setScaleFactor(1.0);
    srand (getpid());

	//--inputs for mutation
	UInt activeChains[] = {0};
    UInt allowedLResidues[] = {A,R,K,L,E,D,I,V,P,G};
    UInt activeResidues[] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30};
    UInt allowedDResidues[] = {dA,dR,dK,dL,dE,dD,dI,dV,dP,G};

    double phi, bestEnergy, pastEnergy, finalEnergy, Energy;
	UInt nobetter = 0, activeResiduesSize = sizeof(activeResidues)/sizeof(activeResidues[0]), activeChainsSize = sizeof(activeChains)/sizeof(activeChains[0]);
	UInt dResidues = sizeof(allowedDResidues)/sizeof(allowedDResidues[0]), lResidues = sizeof(allowedLResidues)/sizeof(allowedLResidues[0]);
    UInt name, mutant = 0, numResidues, plateau = (lResidues*activeResiduesSize*activeChainsSize);
	vector < UInt > mutantPosition, chainSequence, sequencePosition;
	vector < vector < UInt > > proteinSequence, finalSequence;
	stringstream convert;
	string startstr, outFile;
	name = rand() % 1000000;
	convert << name, startstr = convert.str();
	string tempModel = startstr + "_temp.pdb";
	delete thePDB;

	//--Run multiple independent evolution cycles-----------------------------------------------------
	for (UInt a = 1; a < 40; a++)
	{
		PDBInterface* thePDB = new PDBInterface(infile);
		ensemble* theEnsemble = thePDB->getEnsemblePointer();
		molecule* pMol = theEnsemble->getMoleculePointer(0);
		protein* bundle = static_cast<protein*>(pMol);

		//--load in initial pdb and mutate in random starting sequence on active chains and residues
        proteinSequence.clear(), chainSequence.clear(), nobetter = 0;
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
            //randomizeSideChains(bundle, activeChains[i]);
			chainSequence = getChainSequence(bundle, activeChains[i]);
			proteinSequence.push_back(chainSequence);
		}
        bundle->protOpt(false);

		//--Determine next mutation position
        mutantPosition.clear();
        mutantPosition = getMutationPosition(bundle, activeChains, activeResidues, nobetter);
		pdbWriter(bundle, tempModel);

		//--set Energy startpoint
        Energy = bundle->protEnergy();
        pastEnergy = Energy;
        bestEnergy = Energy;
		delete thePDB;
		
		//--Run single evolution loop till hitting plateau
        //#pragma omp parallel while
		do
		{  
			PDBInterface* thePDB = new PDBInterface(infile);
			ensemble* theEnsemble = thePDB->getEnsemblePointer();
			molecule* pMol = theEnsemble->getMoleculePointer(0);
			protein* bundle = static_cast<protein*>(pMol);

			//--Mutate current sequence, new mutant and optimize system
            nobetter++;
			for (UInt i = 0; i < activeChainsSize; i++)
			{
				numResidues = bundle->getNumResidues(activeChains[i]);
				for (UInt j = 0; j < numResidues; j++)
				{
                    bundle->setMoved(activeChains[i], j,1);
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
                //randomizeSideChains(bundle, activeChains[i]);
			}
            bundle->protOpt(false);
			protein* tempBundle = new protein(*bundle);
			
			//--Determine next mutation position
            mutantPosition.clear();
            mutantPosition = getMutationPosition(bundle, activeChains, activeResidues, nobetter);

			//--Energy test and determination of next mutant position
            Energy = bundle->protEnergy();
            if (Energy < (pastEnergy+nobetter))
			{
                //cout << Energy << " " << nobetter << endl;
                if (Energy < bestEnergy)
				{
                    bestEnergy = Energy;
					pdbWriter(tempBundle, tempModel);
				}
                proteinSequence[sequencePosition[0]][sequencePosition[1]] = mutant, pastEnergy = Energy;
                if (nobetter > 0) { nobetter--;
                }
                else{ nobetter = 0;
                }
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
        finalSequence.clear(), chainSequence.clear();
		for (UInt i = 0; i < activeChainsSize; i++)
		{
			chainSequence = getChainSequence(model, activeChains[i]);
			finalSequence.push_back(chainSequence);
		}
        Energy = model->protEnergy();
        finalEnergy = Energy;
        cout << name << " " << finalEnergy << " ";
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
vector <UInt> getMutationPosition(protein* _prot, UInt *_activeChains, UInt *_activeResidues, UInt _nobetter)
{
    //--get median residue energy
    UInt randres, randchain, _activeResiduesSize = sizeof(_activeResidues)/sizeof(_activeResidues[0]), _activeChainsSize = sizeof(_activeChains)/sizeof(_activeChains[0]);
    double posE, medE;
    vector <UInt> _mutantPosition;
    medE = _prot->getMedianResEnergy();

    //--find random position with worse than median energy
    posE = -1E10;
    do
    {
        randchain = _activeChains[rand() % _activeChainsSize];
        randres = _activeResidues[rand() % _activeResiduesSize];
        posE = _prot->resEnergy(randchain,randres);
    }while (posE < (medE/(_nobetter+1)));
    _mutantPosition.push_back(randchain);
    _mutantPosition.push_back(randres);
    return _mutantPosition;
}

