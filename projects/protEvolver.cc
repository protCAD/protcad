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
#include <dirent.h>
#include <sstream>
#include <unistd.h>
#include "ensemble.h"
#include "PDBInterface.h"

void randomizeSideChains(protein* _prot, UInt _chainIndex);
vector <UInt> getChainSequence(protein* _prot, UInt _chainIndex);
vector <UInt> getMutationPosition(protein* _prot, UInt* _activeChains, UInt _activeChainsSize, UInt* _activeResidues, UInt _activeResiduesSize);
UInt getProbabilisticMutation(vector <UInt> _mutantPosition, UInt *_aminoacids, UInt aaSize);
UInt getNumberofSequences();

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
    amberVDW::setRadiusScaleFactor(1.0);
    amberVDW::setLinearRepulsionDampeningOff();
    amberElec::setScaleFactor(1.0);
    srand (getpid());

	//--inputs for mutation
	UInt activeChains[] = {0};
    UInt allowedLResidues[] = {A,R,N,D,Q,E,I,L,K,M,F,P,S,T,W,Y,V,G};
    UInt activeResidues[] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39};
    UInt randomResidues[] = {0,1,2,3,6,9,10,13,16,17,18,21,22,23,25,26,29,32,33,36,39};
    UInt allowedDResidues[] = {dA,dR,dN,dD,dQ,dE,dI,dL,dK,dM,dF,dP,dS,dT,dW,dY,dV,G};

    double phi, bestEnergy, pastEnergy, finalEnergy, Energy;
    UInt nobetter = 0, activeChainsSize = sizeof(activeChains)/sizeof(activeChains[0]), randomResiduesSize = sizeof(randomResidues)/sizeof(randomResidues[0]), activeResiduesSize = sizeof(activeResidues)/sizeof(activeResidues[0]);
    UInt lResidues = sizeof(allowedLResidues)/sizeof(allowedLResidues[0]), dResidues = sizeof(allowedDResidues)/sizeof(allowedDResidues[0]);
    UInt name, mutant = 0, numResidues, plateau = (lResidues*activeChainsSize);
    vector < UInt > mutantPosition, chainSequence, sequencePosition, randomPosition;
	vector < vector < UInt > > proteinSequence, finalSequence;
	stringstream convert;
	string startstr, outFile;
	name = rand() % 1000000;
	convert << name, startstr = convert.str();
	string tempModel = startstr + "_temp.pdb";
	delete thePDB;

	//--Run multiple independent evolution cycles-----------------------------------------------------
    for (UInt a = 1; a < 500; a++)
	{
		PDBInterface* thePDB = new PDBInterface(infile);
		ensemble* theEnsemble = thePDB->getEnsemblePointer();
		molecule* pMol = theEnsemble->getMoleculePointer(0);
		protein* bundle = static_cast<protein*>(pMol);

        //--load in initial pdb and mutate in random starting sequence on active chains and random residues
        proteinSequence.clear(), chainSequence.clear(), nobetter = 0;
		for (UInt i = 0; i < activeChainsSize; i++)
		{
            for (UInt j = 0; j < randomResiduesSize; j++)
			{
                bundle->activateForRepacking(activeChains[i], randomResidues[j]);
                randomPosition.push_back(activeChains[i]);
                randomPosition.push_back(randomResidues[j]);
                phi = bundle->getPhi(activeChains[i], randomResidues[j]);
				if (phi > 0 && phi < 180)
				{
                    mutant = getProbabilisticMutation(randomPosition, allowedDResidues, dResidues);
                    bundle->mutateWBC(activeChains[i], randomResidues[j], mutant);
				}
				if (phi < 0 && phi > -180)
				{
                    mutant = getProbabilisticMutation(randomPosition, allowedLResidues, lResidues);
                    bundle->mutateWBC(activeChains[i], randomResidues[j], mutant);
                }
                randomPosition.clear();
			}
			chainSequence = getChainSequence(bundle, activeChains[i]);
			proteinSequence.push_back(chainSequence);
		}
        bundle->protOpt(false);

		//--Determine next mutation position
        mutantPosition.clear();
        mutantPosition = getMutationPosition(bundle, activeChains, activeChainsSize, activeResidues, activeResiduesSize);
		pdbWriter(bundle, tempModel);

		//--set Energy startpoint
        Energy = bundle->deltaH();
        pastEnergy = Energy;
        bestEnergy = Energy;
		delete thePDB;
		
		//--Run single evolution loop till hitting plateau
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
                            mutant = getProbabilisticMutation(mutantPosition, allowedDResidues, dResidues);
							bundle->mutateWBC(mutantPosition[0],mutantPosition[1], mutant);
						}
						if (phi < 0 && phi > -180)
						{
                            mutant = getProbabilisticMutation(mutantPosition, allowedLResidues, lResidues);
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
            mutantPosition = getMutationPosition(bundle, activeChains, activeChainsSize, activeResidues, activeResiduesSize);

            //--Energy test
            Energy = bundle->deltaH();
            if (Energy < pastEnergy)
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
            UInt seqCount = getNumberofSequences();
            if (seqCount < 200)
            {
                plateau = 0;
            }
            else
            {
                plateau = (lResidues*activeChainsSize);
            }
		}while (nobetter < plateau);

		//--Print final energy and write a pdb file----------------------------------------------------
		PDBInterface* theModelPDB = new PDBInterface(tempModel);
		ensemble* theModelEnsemble = theModelPDB->getEnsemblePointer();
		molecule* modelMol = theModelEnsemble->getMoleculePointer(0);
		protein* model = static_cast<protein*>(modelMol);
        Energy = model->deltaH();
        if (Energy < 0)
        {
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
        }
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

vector <UInt> getMutationPosition(protein* _prot, UInt *_activeChains, UInt _activeChainsSize, UInt *_activeResidues, UInt _activeResiduesSize)
{
    //--get median residue energy
    UInt randres, randchain;
    double posE, medE;
    vector <UInt> _mutantPosition;
    medE = _prot->getMedianDeltaH();

    //--find random position with worse than median energy
    do
    {
        randchain = _activeChains[rand() % _activeChainsSize];
        randres = _activeResidues[rand() % _activeResiduesSize];
        posE = _prot->deltaH(randchain,randres);
    }while (posE < medE);
    _mutantPosition.push_back(randchain);
    _mutantPosition.push_back(randres);
    return _mutantPosition;
}

UInt getProbabilisticMutation(vector <UInt> _mutantPosition, UInt *_aminoacids, UInt aaSize)
{
    UInt mutant, chance;
    double acceptance;
    string inFrame;
    vector <UInt> resFreqs(55,1);
    DIR *pdir;
    struct dirent *pent;
    pdir=opendir(".");
    UInt count = 0;

    //--get sequence evolution results for position
    while ((pent=readdir(pdir)))
    {
        inFrame = pent->d_name;
        if (inFrame.find(".evo.pdb") != std::string::npos)
        {
            count++;
            PDBInterface* theModelPDB = new PDBInterface(inFrame);
            ensemble* theModelEnsemble = theModelPDB->getEnsemblePointer();
            molecule* modelMol = theModelEnsemble->getMoleculePointer(0);
            protein* model = static_cast<protein*>(modelMol);
            model->silenceMessages();

            UInt restype = model->getTypeFromResNum(_mutantPosition[0], _mutantPosition[1]);
            resFreqs[restype] = resFreqs[restype] + 1;
            delete theModelPDB;
        }
    }
    closedir(pdir);

    //--find max frequency for position
    UInt max = 0;
    for (UInt i = 0; i < resFreqs.size(); i++)
    {
        if (resFreqs[i] > max)
        {
            max = resFreqs[i];
        }
    }

    //--determine relative chance of mutation acceptance
    do
    {
        chance = rand() % 100;
        mutant = _aminoacids[rand() % aaSize];
        if (count > 200)
        {
            acceptance = (resFreqs[mutant]/max)*100;
        }
        else
        {
            acceptance = 100;
        }
    }while (chance > acceptance);
    return mutant;
}

UInt getNumberofSequences()
{
    string inFrame;
    DIR *pdir;
    struct dirent *pent;
    pdir=opendir(".");
    UInt count = 0;

    //--get sequence evolution results for position
    while ((pent=readdir(pdir)))
    {
        inFrame = pent->d_name;
        if (inFrame.find(".evo.pdb") != std::string::npos)
        {
            count++;
        }
    }
    closedir(pdir);
    return count;
}
