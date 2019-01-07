//*******************************************************************************************************
//*******************************************************************************************************
//***********************************                          ******************************************
//***********************************     protEvolver 3.0      ******************************************
//***********************************                          ******************************************
//*******************************************************************************************************
//********   -Fold Selective Protein Genetic Algorithm Based Evolution in Implicit Solvent -   **********
//*******************************************************************************************************

/////// Just specify infile structure, active chains and residues indexes, and it will evolve a sequence favorable for folding

//--Included files and functions-------------------------------------------------------------------------
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <unistd.h>
#include "PDBInterface.h"

vector <UInt> getChainSequence(protein* _prot, UInt _chainIndex);
vector <UInt> getMutationPosition(protein* _prot, UIntVec &_activeChains, UIntVec &_activeResidues);
UInt getProbabilisticMutation(vector < vector < UInt > > &_sequencePool, vector < vector < UInt > > &_possibleMutants, UIntVec &_mutantPosition, UIntVec &_activeResidues);
void createPossibleMutantsDatabase(protein* _prot, UIntVec &_activeChains, UIntVec &_activeResidues, UIntVec &_allowedLResidues, UIntVec &_allowedDResidues, bool _homosymmetric);
bool isFrozen(UIntVec _frozenResidues, UInt resIndex);
double calculatePopulationMA();
UInt getSizeofPopulation();
vector < vector < UInt > > buildSequencePool();
vector < vector < UInt > > buildPossibleMutants();

enum aminoAcid {A,R,N,D,Dh,C,Cx,Cf,Q,E,Eh,Hd,He,Hp,I,L,K,M,F,P,O,S,T,W,Y,V,G,dA,dR,dN,dD,dDh,dC,dCx,dCf,dQ,dE,dEh,dHd,dHe,dHp,dI,dL,dK,dM,dF,dP,dO,dS,dT,dW,dY,dV,Csf,Sf4,Hca,Eoc,Oec,Hem};
string aminoAcidString[] = {"A","R","N","D","Dh","C","Cx","Cf","Q","E","Eh","Hd","He","Hp","I","L","K","M","F","P","O","S","T","W","Y","V","G","dA","dR","dN","dD","dDh","dC","dCx","dCf","dQ","dE","dEh","dHd","dHe","dHp","dI","dL","dK","dM","dF","dP","dO","dS","dT","dW","dY","dV","Csf","Sf4","Hca","Eoc","Oec","Hem"};
UInt populationBaseline = 500;

//--Program setup----------------------------------------------------------------------------------------
int main (int argc, char* argv[])
{
	//--Running parameters
	if (argc !=2)
	{
		cout << "protEvolver <inFile.pdb>" << endl;
		exit(1);
	}

  	UInt _activeChains[] = {0};                                                         // chains active for mutation
    UInt _allowedLResidues[] = {A,R,N,D,Q,E,I,L,K,M,F,P,S,T,V,G};                     // amino acids allowed with phi < 0
    UInt _allowedDResidues[] = {A,R,N,D,Q,E,I,L,K,M,F,P,S,T,V,G};                                                     // amino acids allowed with phi > 0
    UInt _activeResidues[] = {0,1,2,3,4,5,7,8,9,10,12,15,16,17,19,21,22,23,24,26,28,29,30,31,33,34,35,36,37,38,39,40,41,42,44,46,47,48,49,51,53,54,55,56,58,60,61,62,63,65,67,68,69,70,71,72,73,74,76,77,78,79,80,81,83,84,85,86,88,91,92,93,95,97,98,99,100,102,103,104,105,106,107,109,110,111,112,113,114,115,116,117,118,120,122,123,124,125,127,129,130,131,132,134,136,137,138,139,141,143,144,145,146,147};                                     // positions active for mutation
    UInt _randomResidues[] = {0,1,2,3,4,5,7,8,9,10,12,15,16,17,19,21,22,23,24,26,28,29,30,31,33,34,35,36,37,38,39,40,41,42,44,46,47,48,49,51,53,54,55,56,58,60,61,62,63,65,67,68,69,70,71,72,73,74,76,77,78,79,80,81,83,84,85,86,88,91,92,93,95,97,98,99,100,102,103,104,105,106,107,109,110,111,112,113,114,115,116,117,118,120,122,123,124,125,127,129,130,131,132,134,136,137,138,139,141,143,144,145,146,147};                                     // positions active for a random start sequence initially
    UInt _frozenResidues[] = {6,11,13,14,18,20,25,27,32,43,45,50,52,57,59,64,66,75,82,87,89,90,94,96,101,108,119,121,126,128,133,135,140,142};                                  // positions that cannot move at all
    bool homoSymmetric = false;                                                          // if true all chains are structurally and sequentially symmetric to desired listed active chain above
    bool backboneRelaxation = false; 
    
	//--running parameters
	residue::setElectroSolvationScaleFactor(1.0);
	residue::setHydroSolvationScaleFactor(1.0);
	amberElec::setScaleFactor(1.0);
	amberVDW::setScaleFactor(1.0);

	//convert input arrays to vectors
	UInt activeChainsSize = sizeof(_activeChains)/sizeof(_activeChains[0]), randomResiduesSize = sizeof(_randomResidues)/sizeof(_randomResidues[0]), activeResiduesSize = sizeof(_activeResidues)/sizeof(_activeResidues[0]);
	UInt allowedLResiduesSize = sizeof(_allowedLResidues)/sizeof(_allowedLResidues[0]), allowedDResiduesSize = sizeof(_allowedDResidues)/sizeof(_allowedDResidues[0]), frozenResiduesSize = sizeof(_frozenResidues)/sizeof(_frozenResidues[0]);
	UIntVec activeChains, allowedLResidues, allowedDResidues, activeResidues, randomResidues, frozenResidues;
	for (UInt i = 0; i < activeChainsSize; i++)		{ activeChains.push_back(_activeChains[i]); }
	for (UInt i = 0; i < allowedLResiduesSize; i++)	{ allowedLResidues.push_back(_allowedLResidues[i]); }
	for (UInt i = 0; i < allowedDResiduesSize; i++)	{ allowedDResidues.push_back(_allowedDResidues[i]); }
	for (UInt i = 0; i < activeResiduesSize; i++)	{ activeResidues.push_back(_activeResidues[i]); }
	for (UInt i = 0; i < randomResiduesSize; i++)	{ randomResidues.push_back(_randomResidues[i]); }
	for (UInt i = 0; i < frozenResiduesSize; i++)	{ frozenResidues.push_back(_frozenResidues[i]); }

	//--set initial variables
	srand (getpid());
	double bestEnergy, pastEnergy, Energy;
	UInt timeid, sec, mutant = 0, numResidues, startingClashes, plateau = 10, nobetter = 0;
	vector < UInt > mutantPosition, chainSequence, sequencePosition, randomPosition;
	vector < vector < UInt > > sequencePool, proteinSequence, finalSequence, possibleMutants;
	stringstream convert;
	string infile = argv[1];
	string startstr, outFile;
	UInt name = rand() % 100000000;
	convert << name, startstr = convert.str();
	string tempModel = startstr + "_temp.pdb";


	//--determine which allowed amino acids are possible for each position from activeResidues and set cutoff energy
	PDBInterface* thePDB = new PDBInterface(infile);
	ensemble* theEnsemble = thePDB->getEnsemblePointer();
	molecule* pMol = theEnsemble->getMoleculePointer(0);
	protein* startProt = static_cast<protein*>(pMol);
	startingClashes = startProt->getNumHardClashes()*2;

	//--mutate all positions starting with a random resdiue to glycine
	for (UInt i = 0; i < activeChains.size(); i++)
	{
		for (UInt j = 0; j < randomResidues.size(); j++)
		{
			startProt->activateForRepacking(activeChains[i], randomResidues[j]);
			startProt->mutateWBC(activeChains[i], randomResidues[j], G);
		}
	}
	possibleMutants = buildPossibleMutants();
	if(possibleMutants.size() < activeResidues.size())
	{
		createPossibleMutantsDatabase(startProt, activeChains, activeResidues, allowedLResidues, allowedDResidues, homoSymmetric);
		possibleMutants = buildPossibleMutants();
	}
	delete thePDB;

	//--Run multiple independent evolution cycles-----------------------------------------------------
	for (UInt a = 1; a < 10000; a++)
	{
		PDBInterface* thePDB = new PDBInterface(infile);
		ensemble* theEnsemble = thePDB->getEnsemblePointer();
		molecule* pMol = theEnsemble->getMoleculePointer(0);
		protein* prot = static_cast<protein*>(pMol);
		sequencePool = buildSequencePool();
		if (homoSymmetric)
		{
			for (UInt i = 1; i < prot->getNumChains(); i++)
			{
				prot->symmetryLinkChainAtoB(i, activeChains[0]);
			}
		}

		//--load in initial pdb and mutate in random starting sequence on active chains and random residues
		nobetter = 0;
		for (UInt i = 0; i < activeChains.size(); i++)
		{
			for (UInt j = 0; j < randomResidues.size(); j++)
			{
				prot->activateForRepacking(activeChains[i], randomResidues[j]);
				randomPosition.push_back(activeChains[i]);
				randomPosition.push_back(randomResidues[j]);
				mutant = getProbabilisticMutation(sequencePool, possibleMutants, randomPosition, randomResidues);
				prot->mutateWBC(activeChains[i], randomResidues[j], mutant);
				randomPosition.clear();
			}
			chainSequence = getChainSequence(prot, activeChains[i]);
			proteinSequence.push_back(chainSequence);
		}
		prot->protOpt(backboneRelaxation, frozenResidues, activeChains);

		//--set Energy startpoint
		Energy = prot->protEnergy();

		//--Determine next mutation position
		mutantPosition.clear();
		mutantPosition = getMutationPosition(prot, activeChains, activeResidues);
		pdbWriter(prot, tempModel);
		pastEnergy = Energy;
		bestEnergy = Energy;

		//--Run through a single evolutionary path (ancestral line) till hitting plateau
		do
		{
			//--Mutate current sequence, new mutant and optimize system
			nobetter++;
			for (UInt i = 0; i < activeChains.size(); i++)
			{
				numResidues = prot->getNumResidues(activeChains[i]);
				for (UInt j = 0; j < numResidues; j++)
				{
					prot->activateForRepacking(activeChains[i],j);
					if (activeChains[i] == mutantPosition[0] && j == mutantPosition[1])
					{
						//--new mutant
						sequencePosition.push_back(i);
						sequencePosition.push_back(j);
						mutant = getProbabilisticMutation(sequencePool, possibleMutants, mutantPosition, activeResidues);
						prot->mutateWBC(mutantPosition[0],mutantPosition[1], mutant);
					}
					bool frozen = isFrozen(frozenResidues,j);
					if (!frozen)
					{
						prot->mutateWBC(activeChains[i],j, proteinSequence[i][j]);
					}
				}
			}
			prot->protOpt(backboneRelaxation, frozenResidues, activeChains);
			protein* tempProt = new protein(*prot);

			//--Determine next mutation position
			mutantPosition.clear();
			mutantPosition = getMutationPosition(prot, activeChains, activeResidues);

			//--Energy test
			Energy = prot->protEnergy();
			if (Energy < pastEnergy)
			{
				if (Energy < bestEnergy)
				{
					bestEnergy = Energy;
					pdbWriter(tempProt, tempModel);
				}
				proteinSequence[sequencePosition[0]][sequencePosition[1]] = mutant, pastEnergy = Energy;
				if (nobetter > 0) { nobetter--;}
				else{ nobetter = 0; }
			}
			sequencePosition.clear();
			delete tempProt;
		}while (nobetter < plateau);
		delete thePDB;

		//--Print final energy and write a pdb file----------------------------------------------------
		PDBInterface* theModelPDB = new PDBInterface(tempModel);
		ensemble* theModelEnsemble = theModelPDB->getEnsemblePointer();
		molecule* modelMol = theModelEnsemble->getMoleculePointer(0);
		protein* model = static_cast<protein*>(modelMol);
		UInt clashes = model->getNumHardClashes();
		if (clashes <= startingClashes)
		{
			model->setMoved(true);
			Energy = model->protEnergy();
			sec = time(NULL);
			timeid = sec;
			stringstream convert;
			string countstr;
			convert << timeid, countstr = convert.str();
			outFile = countstr + "." + startstr + ".evo.pdb";
			pdbWriter(model, outFile);
			finalSequence.clear(), chainSequence.clear();
			for (UInt i = 0; i < activeChains.size(); i++)
			{
				chainSequence = getChainSequence(model, activeChains[i]);
				finalSequence.push_back(chainSequence);
			}
			fstream finalline;
			finalline.open ("results.out", fstream::in | fstream::out | fstream::app);
			finalline << timeid << " " << Energy << " ";
	
			double popMa = calculatePopulationMA();
			fstream fs;
			fs.open ("sequencepool.out", fstream::in | fstream::out | fstream::app);
			for (UInt i = 0; i < activeChains.size(); i++)
			{
				for (UInt j = 0; j < finalSequence[i].size(); j++)
				{
					finalline << aminoAcidString[finalSequence[i][j]] << " ";
					if (Energy < popMa)
					{
						fs << finalSequence[i][j] << ",";
					}
				}
			}
			if (Energy < popMa){fs << endl;}
			fs.close();
			finalline << endl;
			finalline.close();
		}
		delete theModelPDB;
		sequencePool.clear(),proteinSequence.clear(), chainSequence.clear(), mutantPosition.clear(), chainSequence.clear(), sequencePosition.clear(), randomPosition.clear();
		sequencePool.resize(0),proteinSequence.resize(0), chainSequence.resize(0), mutantPosition.resize(0), chainSequence.resize(0), sequencePosition.resize(0), randomPosition.resize(0);
	}
	return 0;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////
//////// functions //////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////

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

vector <UInt> getMutationPosition(protein* _prot, UIntVec &_activeChains, UIntVec &_activeResidues)
{
	//--get median residue energy
	UInt randres, randchain;
	double posE, medE;
	vector <UInt> _mutantPosition;
	medE = _prot->getMedianResidueEnergy(_activeChains, _activeResidues);

	//--find random position with worse than median energy
	do
	{
		randchain = _activeChains[rand() % _activeChains.size()];
		randres = _activeResidues[rand() % _activeResidues.size()];
		posE = _prot->protEnergy(randchain,randres);
	}while (posE < medE);
	_mutantPosition.push_back(randchain);
	_mutantPosition.push_back(randres);
	return _mutantPosition;
}

UInt getProbabilisticMutation(vector < vector < UInt > > &_sequencePool, vector < vector < UInt > > &_possibleMutants, UIntVec &_mutantPosition, UIntVec &_activeResidues)
{
	double acceptance, threshold, resFreqAccept;
	double poolSize = _sequencePool.size();
	vector <UInt> resFreqs(58,1);
	UInt position, entropy, mutant, variance;
	UInt count = getSizeofPopulation();

	//--get sequence evolution results for position
	for (UInt i = 0; i < poolSize; i++)
	{
		UInt restype = _sequencePool[i][_mutantPosition[1]];
		resFreqs[restype] = resFreqs[restype] + 1;
	}

	//--find mutant in vector
	for (UInt i = 0; i < _possibleMutants.size(); i++)
	{
		if (_activeResidues[i] == _mutantPosition[1])
		{
			position = i;
		}
	}
	UInt positionPossibles = _possibleMutants[position].size();

	//--determine population based chance of mutation acceptance or a random mutation
	do
	{
		threshold = (rand() % 100) + 1;
		variance = (rand() % 100) + 1;
		mutant = _possibleMutants[position][rand() % positionPossibles];
		if (count >= ::populationBaseline){
			entropy = 5;  // probabalistically allow 33% random genetic drift once sequence pool is sufficiently large
		}
		else{
			entropy = 100;  // 100% random sequences until sequence pool is built
		}
		if (variance > entropy) //control sequence entropy with probabilty
		{
			resFreqAccept = resFreqs[mutant];
			acceptance = (resFreqAccept/(poolSize-1))*100;  //chance of accepting given amino acid at position is proportional to population
		}
		else
		{
			acceptance = 100;  //random mutation
		}
	}while (threshold > acceptance);
	return mutant;
}

vector < vector < UInt > > buildSequencePool()
{
	ifstream file("sequencepool.out");
	string item, line;
	vector < UInt > sequence;
	vector < vector < UInt > > sequencePool;

	while(getline(file,line))
	{
		stringstream stream(line);
		while(getline(stream,item,','))
		{
			stringstream aaString(item);
			int aaIndex;
			aaString >> aaIndex;
			sequence.push_back(aaIndex);
		}
		sequencePool.push_back(sequence);
		sequence.clear();
	}
	file.close();
	if (sequencePool.size() > ::populationBaseline){
		sequencePool.erase(sequencePool.begin(),sequencePool.end()-::populationBaseline);
	}
	return sequencePool;
}

vector < vector < UInt > > buildPossibleMutants()
{
	ifstream file("possiblemutants.out");
	string item, line;
	vector < UInt > _position;
	vector < vector < UInt > > _possibleMutants;
	while(getline(file,line))
	{
		stringstream stream(line);
		while(getline(stream,item,','))
		{
			stringstream aaString(item);
			int aaIndex;
			aaString >> aaIndex;
			_position.push_back(aaIndex);
		}
		_possibleMutants.push_back(_position);
		_position.clear();
	}
	file.close();
	return _possibleMutants;
}

void createPossibleMutantsDatabase(protein* _prot, UIntVec &_activeChains, UIntVec &_activeResidues, UIntVec &_allowedLResidues, UIntVec &_allowedDResidues, bool _homoSymmetric)
{
	double Energy, phi, startE;
	UInt restype;
	bool added;
	fstream pm;
	pm.open ("possiblemutants.out", fstream::in | fstream::out | fstream::app);
	if (_homoSymmetric)
	{
		for (UInt i = 1; i < _prot->getNumChains(); i++)
		{
			_prot->symmetryLinkChainAtoB(i, _activeChains[0]);
		}
	}

	startE = _prot->protEnergy();
	startE=startE+20; //buffer energy filter of amino acids per position to restrict only by very hard clashes

	for (UInt i = 0; i < _activeChains.size(); i++)
	{
		for (UInt j = 0; j <_activeResidues.size(); j++)
		{
			phi = _prot->getPhi(_activeChains[i], _activeResidues[j]);
			restype = _prot->getTypeFromResNum(_activeChains[i], _activeResidues[j]);
			if ((phi < 0 && phi > -180) || _activeResidues[j] == 0)
			{
				for (UInt k = 0; k <_allowedLResidues.size(); k++)
				{
					added = false;
					_prot->activateForRepacking(_activeChains[i], _activeResidues[j]);
					_prot->mutateWBC(_activeChains[i], _activeResidues[j], _allowedLResidues[k]);
					UIntVec allowedRots = _prot->getAllowedRotamers(_activeChains[i], _activeResidues[j], _allowedLResidues[k], 0);
					if (allowedRots.size() > 0)
					{
						for (UInt l = 0; l < allowedRots.size(); l++)
						{
							_prot->setRotamerWBC(_activeChains[i], _activeResidues[j], 0, allowedRots[l]);
							Energy = _prot->protEnergy();
							if (Energy <= startE)
							{
								if (_allowedLResidues[k] != restype && !added)
								{
									pm << _allowedLResidues[k] << ",";
									added = true;
								}
							}
							else
							{
								_prot->setRotamerNotAllowed(_activeChains[i], _activeResidues[j], _allowedLResidues[k], 0, allowedRots[l]);
							}
						}
					}
					else
					{
						Energy = _prot->protEnergy();
						if (Energy <= startE)
						{
							if (_allowedLResidues[k] != restype && !added)
							{
								pm << _allowedLResidues[k] << ",";
								added = true;
							}
						}
					}
				}
			}
			if (phi > 0 && phi < 180)
			{
				for (UInt k = 0; k <_allowedDResidues.size(); k++)
				{
					added = false;
					_prot->activateForRepacking(_activeChains[i], _activeResidues[j]);
					_prot->mutateWBC(_activeChains[i], _activeResidues[j], _allowedDResidues[k]);
					UIntVec allowedRots = _prot->getAllowedRotamers(_activeChains[i], _activeResidues[j], _allowedDResidues[k], 0);
					if (allowedRots.size() > 0)
					{
						for (UInt l = 0; l < allowedRots.size(); l++)
						{
							_prot->setRotamerWBC(_activeChains[i], _activeResidues[j], 0, allowedRots[l]);
							Energy = _prot->protEnergy();
							if (Energy <= startE)
							{
								if (_allowedDResidues[k] != restype && !added)
								{
									pm << _allowedDResidues[k] << ",";
									added = true;
								}
							}
							else
							{
								_prot->setRotamerNotAllowed(_activeChains[i], _activeResidues[j], _allowedDResidues[k], 0, allowedRots[l]);
							}
						}
					}
					else
					{
						Energy = _prot->protEnergy();
						if (Energy <= startE)
						{
							if (_allowedDResidues[k] != restype && !added)
							{
								pm << _allowedDResidues[k] << ",";
								added = true;
							}
						}
					}
				}
			}
			_prot->mutateWBC(_activeChains[i], _activeResidues[j], restype);
			pm << restype << ",";
			pm << endl;
		}
	}
}

bool isFrozen(UIntVec _frozenResidues, UInt resIndex)
{
	bool frozen = false;
	for (UInt i = 0; i < _frozenResidues.size(); i++)
	{
		if (_frozenResidues[i] == resIndex)
		{
			frozen = true;
		}
	}
	return frozen;
}

double calculatePopulationMA()
{
	ifstream file("results.out");
	string item, line;
	bool secondSpace;
	vector < double > _energy;
	while(getline(file,line))
	{
		secondSpace = false;
		stringstream stream(line);
		while(getline(stream,item,' '))
		{
			if (secondSpace)
			{
				stringstream energyString(item);
				double energy;
				energyString >> energy;
				_energy.push_back(energy);
				break;
			}
			else
			{
				secondSpace = true;
			}
		}
	}
	file.close();
	double cutoff = 0.0;
	if (_energy.size() >= ::populationBaseline){
		double sum = 0.0;
		for (UInt i = _energy.size()-100; i < _energy.size(); i++)
		{
			sum += _energy[i];
		}
		cutoff=sum/100;
	}
	else{
		cutoff = 1E10;
	}
	return cutoff;
}

UInt getSizeofPopulation()
{
	ifstream file("results.out");
	string line;
	UInt counter = 0;
	while(getline(file,line))
	{
		counter++;
	}
	file.close();
	return counter;
}

