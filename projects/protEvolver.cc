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
#include <time.h>
#include <dirent.h>
#include <sstream>
#include <fstream>
#include <unistd.h>
#include "ensemble.h"
#include "PDBInterface.h"

vector <UInt> getChainSequence(protein* _prot, UInt _chainIndex);
vector <UInt> getMutationPosition(protein* _prot, UIntVec &_activeChains, UIntVec &_activeResidues);
UInt getProbabilisticMutation(vector < vector < UInt > > &_sequencePool, vector < vector < UInt > > &_possibleMutants, UIntVec &_mutantPosition, UIntVec &_activeResidues);
void createPossibleMutantsDatabase(protein* bundle, UIntVec &_activeChains, UIntVec &_activeResidues, UIntVec &_allowedLResidues, UIntVec &_allowedDResidues, bool _homosymmetric);
bool isFrozen(UIntVec _frozenResidues, UInt resIndex);
vector < vector < UInt > > buildSequencePool();
vector < vector < UInt > > buildPossibleMutants();
enum aminoAcid {A,R,N,D,Dh,C,Cx,Cf,Q,E,Eh,Hd,He,Hp,I,L,K,M,F,P,O,S,T,W,Y,V,G,dA,dR,dN,dD,dDh,dC,dCx,dCf,dQ,dE,dEh,dHd,dHe,dHp,dI,dL,dK,dM,dF,dP,dO,dS,dT,dW,dY,dV,Csf,Sf4,Hca,Eoc,Oec,Hem};
string aminoAcidString[] = {"A","R","N","D","Dh","C","Cx","Cf","Q","E","Eh","Hd","He","Hp","I","L","K","M","F","P","O","S","T","W","Y","V","G","dA","dR","dN","dD","dDh","dC","dCx","dCf","dQ","dE","dEh","dHd","dHe","dHp","dI","dL","dK","dM","dF","dP","dO","dS","dT","dW","dY","dV","Csf","Sf4","Hca","Eoc","Oec","Hem"};

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
    UInt _allowedLResidues[] = {A,R,N,D,Q,E,I,L,K,F,P,S,T,W,Y,V,G};                   // amino acids allowed with phi < 0
	UInt _allowedDResidues[] = {G};  // amino acids allowed with phi > 0
    UInt _activeResidues[] = {6,8,13,14};                  // positions active for mutation
    UInt _randomResidues[] = {6,8,13,14};                  // positions active for a random start sequence initially
    UInt _frozenResidues[] = {7,10,15,43,19,54,59,60}; //13,17,61,65                                              // positions that cannot move at all
	bool homoSymmetric = false;                                                          // if true all chains are structurally symmetrical to the one listed active chain above
	bool backboneRelaxation = false;                                             // if true allow backrub relaxation in structural optimization

	//--running parameters
	residue::setCutoffDistance(8.0);
	residue::setTemperature(300);
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
	double bestEnergy, pastEnergy, Energy, startEnergy;
	UInt timeid, sec, mutant = 0, numResidues, plateau = 10, nobetter = 0;
	vector < UInt > mutantPosition, chainSequence, sequencePosition, randomPosition;
	vector < vector < UInt > > sequencePool, proteinSequence, finalSequence, possibleMutants;
	vector < double > bindingEnergy;
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
	protein* start_bundle = static_cast<protein*>(pMol);
	if (homoSymmetric)
	{
		startEnergy = start_bundle->intraSoluteEnergy(true, _activeChains[0]);
	}
	else
	{
		startEnergy = start_bundle->intraSoluteEnergy(true);
	}
	possibleMutants = buildPossibleMutants();
	if(possibleMutants.size() < activeResidues.size())
	{
		createPossibleMutantsDatabase(start_bundle, activeChains, activeResidues, allowedLResidues, allowedDResidues, homoSymmetric);
		possibleMutants = buildPossibleMutants();
	}

	//--Run multiple independent evolution cycles-----------------------------------------------------
	for (UInt a = 1; a < 10000; a++)
	{
		sequencePool = buildSequencePool();
		if (homoSymmetric)
		{
			for (UInt i = 1; i < start_bundle->getNumChains(); i++)
			{
				start_bundle->symmetryLinkChainAtoB(i, activeChains[0]);
			}
		}

		//--load in initial pdb and mutate in random starting sequence on active chains and random residues
		nobetter = 0;
		for (UInt i = 0; i < activeChains.size(); i++)
		{
			for (UInt j = 0; j < randomResidues.size(); j++)
			{
				start_bundle->activateForRepacking(activeChains[i], randomResidues[j]);
				randomPosition.push_back(activeChains[i]);
				randomPosition.push_back(randomResidues[j]);
				mutant = getProbabilisticMutation(sequencePool, possibleMutants, randomPosition, randomResidues);
				start_bundle->mutateWBC(activeChains[i], randomResidues[j], mutant);
				randomPosition.clear();
			}
			chainSequence = getChainSequence(start_bundle, activeChains[i]);
			proteinSequence.push_back(chainSequence);
		}

		start_bundle->protOpt(backboneRelaxation, frozenResidues, activeChains);

		//--set Energy startpoint
		if (homoSymmetric)
		{
			Energy = start_bundle->intraSoluteEnergy(true, activeChains[0]);
		}
		else
		{
			Energy = start_bundle->protEnergy();
		}
		//--Determine next mutation position
		mutantPosition.clear();
		mutantPosition = getMutationPosition(start_bundle, activeChains, activeResidues);
		pdbWriter(start_bundle, tempModel);
		pastEnergy = Energy;
		bestEnergy = Energy;

		//--Run through a single evolutionary path (ancestral line) till hitting plateau
		do
		{
			//--Mutate current sequence, new mutant and optimize system
			nobetter++;
			for (UInt i = 0; i < activeChains.size(); i++)
			{
				numResidues = start_bundle->getNumResidues(activeChains[i]);
				for (UInt j = 0; j < numResidues; j++)
				{
					start_bundle->activateForRepacking(activeChains[i],j);
					if (activeChains[i] == mutantPosition[0] && j == mutantPosition[1])
					{
						//--new mutant
						sequencePosition.push_back(i);
						sequencePosition.push_back(j);
						mutant = getProbabilisticMutation(sequencePool, possibleMutants, mutantPosition, activeResidues);
						start_bundle->mutateWBC(mutantPosition[0],mutantPosition[1], mutant);
					}
					bool frozen = isFrozen(frozenResidues,j);
					if (!frozen)
					{
						start_bundle->mutateWBC(activeChains[i],j, proteinSequence[i][j]);
					}
				}
			}
			start_bundle->protOpt(backboneRelaxation, frozenResidues, activeChains);
			protein* tempBundle = new protein(*start_bundle);

			//--Determine next mutation position
			mutantPosition.clear();
			mutantPosition = getMutationPosition(start_bundle, activeChains, activeResidues);

			//--Energy test
			if (homoSymmetric)
			{
				Energy = start_bundle->intraSoluteEnergy(true, activeChains[0]);
			}
			else
			{
				Energy = start_bundle->protEnergy();
			}
			if (Energy < pastEnergy)
			{
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
			delete tempBundle;
		}while (nobetter < plateau);

		//--Print final energy and write a pdb file----------------------------------------------------
		PDBInterface* theModelPDB = new PDBInterface(tempModel);
		ensemble* theModelEnsemble = theModelPDB->getEnsemblePointer();
		molecule* modelMol = theModelEnsemble->getMoleculePointer(0);
		protein* model = static_cast<protein*>(modelMol);
		if (homoSymmetric)
		{
			Energy = model->intraSoluteEnergy(true, activeChains[0]);
		}
		else
		{
			Energy = model->protEnergy();
		}
		bindingEnergy.clear();
		bindingEnergy = model->chainBindingEnergy();
        if (Energy < startEnergy)
		{
			name = rand() % 100;
			sec = time(NULL);
			timeid = name + sec;
			stringstream convert;
			string countstr;
			convert << timeid, countstr = convert.str();
			outFile = countstr + ".evo.pdb";
			pdbWriter(model, outFile);
			finalSequence.clear(), chainSequence.clear();
			for (UInt i = 0; i < activeChains.size(); i++)
			{
				chainSequence = getChainSequence(model, activeChains[i]);
				finalSequence.push_back(chainSequence);
			}
			fstream finalline;
			finalline.open ("results.out", fstream::in | fstream::out | fstream::app);
			finalline << timeid << " " << bindingEnergy[0] << " " << bindingEnergy[1] << " ";

			fstream fs;
			fs.open ("sequencepool.out", fstream::in | fstream::out | fstream::app);
			for (UInt i = 0; i < activeChains.size(); i++)
			{
				for (UInt j = 0; j < finalSequence[i].size(); j++)
				{
					finalline << aminoAcidString[finalSequence[i][j]] << " ";
					fs << finalSequence[i][j] << ",";
				}
			}
			fs << endl;
			finalline << endl;
			finalline.close();
			fs.close();
		}
		delete theModelPDB;
		bindingEnergy.clear(),sequencePool.clear(),proteinSequence.clear(), chainSequence.clear(), mutantPosition.clear(), chainSequence.clear(), sequencePosition.clear(), randomPosition.clear();
		bindingEnergy.resize(0),sequencePool.resize(0),proteinSequence.resize(0), chainSequence.resize(0), mutantPosition.resize(0), chainSequence.resize(0), sequencePosition.resize(0), randomPosition.resize(0);
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
	medE = _prot->getMedianResEnergy(_activeChains, _activeResidues);

	//--find random position with worse than median energy
	do
	{
		randchain = _activeChains[rand() % _activeChains.size()];
		randres = _activeResidues[rand() % _activeResidues.size()];
		posE = _prot->resEnergy(randchain,randres);
	}while (posE < medE);
	_mutantPosition.push_back(randchain);
	_mutantPosition.push_back(randres);
	return _mutantPosition;
}

UInt getProbabilisticMutation(vector < vector < UInt > > &_sequencePool, vector < vector < UInt > > &_possibleMutants, UIntVec &_mutantPosition, UIntVec &_activeResidues)
{
	float mutant, chance, entropy, acceptance, pooling, resFreqAccept;
	vector <UInt> resFreqs(58,1);
	UInt position;
	float count = _sequencePool.size();

	//--get sequence evolution results for position
	for (UInt i = 0; i < _sequencePool.size(); i++)
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

	//--determine population based chance of mutation acceptance or a random mutation, via linear regression of sequence entropy
	do
	{
		chance = (rand() % 100) + 1;
		entropy = (rand() % 100) + 1;
		mutant = _possibleMutants[position][rand() % positionPossibles];
		pooling = 5;//-0.3166*count+195; //At 300 sequences start pooling, at 600+ 95%
		if (pooling < 5){
			pooling = 5;
		}
		if (entropy > pooling)  //sequence entropy determined by pooling linear decline to resolve minima after suitable diversity
		{
			resFreqAccept = resFreqs[mutant];
			acceptance = (resFreqAccept/(count-1))*100; //chance of accepting given amino acid at position is proportional to population
		}
		else
		{
			acceptance = 100;  //random mutation
		}
	}while (chance > acceptance);
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

void createPossibleMutantsDatabase(protein* _bundle, UIntVec &_activeChains, UIntVec &_activeResidues, UIntVec &_allowedLResidues, UIntVec &_allowedDResidues, bool _homoSymmetric)
{
	double Energy, phi, startE, totEnergy, refEnergy;
	UInt restype;
	bool added;
	fstream pm;
	pm.open ("possiblemutants.out", fstream::in | fstream::out | fstream::app);
	if (_homoSymmetric)
	{
		for (UInt i = 1; i < _bundle->getNumChains(); i++)
		{
			_bundle->symmetryLinkChainAtoB(i, _activeChains[0]);
		}
		startE = _bundle->intraSoluteEnergy(true, _activeChains[0]);
	}
	else
	{
		startE = _bundle->protEnergy();
	}
    startE=startE+20; //buffer energy filter of amino acids per position to restrict only by very hard clashes

	for (UInt i = 0; i < _activeChains.size(); i++)
	{
		for (UInt j = 0; j <_activeResidues.size(); j++)
		{
			phi = _bundle->getPhi(_activeChains[i], _activeResidues[j]);
			restype = _bundle->getTypeFromResNum(_activeChains[i], _activeResidues[j]);
			if ((phi < 0 && phi > -180) || _activeResidues[j] == 0)
			{
				for (UInt k = 0; k <_allowedLResidues.size(); k++)
				{
					added = false;
					_bundle->activateForRepacking(_activeChains[i], _activeResidues[j]);
					_bundle->mutateWBC(_activeChains[i], _activeResidues[j], _allowedLResidues[k]);
					UIntVec allowedRots = _bundle->getAllowedRotamers(_activeChains[i], _activeResidues[j], _allowedLResidues[k], 0);
					if (allowedRots.size() > 0)
					{
						for (UInt l = 0; l < allowedRots.size(); l++)
						{
							_bundle->setRotamerWBC(_activeChains[i], _activeResidues[j], 0, allowedRots[l]);
							if (_homoSymmetric)
							{
								totEnergy = _bundle->intraSoluteEnergy(true, _activeChains[0]);
								refEnergy = _bundle->getFreeAminoAcidEnergy(_activeChains[i],_activeResidues[j]);
								Energy = totEnergy-refEnergy;
							}
							else
							{
								Energy = _bundle->protEnergy();
							}
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
								_bundle->setRotamerNotAllowed(_activeChains[i], _activeResidues[j], _allowedLResidues[k], 0, allowedRots[l]);
							}
						}
					}
					else
					{
						if (_homoSymmetric)
						{
							totEnergy = _bundle->intraSoluteEnergy(true, _activeChains[0]);
							refEnergy = _bundle->getFreeAminoAcidEnergy(_activeChains[i],_activeResidues[j]);
							Energy = totEnergy-refEnergy;
						}
						else
						{
							Energy = _bundle->protEnergy();
						}
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
					_bundle->activateForRepacking(_activeChains[i], _activeResidues[j]);
					_bundle->mutateWBC(_activeChains[i], _activeResidues[j], _allowedDResidues[k]);
					UIntVec allowedRots = _bundle->getAllowedRotamers(_activeChains[i], _activeResidues[j], _allowedDResidues[k], 0);
					if (allowedRots.size() > 0)
					{
						for (UInt l = 0; l < allowedRots.size(); l++)
						{
							_bundle->setRotamerWBC(_activeChains[i], _activeResidues[j], 0, allowedRots[l]);
							if (_homoSymmetric)
							{
								totEnergy = _bundle->intraSoluteEnergy(true, _activeChains[0]);
								refEnergy = _bundle->getFreeAminoAcidEnergy(_activeChains[i],_activeResidues[j]);
								Energy = totEnergy-refEnergy;
							}
							else
							{
								Energy = _bundle->protEnergy();
							}
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
								_bundle->setRotamerNotAllowed(_activeChains[i], _activeResidues[j], _allowedDResidues[k], 0, allowedRots[l]);
							}
						}
					}
					else
					{
						if (_homoSymmetric)
						{
							totEnergy = _bundle->intraSoluteEnergy(true, _activeChains[0]);
							refEnergy = _bundle->getFreeAminoAcidEnergy(_activeChains[i],_activeResidues[j]);
							Energy = totEnergy-refEnergy;
						}
						else
						{
							Energy = _bundle->protEnergy();
						}
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
			_bundle->mutateWBC(_activeChains[i], _activeResidues[j], restype);
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
