//*******************************************************************************************************
//*******************************************************************************************************
//***********************************                          ******************************************
//***********************************      protFolder 1.0      ******************************************
//***********************************                          ******************************************
//*******************************************************************************************************
//********  -Sequence Selective Protein Genetic Algorithm Based Folding in Implicit Solvent -  **********
//*******************************************************************************************************

/////// Just specify infile structure, active chains and residues indexes, and it will evolve a fold

//--Included files and functions-------------------------------------------------------------------------
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <unistd.h>
#include "PDBInterface.h"

vector <UInt> getChainSequence(protein* _prot, UInt _chainIndex);
vector <UInt> getMutationPosition(UIntVec &_activeChains, UIntVec &_activeResidues);
UInt getProbabilisticMutation(vector < vector < UInt > > &_sequencePool, vector < vector < UInt > > &_possibleMutants, UIntVec &_mutantPosition, UIntVec &_activeResidues);
void createPossibleMutantsDatabase(protein* &_prot, UIntVec &_activeChains, UIntVec &_activeResidues, UIntVec &_allowedTypes);
UInt getSizeofPopulation();
vector < vector < UInt > > buildSequencePool();
vector < vector < UInt > > buildPossibleMutants();

enum structure {Z,M,C,L,P,B,E,Y,A,I,G,N,D,Q,R,F,H,W,K,S,T};
string backboneSeq[] =   {"", "M", "C", "L", "P", "B","E","Y","A","I","G",  "N",  "D",  "Q",  "R",  "F", "H", "W", "K", "S", "T"};
string backboneTypes[] = {"","-γ","-π","-α","-ρ","-β","β","ρ","α","π","γ","-γl","-πl","-αl","-ρl","-βl","βl","ρl","αl","πl","γl"};
UInt populationBaseline = 1000;

//--Program setup----------------------------------------------------------------------------------------
int main (int argc, char* argv[])
{
	//--Running parameters
	if (argc !=2)
	{
		cout << "protFolder <inFile.pdb>" << endl;
		exit(1);
	}

	UInt _activeChains[] = {0};                                                         // chains active for mutation
	UInt _allowedTypes[] = {M,C,L,P,B,E,Y,A,I,G,N,D,Q,R,F,H,W,K,S,T};                     // backbone types allowable
	UInt _activeResidues[] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19};                                     // positions active for mutation
	UInt _randomResidues[] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19};                                     // positions active for a random start sequence initially

	//--running parameters
	residue::setElectroSolvationScaleFactor(1.0);
	residue::setHydroSolvationScaleFactor(1.0);
	amberElec::setScaleFactor(1.0);
	amberVDW::setScaleFactor(1.0);
	residue::setTemperature(300);

	//convert input arrays to vectors
	UInt activeChainsSize = sizeof(_activeChains)/sizeof(_activeChains[0]), randomResiduesSize = sizeof(_randomResidues)/sizeof(_randomResidues[0]);
	UInt allowedTypesSize = sizeof(_allowedTypes)/sizeof(_allowedTypes[0]),activeResiduesSize = sizeof(_activeResidues)/sizeof(_activeResidues[0]);
	UIntVec activeChains, allowedTypes, activeResidues, randomResidues;
	for (UInt i = 0; i < activeChainsSize; i++)		{ activeChains.push_back(_activeChains[i]); }
	for (UInt i = 0; i < allowedTypesSize; i++)	{ allowedTypes.push_back(_allowedTypes[i]); }
	for (UInt i = 0; i < activeResiduesSize; i++)	{ activeResidues.push_back(_activeResidues[i]); }
	for (UInt i = 0; i < randomResiduesSize; i++)	{ randomResidues.push_back(_randomResidues[i]); }

	//--set initial variables
	srand (getpid());
	double startEnergy = 1E10, pastEnergy, Energy, Entropy, sPhi, sPsi, PiPj, KT = KB*residue::getTemperature();
	vector <double> backboneAngles(2);
	UInt timeid, sec, startNumClashes, numClashes, mutant = 0, plateau = 50, nobetter = 0;
	vector < UInt > mutantPosition, chainSequence, sequencePosition, randomPosition;
	vector < vector < UInt > > sequencePool, proteinSequence, finalSequence, possibleMutants;
	stringstream convert;
	string infile = argv[1];
	string startstr, outFile;
	UInt count, name = rand() % 100000000;
	convert << name, startstr = convert.str();
	string tempModel = startstr + "_temp.pdb";

	//-build possible sequence database per position
	PDBInterface* sthePDB = new PDBInterface(infile);
	ensemble* stheEnsemble = sthePDB->getEnsemblePointer();
	molecule* spMol = stheEnsemble->getMoleculePointer(0);
	protein* sprot = static_cast<protein*>(spMol);
	possibleMutants = buildPossibleMutants();
	if(possibleMutants.size() < activeResidues.size())
	{
		createPossibleMutantsDatabase(sprot, activeChains, activeResidues, allowedTypes);
		possibleMutants = buildPossibleMutants();
	}
	delete sthePDB;
	//--Run multiple independent evolution cycles-----------------------------------------------------
	while(true)
	{
		PDBInterface* thePDB = new PDBInterface(infile);
		ensemble* theEnsemble = thePDB->getEnsemblePointer();
		molecule* pMol = theEnsemble->getMoleculePointer(0);
		protein* prot = static_cast<protein*>(pMol);
		sequencePool = buildSequencePool();

		//--load in initial pdb and mutate in random starting structure on active chains and random residues
		nobetter = 0;
		for (UInt i = 0; i < activeChains.size(); i++)
		{
			for (UInt j = 0; j < randomResidues.size(); j++)
			{
				prot->activateForRepacking(activeChains[i], randomResidues[j]);
				randomPosition.push_back(activeChains[i]);
				randomPosition.push_back(randomResidues[j]);
				mutant = getProbabilisticMutation(sequencePool, possibleMutants, randomPosition, randomResidues);
				backboneAngles = prot->getRandPhiPsifromBackboneSequenceType(mutant);
				prot->setDihedral(activeChains[i], randomResidues[j],backboneAngles[0],0,0);
				prot->setDihedral(activeChains[i], randomResidues[j],backboneAngles[1],1,0);
				randomPosition.clear();
			}
			chainSequence = getChainSequence(prot, activeChains[i]);
			proteinSequence.push_back(chainSequence);
		}
		prot->protMin();
		startNumClashes = prot->getNumHardClashes();
		prot->setMoved(true);

		//--set Energy startpoint
		Energy = prot->protEnergy();

		//--Determine next mutation position
		mutantPosition.clear();
		mutantPosition = getMutationPosition(activeChains, activeResidues);
		pdbWriter(prot, tempModel);
		pastEnergy = Energy;

		//--Run through a single evolutionary path (ancestral line) till hitting plateau
		do
		{
			//--Mutate current sequence, new mutant and optimize system
			nobetter++;
			sequencePosition.push_back(mutantPosition[0]);
			sequencePosition.push_back(mutantPosition[1]);
			sPhi = prot->getPhi(mutantPosition[0],mutantPosition[1]);
			sPsi = prot->getPsi(mutantPosition[0],mutantPosition[1]);
			do{
				mutant = getProbabilisticMutation(sequencePool, possibleMutants, mutantPosition, activeResidues);
				backboneAngles = prot->getRandPhiPsifromBackboneSequenceType(mutant);
				prot->setDihedral(mutantPosition[0],mutantPosition[1],backboneAngles[0],0,0);
				prot->setDihedral(mutantPosition[0],mutantPosition[1],backboneAngles[1],1,0);
				prot->protRelax();
				numClashes = prot->getNumHardClashes();
				if (numClashes > startNumClashes){
					prot->setDihedral(mutantPosition[0],mutantPosition[1],sPhi,0,0);
					prot->setDihedral(mutantPosition[0],mutantPosition[1],sPsi,1,0);
				}
			}while (numClashes > startNumClashes);
			prot->protOpt(true);
			protein* tempProt = new protein(*prot);

			//--Determine next mutation position
			mutantPosition.clear();
			mutantPosition = getMutationPosition(activeChains, activeResidues);

			//--Energy test
			Energy = prot->protEnergy();
			Entropy = (1000000/((rand() % 1000000)+1))-1;
			PiPj = pow(EU,((Energy-pastEnergy)/KT));
			if (PiPj < Entropy){
				pdbWriter(tempProt, tempModel);
				proteinSequence[sequencePosition[0]][sequencePosition[1]] = mutant, pastEnergy = Energy;
				nobetter = 0;
			}
			else{
				prot->setDihedral(mutantPosition[0],mutantPosition[1],sPhi,0,0);
				prot->setDihedral(mutantPosition[0],mutantPosition[1],sPsi,1,0);
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
		
		//-Determine probability of being accepted into pool
		Energy = model->protEnergy();
		Entropy = (1000000/((rand() % 1000000)+1))-1;
		PiPj = pow(EU,((Energy-startEnergy)/KT));
		startEnergy = Energy;
		
		//-generate pdb output
		sec = time(NULL);
		timeid = sec;
		stringstream convert;
		string countstr;
		convert << timeid, countstr = convert.str();
		outFile = countstr + "." + startstr + ".fold.pdb";
		pdbWriter(model, outFile);
		
		//-write to data files
		count = getSizeofPopulation();
		finalSequence.clear(), chainSequence.clear();
		for (UInt i = 0; i < activeChains.size(); i++)
		{
			chainSequence = getChainSequence(model, activeChains[i]);
			finalSequence.push_back(chainSequence);
		}
		fstream finalline;
		finalline.open ("results.out", fstream::in | fstream::out | fstream::app);
		finalline << timeid << " " << Energy << " ";

		fstream fs;
		fs.open ("sequencepool.out", fstream::in | fstream::out | fstream::app);
		for (UInt i = 0; i < activeChains.size(); i++)
		{
			for (UInt j = 0; j < finalSequence[i].size(); j++)
			{
				finalline << backboneSeq[finalSequence[i][j]] << " ";
				if (PiPj < Entropy || count < ::populationBaseline){
					fs << finalSequence[i][j] << ",";
				}
			}
		}
		if (PiPj < Entropy || count < ::populationBaseline){
			finalline << "p";
			fs << endl;
		}
		fs.close();
		finalline << endl;
		finalline.close();
		
		//-clear variables for next iteration
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
	UInt type, numResidues;
	vector < UInt > sequence;
	numResidues = _prot->getNumResidues(_chainIndex);
	for (UInt j = 0; j < numResidues; j++)
	{
		type = _prot->getBackboneSequenceType(_chainIndex, j);
		sequence.push_back(type);
	}
	return sequence;
}

vector <UInt> getMutationPosition(UIntVec &_activeChains, UIntVec &_activeResidues)
{
	UInt randres, randchain;
	vector <UInt> _mutantPosition;
	randchain = _activeChains[rand() % _activeChains.size()];
	randres = _activeResidues[rand() % _activeResidues.size()];
	_mutantPosition.push_back(randchain);
	_mutantPosition.push_back(randres);
	return _mutantPosition;
}

UInt getProbabilisticMutation(vector < vector < UInt > > &_sequencePool, vector < vector < UInt > > &_possibleMutants, UIntVec &_mutantPosition, UIntVec &_activeResidues)
{
	double acceptance, threshold, FreqAccept;
	double poolSize = _sequencePool.size();
	vector <UInt> Freqs(21,1);
	UInt position, entropy, mutant, variance;
	UInt count = getSizeofPopulation();

	//--get sequence evolution results for position
	for (UInt i = 0; i < poolSize; i++)
	{
		UInt type = _sequencePool[i][_mutantPosition[1]];
		Freqs[type] = Freqs[type] + 1;
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
			entropy = 5;  // probabalistically allow 5% random genetic drift once sequence pool is sufficiently large
		}
		else{
			entropy = 100;  // 100% random sequences until sequence pool is built
		}
		if (variance > entropy) //control sequence entropy with probabilty
		{
			FreqAccept = Freqs[mutant];
			acceptance = (FreqAccept/(poolSize-1))*100;  //chance of accepting given fold at position is proportional to population
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
			stringstream seqString(item);
			int seqIndex;
			seqString >> seqIndex;
			sequence.push_back(seqIndex);
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
			stringstream seqString(item);
			int seqIndex;
			seqString >> seqIndex;
			_position.push_back(seqIndex);
		}
		_possibleMutants.push_back(_position);
		_position.clear();
	}
	file.close();
	return _possibleMutants;
}

void createPossibleMutantsDatabase(protein* &_prot, UIntVec &_activeChains, UIntVec &_activeResidues, UIntVec &_allowedTypes)
{
	fstream pm;
	pm.open ("possiblemutants.out", fstream::in | fstream::out | fstream::app);
	
	for (UInt i = 0; i < _activeChains.size(); i++)
	{
		for (UInt j = 0; j <_activeResidues.size(); j++)
		{
			for (UInt k = 0; k <_allowedTypes.size(); k++)
			{
				if (_allowedTypes[k] < 11 || (_allowedTypes[k] > 10 && _prot->getTypeFromResNum(_activeChains[i],_activeResidues[j]) > 25 )){
					pm << _allowedTypes[k] << ",";
				}
			}
			pm << endl;
		}
	}
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

