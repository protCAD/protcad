//*******************************************************************************************************
//*******************************************************************************************************
//***********************************                          ******************************************
//***********************************       protFolder 1       ******************************************
//***********************************                          ******************************************
//*******************************************************************************************************
//********  -Confirmation Selective Protein Genetic Algorithm Based Folding in Implicit Solvent -  **********
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
UInt getProbabilisticMutation(protein* _prot, vector < vector < UInt > > &_sequencePool, vector < vector < UInt > > &_possibleMutants, UIntVec &_mutantPosition, UIntVec &_activeResidues);
void createPossibleMutantsDatabase(protein* &_prot, UIntVec &_activeChains, UIntVec &_activeResidues, UIntVec &_allowedTypes);
UInt getSizeofPopulation();
vector < vector < UInt > > readSequencePool();
vector < vector < UInt > > readPossibleMutants();

enum structure {Z,M,C,L,P,B,E,Y,A,I,G,N,D,Q,R,F,H,W,K,S,T};
string backboneSeq[] =   {"", "M", "C", "L", "P", "B","E","Y","A","I","G",  "N",  "D",  "Q",  "R",  "F", "H", "W", "K", "S", "T"};
string backboneTypes[] = {"","-γ","-π","-α","-ρ","-β","β","ρ","α","π","γ","-γl","-πl","-αl","-ρl","-βl","βl","ρl","αl","πl","γl"};
UInt populationBaseline = 1000;
double KT = KB*residue::getTemperature();

//--Program setup----------------------------------------------------------------------------------------
int main (int argc, char* argv[])
{
	//--Running parameters
	if (argc !=2)
	{
		cout << "protFolder <inFile.pdb>" << endl;
		exit(1);
	}

	//--input
	UInt _activeChains[] = {0};                                                         // chains active for mutation
	UInt _allowedTypes[] = {C,L,P,B,E,Y,A,I,D,Q,R,F,H,W,K,S};                     // backbone types allowable
	UInt _activeResidues[] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19};                                     // positions active for mutation
	UInt _randomResidues[] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19};                                     // positions active for a random start sequence initially

	//--running parameters
	residue::setElectroSolvationScaleFactor(1.0);
	residue::setHydroSolvationScaleFactor(1.0);
	residue::setPolarizableElec(true);
	amberElec::setScaleFactor(1.0);
	amberVDW::setScaleFactor(1.0);

	//convert input arrays to vectors
	UInt activeChainsSize = sizeof(_activeChains)/sizeof(_activeChains[0]), randomResiduesSize = sizeof(_randomResidues)/sizeof(_randomResidues[0]);
	UInt allowedTypesSize = sizeof(_allowedTypes)/sizeof(_allowedTypes[0]),activeResiduesSize = sizeof(_activeResidues)/sizeof(_activeResidues[0]);
	UIntVec activeChains, allowedTypes, activeResidues, randomResidues;
	for (UInt i = 0; i < activeChainsSize; i++)		{ activeChains.push_back(_activeChains[i]); }
	for (UInt i = 0; i < allowedTypesSize; i++)	{ allowedTypes.push_back(_allowedTypes[i]); }
	for (UInt i = 0; i < activeResiduesSize; i++)	{ activeResidues.push_back(_activeResidues[i]); }
	for (UInt i = 0; i < randomResiduesSize; i++)	{ randomResidues.push_back(_randomResidues[i]); }

	//--set initial variables
	int seed = (int)getpid()*(int)gethostid();
	srand (seed);
	double startEnergy = 1E10, pastEnergy, Energy, deltaEnergy;
	vector <double> backboneAngles(2);
	UInt timeid, sec, numClashes, startNumBackboneClashes, mutant = 0, plateau = 1000, nobetter = 0;
	vector < UInt > mutantPosition, chainSequence, randomPosition;
	vector < vector < UInt > > sequencePool, finalSequence, possibleMutants;
	stringstream convert;
	string infile = argv[1];
	string startstr, outFile;
	UInt count, name = rand() % 100000000;
	convert << name, startstr = convert.str();
	string tempModel = startstr + "_temp.pdb";
	bool revert, boltzmannAcceptance;

	//-build possible fold sequence database per position
	PDBInterface* thePDB = new PDBInterface(infile);
	ensemble* theEnsemble = thePDB->getEnsemblePointer();
	molecule* pMol = theEnsemble->getMoleculePointer(0);
	protein* prot = static_cast<protein*>(pMol);
	possibleMutants = readPossibleMutants();
	if(possibleMutants.size() < activeResidues.size())
	{
		createPossibleMutantsDatabase(prot, activeChains, activeResidues, allowedTypes);
		possibleMutants = readPossibleMutants();
	}
	delete thePDB;
	
	//--Run multiple independent fold cycles-----------------------------------------------------
	while(true)
	{
		thePDB = new PDBInterface(infile);
		theEnsemble = thePDB->getEnsemblePointer();
		pMol = theEnsemble->getMoleculePointer(0);
		prot = static_cast<protein*>(pMol);
		sequencePool = readSequencePool();

		//--load in initial pdb with a random starting structure on active chains and random residues
		nobetter = 0;
		for (UInt i = 0; i < activeChains.size(); i++)
		{
			for (UInt j = 0; j < randomResidues.size(); j++)
			{
				randomPosition.push_back(activeChains[i]);
				randomPosition.push_back(randomResidues[j]);
				mutant = getProbabilisticMutation(prot, sequencePool, possibleMutants, randomPosition, randomResidues);
				backboneAngles = prot->getRandPhiPsifromBackboneSequenceType(mutant);
				prot->setDihedral(activeChains[i], randomResidues[j],backboneAngles[0],0,0);
				prot->setDihedral(activeChains[i], randomResidues[j],backboneAngles[1],1,0);
				randomPosition.clear();
			}
		}
		prot->protMin(true);
		pdbWriter(prot, tempModel);
		Energy = prot->protEnergy();
		pastEnergy = Energy;

		//--Run through a single folding path till hitting plateau
		do
		{
			//--Try new confirmation
			nobetter++; revert = true;
			startNumBackboneClashes = prot->getNumHardBackboneClashes(); mutantPosition.clear();
			mutantPosition = getMutationPosition(activeChains, activeResidues);
			mutant = getProbabilisticMutation(prot, sequencePool, possibleMutants, mutantPosition, activeResidues);
			backboneAngles = prot->getRandPhiPsifromBackboneSequenceType(mutant);
			prot->setDihedral(mutantPosition[0],mutantPosition[1],backboneAngles[0],0,0);
			prot->setDihedral(mutantPosition[0],mutantPosition[1],backboneAngles[1],1,0);
			
			//--Prior to full optimization perform a clash test for an acceptable confirmation
			numClashes = prot->getNumHardBackboneClashes();
			if (numClashes <= startNumBackboneClashes){
				revert = false;
			}
			
			//--Energy test
			if(!revert){
				prot->protMin(true);
				Energy = prot->protEnergy();
				deltaEnergy = Energy-pastEnergy;
				boltzmannAcceptance = prot->boltzmannEnergyCriteria(deltaEnergy,::KT);
				if (boltzmannAcceptance){
					pdbWriter(prot, tempModel);
					pastEnergy = Energy; nobetter = 0;
				}
				else{revert = true;}
			}
			
			//--Revert to previously saved structure
			if (revert){
				delete thePDB;
				thePDB = new PDBInterface(tempModel);
				theEnsemble = thePDB->getEnsemblePointer();
				pMol = theEnsemble->getMoleculePointer(0);
				prot = static_cast<protein*>(pMol);
			}
		}while (nobetter < plateau);
		delete thePDB;

		//--Print final energy and write a pdb file----------------------------------------------------
		PDBInterface* theModelPDB = new PDBInterface(tempModel);
		ensemble* theModelEnsemble = theModelPDB->getEnsemblePointer();
		molecule* modelMol = theModelEnsemble->getMoleculePointer(0);
		protein* model = static_cast<protein*>(modelMol);
		
		//-Determine probability of being accepted into pool
		Energy = model->protEnergy();
		deltaEnergy = Energy-startEnergy;
		boltzmannAcceptance = model->boltzmannEnergyCriteria(deltaEnergy,::KT);
		startEnergy = Energy;
		
		//-generate pdb output
		sec = time(NULL); timeid = sec;
		stringstream convert; string countstr;
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
				if (boltzmannAcceptance || count < ::populationBaseline){
					fs << finalSequence[i][j] << ",";
				}
			}
		}
		if (boltzmannAcceptance || count < ::populationBaseline){
			finalline << " pool";
			fs << endl;
		}
		fs.close();
		finalline << endl;
		finalline.close();
		
		//-clear variables for next iteration
		delete theModelPDB;
		sequencePool.clear(), mutantPosition.clear(), chainSequence.clear(), randomPosition.clear();
		sequencePool.resize(0),  mutantPosition.resize(0), chainSequence.resize(0), randomPosition.resize(0);
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

UInt getProbabilisticMutation(protein* _prot, vector < vector < UInt > > &_sequencePool, vector < vector < UInt > > &_possibleMutants, UIntVec &_mutantPosition, UIntVec &_activeResidues)
{
	bool acceptance = false, Cterm = false;
	double prob0, prob1, prob2, Pi, Pj, poolSize = _sequencePool.size();
	vector <UInt> Freqs(20,1);
	UInt mutant, type0, type1, type2;
	UInt count = getSizeofPopulation();
	
	// determine whether probability will consist of two postitions at terminus or three positions
	if (_mutantPosition[1] == _prot->getNumResidues(_mutantPosition[0])-1){Cterm = true;}
	UInt positionPossibles = _possibleMutants[(_mutantPosition[0]+1)*_mutantPosition[1]].size();

	//--determine boltzmann probability based chance of conformation acceptance
	do
	{
		mutant = _possibleMutants[(_mutantPosition[0]+1)*_mutantPosition[1]][rand() % positionPossibles];
		if (count >= ::populationBaseline){
			//--get population conformations for locality of position
			for (UInt i = 0; i < poolSize; i++)
			{
				if (!Cterm){
					type1 = _sequencePool[i][(_mutantPosition[0]+1)*_mutantPosition[1]];
					Freqs[type1] += 1;
					type2 = _sequencePool[i][((_mutantPosition[0]+1)*_mutantPosition[1])+1];
					Freqs[type2] += 1;
				}
				else{
					type0 = _sequencePool[i][((_mutantPosition[0]+1)*_mutantPosition[1])-1];
					Freqs[type0] += 1;
					type1 = _sequencePool[i][(_mutantPosition[0]+1)*_mutantPosition[1]];
					Freqs[type1] += 1;
				}
			}
			if (!Cterm){
				prob1 = Freqs[mutant]/(poolSize-1);
				prob2 = Freqs[_prot->getBackboneSequenceType(_mutantPosition[0],_mutantPosition[1]+1)]/(poolSize-1);
				Pi = prob1*prob2;
				prob1 = 1/_possibleMutants[(_mutantPosition[0]+1)*_mutantPosition[1]].size(), prob2 = 1/_possibleMutants[((_mutantPosition[0]+1)*_mutantPosition[1])+1].size();
				Pj = prob1*prob2;
			}
			else{
				prob0 = Freqs[_prot->getBackboneSequenceType(_mutantPosition[0],_mutantPosition[1]-1)]/(poolSize-1);
				prob1 = Freqs[mutant]/(poolSize-1);
				Pi = prob0*prob1;
				prob1 = 1/_possibleMutants[(_mutantPosition[0]+1)*_mutantPosition[1]].size(), prob0 = 1/_possibleMutants[((_mutantPosition[0]+1)*_mutantPosition[1])-1].size();
				Pj = prob0*prob1;
			}
			acceptance = _prot->boltzmannProbabilityCriteria(Pi, Pj, ::KT);
		}
		else{acceptance = true;}  //random conformation
	}while (!acceptance);
	return mutant;
}

vector < vector < UInt > > readSequencePool()
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

vector < vector < UInt > > readPossibleMutants()
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
	fstream pm; bool active;
	pm.open ("possiblemutants.out", fstream::in | fstream::out | fstream::app);
	
	for (UInt i = 0; i < _activeChains.size(); i++)
	{
		for (UInt j = 0; j <_prot->getNumResidues(i); j++)
		{
			active = false;
			for (UInt k = 0; k < _activeResidues.size(); k++)
			{
				if(j == _activeResidues[k]){active = true;}
			}
			if (active){
				for (UInt l = 0; l <_allowedTypes.size(); l++)
				{
					if (_allowedTypes[l] < 11 || (_allowedTypes[l] > 10 && _prot->getTypeFromResNum(_activeChains[i],j) > 25 )){
						pm << _allowedTypes[l] << ",";
					}
				}
			}
			else{pm << _prot->getBackboneSequenceType(_activeChains[i],j) << ",";}
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

