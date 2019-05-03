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
UInt getProbabilisticMutation(vector < vector < UInt > > &_sequencePool, vector < vector < UInt > > &_possibleMutants, UIntVec &_mutantPosition);
void createPossibleMutantsDatabase(protein* &_prot, UIntVec &_activeChains, UIntVec &_activeResidues, UIntVec &_allowedTypes);
UInt getSizeofPopulation();
vector < vector < UInt > > readSequencePool();
vector < vector < UInt > > readPossibleMutants();

enum structure {C,L,P,T,E,Y,A,I,D,Q,R,F,H,W,K,S};
string backboneSeq[] =   { "C", "L", "P", "T","E","Y","A","I",  "D",  "Q",  "R",  "F", "H", "W", "K", "S"};
string backboneTypes[] = {"-π","-α","-ρ","-β","β","ρ","α","π","-πi","-αi","-ρi","-βi","βi","ρi","αi","πi"};
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

	//--input
	UInt _activeChains[] = {0};                                                         // chains active for mutation
	UInt _allowedTypes[] = {C,L,P,T,E,Y,A,I,D,Q,R,F,H,W,K,S};                     // backbone types allowable
	UInt _activeResidues[] = {0,1,2,3,4,5,6,7,8,9,10,11};                                     // positions active for mutation
	UInt _randomResidues[] = {0,1,2,3,4,5,6,7,8,9,10,11};                                     // positions active for a random start sequence initially

	//--running parameters
	residue::setElectroSolvationScaleFactor(1.0);
	residue::setHydroSolvationScaleFactor(1.0);
	residue::setPolarizableElec(true);
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
	int seed = (int)getpid()*(int)gethostid(); srand (seed);
	double startEnergy = 1E10, pastEnergy, Energy, deltaEnergy;
	vector <double> backboneAngles(2);
	UInt timeid, sec, numClashes, startNumBackboneClashes;
	UInt mutant = 0, plateau = 50, nobetter = 0;
	vector < UInt > mutantPosition, chainSequence, randomPosition;
	vector < vector < UInt > > sequencePool, finalSequence, possibleMutants;
	stringstream convert; string infile = argv[1], startstr, outFile;
	UInt count, name = rand() % 100000000;
	convert << name, startstr = convert.str();
	string tempModel = startstr + "_temp.pdb";
	bool revert, acceptance;

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
				mutant = getProbabilisticMutation(sequencePool, possibleMutants, randomPosition);
				backboneAngles = prot->getRandPhiPsifromBackboneSequenceType(mutant);
				prot->setDihedral(activeChains[i], randomResidues[j],backboneAngles[0],0,0);
				prot->setDihedral(activeChains[i], randomResidues[j],backboneAngles[1],1,0);
				randomPosition.clear();
			}
		}
		prot->protRelax(1000);
		pdbWriter(prot, tempModel);
		Energy = prot->protEnergy();
		pastEnergy = Energy;

		//--Run through a single folding path till hitting plateau
		do
		{
			//--Try new confirmation
			revert = true;
			startNumBackboneClashes = prot->getNumHardBackboneClashes(); mutantPosition.clear();
			mutantPosition = getMutationPosition(activeChains, activeResidues);
			mutant = getProbabilisticMutation(sequencePool, possibleMutants, mutantPosition);
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
				nobetter++;
				prot->protRelax(1000);
				Energy = prot->protEnergy();
				deltaEnergy = Energy-pastEnergy;
				acceptance = prot->boltzmannEnergyCriteria(deltaEnergy);
				if (acceptance){
					pdbWriter(prot, tempModel);
					pastEnergy = Energy;
					if (deltaEnergy < (residue::getKT()*-1)){nobetter = 0;}
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
		model->protMin(true);
		Energy = model->protEnergy();
		deltaEnergy = Energy-startEnergy;
		acceptance = prot->boltzmannEnergyCriteria(deltaEnergy);
		startEnergy = Energy;
		
		//-generate pdb output
		sec = time(NULL); timeid = sec;
		stringstream convert; string countstr;
		convert << timeid, countstr = convert.str();
		outFile = countstr + "." + startstr + ".fold.pdb";
		pdbWriter(model, outFile);
		
		//-write to data files
		count = getSizeofPopulation(); finalSequence.clear(), chainSequence.clear();
		for (UInt i = 0; i < activeChains.size(); i++)
		{
			chainSequence = getChainSequence(model, activeChains[i]);
			finalSequence.push_back(chainSequence);
		}
		fstream finalline; finalline.open ("results.out", fstream::in | fstream::out | fstream::app);
		finalline << timeid << " " << Energy << " ";
		fstream fs; fs.open ("sequencepool.out", fstream::in | fstream::out | fstream::app);
		for (UInt i = 0; i < activeChains.size(); i++)
		{
			for (UInt j = 0; j < finalSequence[i].size(); j++)
			{
				finalline << backboneSeq[finalSequence[i][j]] << " ";
				if (acceptance || count < ::populationBaseline){
					fs << finalSequence[i][j] << ",";
				}
			}
		}
		if (acceptance || count < ::populationBaseline){
			finalline << " pool"; fs << endl;
		}
		fs.close(); finalline << endl; finalline.close();
		
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

UInt getProbabilisticMutation(vector < vector < UInt > > &_sequencePool, vector < vector < UInt > > &_possibleMutants, UIntVec &_mutantPosition)
{
	double Pi, entropy, poolSize = _sequencePool.size();
	vector <UInt> Freqs(16,1);
	UInt mutant, type, count = getSizeofPopulation();
	UInt positionPossibles = _possibleMutants[(_mutantPosition[0]+1)*_mutantPosition[1]].size();

	//--determine frequency based chance of mutation acceptance (statistical potential)
	do
	{
		entropy = rand() % 100+1;
		mutant = _possibleMutants[(_mutantPosition[0]+1)*_mutantPosition[1]][rand() % positionPossibles];
		if (count >= ::populationBaseline){
			for (UInt i = 0; i < poolSize; i++)
			{
				type = _sequencePool[i][(_mutantPosition[0]+1)*_mutantPosition[1]];
				Freqs[type] += 1;
			}
			Pi = (Freqs[mutant]/(poolSize-1))*100;
		}
		else{Pi = 100;}  //random mutant
	}while (entropy > Pi);
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
					pm << _allowedTypes[l] << ",";
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

