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
UInt convertSeqStringtoInt(string Seq, string backboneSeq[], UInt size);
vector < vector < UInt > > readSequencePool();
vector < vector < UInt > > readPossibleMutants();

enum structure {C,L,P,T,E,Y,A,I,D,Q,R,F,H,W,K,S};
string backboneSeq[] =   { "C", "L", "P", "T","E","Y","A","I",  "D",  "Q",  "R",  "F", "H", "W", "K", "S"};
string backboneTypes[] = {"-π","-α","-ρ","-β","β","ρ","α","π","-πi","-αi","-ρi","-βi","βi","ρi","αi","πi"};
UInt seqSize = sizeof(backboneSeq)/sizeof(backboneSeq[0]);
UInt populationBaseline = 1000;

//--Program setup----------------------------------------------------------------------------------------
int main (int argc, char* argv[])
{
	//--Running parameters
	if (argc !=1)
	{
		cout << "protFolder" << endl;
		exit(1);
	}
	
	UIntVec activeChains, allowedTypes, activeResidues, randomResidues;
	string infile;
	ifstream file("folder.in");
	if (file){
		string item, line;
		UInt delimitercounter, linecounter = 0;
		while(getline(file,line))
		{
			delimitercounter = 0;
			stringstream stream(line);
			while(getline(stream,item,','))
			{
				if (delimitercounter > 0){
					if (linecounter == 0){
						infile = item;
					}
					if (linecounter == 1){
						stringstream inputString(item);
						UInt index;
						inputString >> index;
						activeChains.push_back(index);
					}
					if (linecounter == 2){
						stringstream inputString(item);
						UInt index;
						inputString >> index;
						activeResidues.push_back(index);
					}
					if (linecounter == 3){
						stringstream inputString(item);
						UInt index;
						inputString >> index;
						randomResidues.push_back(index);
					}
					if (linecounter == 4){
						UInt index = convertSeqStringtoInt(item, backboneTypes, seqSize);
						allowedTypes.push_back(index);
					}
				}
				delimitercounter++;
			}
			linecounter++;
		}
		file.close();
	}
	else{
		fstream inf;
		inf.open ("folder.in", fstream::in | fstream::out | fstream::app);
		inf << "Input PDB File,xyz.pdb," << endl;
		inf << "Active Chains,0,1,2," << endl;
		inf << "Active Positions,0,1,2,3,5,6,7,9,10," << endl;
		inf << "Random Positions,0,2,5,6,10," << endl;
		inf << "Backbone Types,C,L,P,T,E,Y,A,I,D,Q,R,F,H,W,K,S" << endl;
		cout << "Error: Required input file doesn't exist." << endl << "Template input file has been generated, please fill it out and rerun." << endl;
		exit(1);
	}
	
	//--running parameters
	residue::setElectroSolvationScaleFactor(1.0);
	residue::setHydroSolvationScaleFactor(1.0);
	residue::setPolarizableElec(true);
	amberElec::setScaleFactor(1.0);
	amberVDW::setScaleFactor(1.0);
	residue::setTemperature(300);

	//--set initial variables
	int seed = (int)getpid()*(int)gethostid(); srand (seed);
	double startEnergy = 1E10, pastEnergy, Energy, deltaEnergy, RMSD;
	vector <double> backboneAngles(2);
	UInt timeid, sec, numClashes, startNumBackboneClashes;
	UInt mutant = 0, plateau = 100, nobetter = 0;
	vector < UInt > mutantPosition, chainSequence, randomPosition;
	vector < vector < UInt > > sequencePool, finalSequence, possibleMutants;
	stringstream convert; string startstr, outFile;
	UInt name = rand() % 100000000;
	convert << name, startstr = convert.str();
	string tempModel = startstr + "_temp.pdb";
	bool revert, acceptance;

	//-build possible fold sequence database per position
	PDBInterface* thePDBi = new PDBInterface(infile);
	ensemble* theEnsemblei = thePDBi->getEnsemblePointer();
	molecule* pMoli = theEnsemblei->getMoleculePointer(0);
	protein* proti = static_cast<protein*>(pMoli);
	possibleMutants = readPossibleMutants();
	if(possibleMutants.size() < activeResidues.size())
	{
		createPossibleMutantsDatabase(proti, activeChains, activeResidues, allowedTypes);
		possibleMutants = readPossibleMutants();
	}
	
	//--Run multiple independent fold cycles-----------------------------------------------------
	while(true)
	{
		PDBInterface* thePDB = new PDBInterface(infile);
		ensemble* theEnsemble = thePDB->getEnsemblePointer();
		molecule* pMol = theEnsemble->getMoleculePointer(0);
		protein* prot = static_cast<protein*>(pMol);
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
			startNumBackboneClashes = prot->getNumHardBackboneClashes();
			mutantPosition.clear(); mutantPosition = getMutationPosition(activeChains, activeResidues);
			mutant = getProbabilisticMutation(sequencePool, possibleMutants, mutantPosition);
			backboneAngles = prot->getRandPhiPsifromBackboneSequenceType(mutant);
			prot->setDihedral(mutantPosition[0],mutantPosition[1],backboneAngles[0],0,0);
			prot->setDihedral(mutantPosition[0],mutantPosition[1],backboneAngles[1],1,0);
			
			//--backbone clash check
			numClashes = prot->getNumHardBackboneClashes();
			if (numClashes <= startNumBackboneClashes){
				revert = false;
			}
			
			//--Sidechain relaxation and energy calculation
			if(!revert){
				prot->protRelax(1000); nobetter++;
				Energy = prot->protEnergy();
				deltaEnergy = Energy-pastEnergy;
				acceptance = prot->boltzmannEnergyCriteria(deltaEnergy);
				
				//--Boltzmann energy acceptance of structure and save
				if (acceptance){
					prot->saveCurrentState();
					pdbWriter(prot, tempModel);
					pastEnergy = Energy;
					if (deltaEnergy < (residue::getKT()*-1)){nobetter = 0;}
				}
			}
			//--Revert to previous structure state
			if (revert){prot->undoState();}
		}while (nobetter < plateau);
		delete thePDB;

		//--Print final energy and write a pdb file----------------------------------------------------
		PDBInterface* theModelPDB = new PDBInterface(tempModel);
		ensemble* theModelEnsemble = theModelPDB->getEnsemblePointer();
		molecule* modelMol = theModelEnsemble->getMoleculePointer(0);
		protein* model = static_cast<protein*>(modelMol);
		
		//-Determine probability of being accepted into pool
		model->protMin(true);
		RMSD = proti->getRMSD(model);
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
		finalSequence.clear(), chainSequence.clear();
		for (UInt i = 0; i < activeChains.size(); i++)
		{
			chainSequence = getChainSequence(model, activeChains[i]);
			finalSequence.push_back(chainSequence);
		}
		fstream finalline; finalline.open ("results.out", fstream::in | fstream::out | fstream::app);
		finalline << timeid << " " << RMSD << " " << Energy << " ";
		fstream fs; fs.open ("sequencepool.out", fstream::in | fstream::out | fstream::app);
		for (UInt i = 0; i < activeChains.size(); i++)
		{
			for (UInt j = 0; j < finalSequence[i].size(); j++)
			{
				finalline << backboneSeq[finalSequence[i][j]] << " ";
				if (acceptance){
					fs << finalSequence[i][j] << ",";
				}
			}
		}
		if (acceptance){
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
	UInt mutant, type;
	UInt positionPossibles = _possibleMutants[(_mutantPosition[0]+1)*_mutantPosition[1]].size();

	//--determine frequency based chance of mutation acceptance (statistical potential)
	do
	{
		entropy = rand() % 100+1;
		mutant = _possibleMutants[(_mutantPosition[0]+1)*_mutantPosition[1]][rand() % positionPossibles];
		if (poolSize >= ::populationBaseline){
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

UInt convertSeqStringtoInt(string Seq, string backboneSeq[], UInt size)
{
	UInt index;
	for (UInt i = 0; i < size; i++)
	{
		if (Seq.compare(backboneSeq[i]) == 0){index = i;}
	}
	return index;
}
