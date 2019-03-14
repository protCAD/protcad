#include <iostream>
#include <string>
#include <vector>
#include "ensemble.h"
#include "PDBInterface.h"
#include "ran.h"
#include "generalio.h"


void optimizeRotations(protein* _prot, vector <UIntVec> _activePositions, UInt _numSteps, double _stepSize, ran &_ranNumber);
void mapSequenceOntoProtein (protein* _prot, vector <UInt> _tempSequence);
vector <string> parse(string& _currentLine);


int main (int argc, char* argv[])
{


	enum aminoAcid {A, R, N, D, C, Q, E, G, H, I, L, K, M, F, P, S, T, W, Y, V};


	//read command line
	if (argc < 2)
	{
		cout << "glycoDocker glycoinputfile.inp" << endl;
		exit(1);
	}

	//read inputfile

	string inputFileName = argv[1];
	ifstream inFile;
	inFile.open(inputFileName.c_str());
	if (!inFile)
	{
		cout << "Unable to find or open file" << endl;
		exit(1);
	}
	string currentLine;
	vector <string> parsedStrings;
	parsedStrings.resize(0);


	// read name of file from inputfile
	getline(inFile, currentLine, '\n'); parsedStrings = parse(currentLine);
	 string pdb1 = parsedStrings[0];
	 

	//read random seed from inputfile
	ran ranNumber;
	int randomSeed;
	getline(inFile, currentLine, '\n');
 	 sscanf(currentLine.c_str(), "%d", &randomSeed);

	ranNumber.setSeed((UInt)randomSeed);


	ofstream fOutput ("output.log");

	// read in first prot structure

	PDBInterface* thePDB = new PDBInterface(pdb1);
	ensemble* theEnsemble = thePDB->getEnsemblePointer();
	molecule* theMol = theEnsemble->getMoleculePointer(0);
	protein* prot = static_cast<protein*>(theMol);

	prot-> silenceMessages();
	prot->symmetryLinkChainAtoB(1,0);//must symmetry link chains first!!
	prot-> activateAllForRepacking(0);
	prot-> setCanonicalHelixRotamersOnly(0);
	prot->mutate(0,0,9);
	pdbWriter(prot,"mutate.pdb");
	prot->setRotamer(0,0,0,7);
	pdbWriter(prot, "rotate7.pdb");
	prot->setRotamer(0,0,0,8);
	pdbWriter(prot,"rotate8.pdb");

	residue::setCutoffDistance(12.0);
	pmf::setScaleFactor(0.0);
	rotamer::setScaleFactor(1.0);
	microEnvironment::setScaleFactor(0.0);
	amberVDW::setScaleFactor(1.0);
	amberVDW::setRadiusScaleFactor(1.0);
	amberVDW::setLinearRepulsionDampeningOn();
	amberElec::setScaleFactor(0.0);
	amberElec::setDielectricConstant(10.0);
	amberElec::distanceDependance = true;
	solvation::setItsScaleFactor(0.0);
	//double hydrogenBondScaleFactor = 0.0;


	//...for the optimize rotations function

	//makes matrix of residues

	//UInt numchain = 0;//assumes symmetry of chains 0 and 1

	/*vector <vector <UInt> > activePositions;

	activePositions.resize(0);
	for (UInt numres = 0; numres < 22; numres++)// outside of previous function... assumes all sequences same size as that of wild type sequence.
	{
		vector <UInt> position;
		position.resize(0);
		position.push_back(0);
		position.push_back(numres);
		activePositions.push_back(position);
	}


	prot->optimizeRotamers();
	optimizeRotations(prot, activePositions, 2, 15, ranNumber);


	//calculate energy for nmr structure
	double dimerEnergy = prot->intraEnergy();//no H bond energy factor
	double refEnergy = prot->intraEnergy(0);
	double energy = dimerEnergy - 2*refEnergy;

	fOutput << "NMR structure: energy of dimer:" << dimerEnergy << endl;
	fOutput << "NMR structure: energy of monomer:" << refEnergy << endl;
	fOutput << "NMR structure: association energy:" << energy << endl;
	*/
	prot->optimizeRotamers();

	pdbWriter(prot, "out.pdb");



	fOutput.close();
	return 0;
}

void optimizeRotations(protein* _prot, vector <vector <UInt> > _activePositions, UInt _numSteps, double _stepSize, ran &_ranNumber)
{
	for (UInt n = 0; n < 3*_activePositions.size(); n++)
	{
		UInt i = (UInt)(_ranNumber.getNext()*_activePositions.size());
		UInt numChis = _prot->getNumChis(_activePositions[i][0], _activePositions[i][1], 0);
		if (numChis > 1)
		{
			UIntVec rotamer = _prot->getCurrentRotamer(_activePositions[i][0], _activePositions[i][1]);
			_prot->setRotamer(_activePositions[i][0], _activePositions[i][1], 0, rotamer[0]);
			_prot->optimizeSmallRotations(_activePositions[i], _numSteps, _stepSize);
		}
	}
	return;
}

void mapSequenceOntoProtein (protein* _prot, vector <UInt> _tempSequence)
{
	for (UInt aaNum = 0; aaNum < _tempSequence.size(); aaNum++)//go thru all aa's in each sequence
	{
		_prot->mutate(0,aaNum,_tempSequence[aaNum]);//by mutating one chain the other is mutated as well
	}
	return;
}

vector <string> parse(string& _currentLine)
{
	StrVec parsedStrings;
	string tmpStrChi;
	string tmpStr;
	parsedStrings.resize(0);
	tmpStrChi.resize(1);
	tmpStr.resize(0);

	for (UInt i = 0; i < _currentLine.size(); i++)
	{
		if (_currentLine[i] != ' ' && _currentLine[i] != '\t' && _currentLine[i] != '\n' && i != (_currentLine.size()-1))
		{
			tmpStrChi[0] = _currentLine[i];
			tmpStr.append(tmpStrChi);
		}
		else if (tmpStr.size() != 0)
		{
			parsedStrings.push_back(tmpStr);
			tmpStr.resize(0);
		}
	}

	for (UInt i = 0; i < parsedStrings.size(); i++)
	{
		cout << parsedStrings[i] << " ";
	}
	cout << endl;
	return parsedStrings;
}
