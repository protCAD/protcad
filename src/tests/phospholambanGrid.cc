#include <iostream>
#include <string>
#include <vector>
#include "ensemble.h"
#include "PDBInterface.h"
#include "ran.h"
#include "generalio.h"

vector <string> parse(string& _currentLine);
void createBackboneGrid(protein* _prot, double _zrotGrid, double _rGrid, double _yrotGrid, double __offsetGrid);
void undoBackboneGrid(protein* _prot, double _zrotGrid, double _rGrid, double _yrotGrid, double __offsetGrid);
double calcEnsembleProb (protein* _prot, vector <vector <UInt> > _seqArray, vector <string> _typeSeq, vector <double> _oldEnergyArray, vector <double> &_newEnergyArray, double _temp, ran &_ranNumber, vector <vector <UInt> > _activePositions, double &_Fitness, double &_EnergyGrid, vector <double> _monEnergyArray);
void mapSequenceOntoProtein (protein* _prot, vector <vector <UInt> > _seqArray, UInt _m);
void optimizeRotations(protein* _prot, vector <UIntVec> _activePositions, UInt _numSteps, double _stepSize, ran &_ranNumber);
void insertionSort (vector <vector<double> > &_bestEnergy, int _rows);

int main (int argc, char* argv[])
{

	enum aminoAcid {A, R, N, D, C, Q, E, G, H, I, L, K, M, F, P, S, T, W, Y, V};
	if (argc < 2)
	{
		cout << "glycoDocker glycoinputfile.inp" << endl;
		exit(1);
	}
	string inputFileName = argv[1];
	ifstream inFile;
	inFile.open(inputFileName.c_str());
	if (!inFile)
	{
		cout << "Unable to find or open file" << endl;
		exit(1);
	}
	ran ranNumber;
	string currentLine;
	vector <string> parsedStrings;
	parsedStrings.resize(0);

	// read in file names
	getline(inFile, currentLine, '\n'); parsedStrings = parse(currentLine);
	 string inputFile = parsedStrings[0];
	    //string outputFile = parsedStrings[1];   note that I am commenting out outputFile because it is not used for this application.
	 string bestFileOne = parsedStrings[1];
	 string bestFileTwo = parsedStrings[2];
	 string bestFileThree = parsedStrings[3];
	 string bestFileFour = parsedStrings[4];
	 string bestFileFive = parsedStrings[5];

	// have rasmol active?
	getline(inFile, currentLine, '\n'); parsedStrings = parse(currentLine);
	bool rasmol = false;
	if (parsedStrings[0] == "yes_rasmol") rasmol = true;
	if (rasmol) system ("rasmol&");

	// read in backbone limits *for energy min grid* : lowest parameter value, highest parameter value, step size/resolution.
	double zrotminGrid, zrotmaxGrid, rminGrid, rmaxGrid, yrotminGrid, yrotmaxGrid, offsetminGrid, offsetmaxGrid;
	double rdiffGrid,zrotdiffGrid,yrotdiffGrid,offsetdiffGrid;
	getline(inFile, currentLine, '\n'); parsedStrings = parse(currentLine);
	 sscanf(parsedStrings[0].c_str(), "%lf", &rminGrid);  sscanf(parsedStrings[1].c_str(), "%lf", &rmaxGrid);
	 sscanf(parsedStrings[2].c_str(), "%lf", &rdiffGrid);
	getline(inFile, currentLine, '\n'); parsedStrings = parse(currentLine);
	 sscanf(parsedStrings[0].c_str(), "%lf", &zrotminGrid);  sscanf(parsedStrings[1].c_str(), "%lf", &zrotmaxGrid);
	 sscanf(parsedStrings[2].c_str(), "%lf", &zrotdiffGrid);
	getline(inFile, currentLine, '\n'); parsedStrings = parse(currentLine);
	 sscanf(parsedStrings[0].c_str(), "%lf", &yrotminGrid); sscanf(parsedStrings[1].c_str(), "%lf", &yrotmaxGrid);
	 sscanf(parsedStrings[2].c_str(), "%lf", &yrotdiffGrid);
	getline(inFile, currentLine, '\n'); parsedStrings = parse(currentLine);
	 sscanf(parsedStrings[0].c_str(), "%lf", &offsetminGrid);  sscanf(parsedStrings[1].c_str(), "%lf", &offsetmaxGrid);
	 sscanf(parsedStrings[2].c_str(), "%lf", &offsetdiffGrid);


	//set initial backbone

	double zrotGrid, rGrid, yrotGrid, offsetGrid;
	   zrotGrid = zrotminGrid;
	   rGrid = rminGrid;
	   yrotGrid = yrotminGrid;
	   offsetGrid = offsetminGrid;

	cout << "zrotGrid:" << zrotGrid << endl;
	cout << "rGrid:" << rGrid << endl;
	cout << "yrotGrid:" << yrotGrid << endl;
	cout << "offsetGrid:" << offsetGrid << endl;



	// read in random seed from inputfile or the command line arguments
	int randomSeed;
	if (argc == 3)
	{
		getline(inFile, currentLine, '\n'); // read in but skip
		sscanf(argv[2], "%d", &randomSeed); // use second command line arg instead
	}
	else
	{
		getline(inFile, currentLine, '\n');
	 	 sscanf(currentLine.c_str(), "%d", &randomSeed);
	}
	cout << "Randomseed = " << randomSeed << endl;

	ranNumber.setSeed((UInt)randomSeed);

	// read in prot structure
	PDBInterface* thePDB = new PDBInterface(inputFile);
	ensemble* theEnsemble = thePDB->getEnsemblePointer();
	molecule* theMol = theEnsemble->getMoleculePointer(0);
	protein* prot = static_cast<protein*>(theMol);

	// set parameters for force field

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


	// set 2-fold symmetry

	prot->symmetryLinkChainAtoB(1,0);
	prot-> activateAllForRepacking(0);
	prot-> setCanonicalHelixRotamersOnly(0);


	//read in array of wild type, unfavorable, and favorable sequences and store them in a 2-d vector, seqArray (row-seq #, col-aa #)

	vector <vector<UInt> > seqArray;//initialize array of sequences
	seqArray.resize(0);
	vector <string> typeSeq; //initialize vector of sequence type (wild type, fav, unfav)-- determining if it lowers or increases newProb
	typeSeq.resize(0);


	while (getline(inFile, currentLine, '\n'))
	{
		parsedStrings = parse(currentLine);
		vector <UInt> tempSequence; // one peptide sequence
		tempSequence.resize(0);
		if (parsedStrings[0] == "w" || parsedStrings[0] == "f" || parsedStrings[0] == "u")
		{
			typeSeq.push_back(parsedStrings[0]);
		}
		else
		{
			cout << "ERROR - sequence designation undefined... quitting" <<endl;
			exit(1);
		}

		for (UInt j = 1; j < parsedStrings.size(); j++)
		{
			UInt resId = 1000;
			if (parsedStrings[j] == "A") resId = A;
			if (parsedStrings[j] == "R") resId = R;
			if (parsedStrings[j] == "N") resId = N;
			if (parsedStrings[j] == "D") resId = D;
			if (parsedStrings[j] == "C") resId = C;
			if (parsedStrings[j] == "Q") resId = Q;
			if (parsedStrings[j] == "E") resId = E;
			if (parsedStrings[j] == "G") resId = G;
			if (parsedStrings[j] == "H") resId = H;
			if (parsedStrings[j] == "I") resId = I;
			if (parsedStrings[j] == "L") resId = L;
			if (parsedStrings[j] == "K") resId = K;
			if (parsedStrings[j] == "M") resId = M;
			if (parsedStrings[j] == "F") resId = F;
			if (parsedStrings[j] == "P") resId = P;
			if (parsedStrings[j] == "S") resId = S;
			if (parsedStrings[j] == "T") resId = T;
			if (parsedStrings[j] == "W") resId = W;
			if (parsedStrings[j] == "Y") resId = Y;
			if (parsedStrings[j] == "V") resId = V;
			if (resId == 1000)
			{
				cout << "ERROR - bad amino acid designation ... quitting." << endl;
				exit(1);
			}
			tempSequence.push_back(resId);
		}
		seqArray.push_back(tempSequence);
	}





	// create backbone
	prot->silenceMessages();
	createBackboneGrid(prot, zrotGrid, rGrid, yrotGrid, offsetGrid);
	mapSequenceOntoProtein(prot, seqArray, 0);
	prot->optimizeRotamers();
	UIntVec rotamer=prot->getCurrentRotamer(0,20);
	cout << "TYR rotamer: " << rotamer[0] << endl;


	pdbWriter(prot, "start.pdb");
	ofstream fStart("startParamsGrid.txt");
	fStart << "zrotGrid " << zrotGrid << " r " << rGrid << " yrot " << yrotGrid << " offsetGrid " << offsetGrid << endl;
	fStart.close();
	// undo backbone
	undoBackboneGrid(prot, zrotGrid, rGrid, yrotGrid, offsetGrid);


	//declare variables needed for calcEnsembleProb function

	vector <double> oldEnergyArray;
	vector <double> energyHighestProb;//energies of sequences in highest probable configuration.
	vector <double> newEnergyArray;
	oldEnergyArray.resize(0);
	energyHighestProb.resize(0);
	newEnergyArray.resize(0);
	for (UInt seqNumber = 0; seqNumber < seqArray.size();seqNumber++)
	{

		newEnergyArray.push_back (0.0);
		energyHighestProb.push_back(0.0);
	}

	double newProb = 0.0;




	//...for the optimize rotations function

	//first "for" loop makes '2' row by 'numres' column array of residues (first column = chain #, second column = res #);
	// second "for" loop prints this matrix for confirmation.

	//UInt numchain = 0;//assumes symmetry of chains 0 and 1

	vector <vector <UInt> > activePositions;

	activePositions.resize(0);
	for (UInt numres = 0; numres < seqArray[0].size(); numres++)// outside of previous function... assumes all sequences same size as that of wild type sequence.
	{
		vector <UInt> position;
		position.resize(0);
		position.push_back(0);// where 0 = chain number 0
		position.push_back(numres);
		activePositions.push_back(position);
	}

	for (UInt i = 0; i < activePositions.size(); i++)
	{
		for (UInt j = 0; j < activePositions[i].size(); j++)
		{
			cout << " " << activePositions[i][j];
		}
		cout << endl;
	}


	// For structures (of the wild type) corresponding to a grid of values for r, yrot, zrot, and offset, find (a) the energy,
	// 	(b) the probability, and (c) the fitness. Write these values to a text file called GridEnergy.log.
	//      Store the values of the parameters for the structures with the best energy in an array.


	ofstream fBestGridLog ("GridBest.log");
	ofstream fGridLog("GridEnergy.log");
	double Fitness = 0.0;
	double EnergyGrid = 0.0;

	int counter = 0;
	vector <vector<double> > bestEnergy;
	bestEnergy.resize(0);
	vector <double> parameters;
	parameters.resize(0);

	//optimize rotamers for monomeric state, for each sequence;
	//calculate energies for this monomer with optimized rotamers
	//store those energies in a vector -> they will become energies for unfolded state.

	createBackboneGrid(prot, 0.0, 50.0, 0.0, 0.0);// radius large, helices parallel

	vector <double> monEnergyArray;
	monEnergyArray.resize(0);

	for (UInt m = 0; m < seqArray.size(); m++)
	{
		//optimizing rotamers
		mapSequenceOntoProtein (prot, seqArray, m);
		prot->optimizeRotamers();
		optimizeRotations(prot, activePositions, 2, 15, ranNumber);
		//storing energies of monomer (rotamers optimized) in a vector
		double intraEnergyMon = prot->intraEnergy(0);//no H bond energy factor
		monEnergyArray.push_back(intraEnergyMon);
	}


	ofstream fMonomer ("monomerEnergy.log");
	fMonomer << " Energies of Monomer for Each Sequence, with rotamers optimized for the monomeric state " << endl;

	for (UInt num = 0; num < monEnergyArray.size(); num++)
	{
		fMonomer << " " << monEnergyArray[num] << " ";
	}
	fMonomer.close();

	undoBackboneGrid(prot, 0.0, 50.0, 0.0, 0.0);


	cout << endl << rGrid << " " << zrotGrid << " " << yrotGrid << " " << offsetGrid << endl;
	cout << endl << rminGrid << " " << zrotminGrid << " " << yrotminGrid << " " << offsetminGrid << endl;

	for (rGrid = rminGrid; rGrid <= rmaxGrid; rGrid = rGrid + rdiffGrid)
	{
		for (zrotGrid = zrotminGrid; zrotGrid <= zrotmaxGrid; zrotGrid = zrotGrid + zrotdiffGrid)
		{
			for (yrotGrid = yrotminGrid; yrotGrid <= yrotmaxGrid; yrotGrid = yrotGrid + yrotdiffGrid)
			{
				for (offsetGrid = offsetminGrid; offsetGrid <= offsetmaxGrid; offsetGrid = offsetGrid + offsetdiffGrid)
				{



					cout << "Made it to for loop" << endl;
					cout << "offsetGrid: " << offsetGrid << "offsetdiffGrid: " << offsetdiffGrid << endl;
					createBackboneGrid(prot, zrotGrid, rGrid, yrotGrid, offsetGrid);

						newProb = calcEnsembleProb (prot, seqArray, typeSeq, oldEnergyArray, newEnergyArray, 5000.0, ranNumber, activePositions, Fitness, EnergyGrid, monEnergyArray);


						fGridLog << " " << rGrid << " " << zrotGrid << " " << yrotGrid << " " << offsetGrid << " ";
						fGridLog << " " << EnergyGrid << " " << Fitness << " " << newProb << " " << endl;


					//for first five structures, add all parameters and energies to the top 5 list.
					if (counter < 5)
					{
						parameters.push_back(rGrid);
						parameters.push_back(zrotGrid);
						parameters.push_back(yrotGrid);
						parameters.push_back(offsetGrid);
						parameters.push_back(EnergyGrid);
						parameters.push_back(Fitness);
						parameters.push_back(newProb);
						bestEnergy.push_back(parameters);
						parameters.resize(0);
						counter++;
					}





					//after first five structures... first sort from lowest to highest energies; then
					//check if energy is among 5 lowest; if energy is among 5 lowest,add this structure to bestEnergy
					//array and resort the array.
					// ... if energy not among 5 lowest, do nothing.



					// at first, after building vector, sort...
					if (counter == 5)
					{
						insertionSort(bestEnergy, 5);
					}





					if (counter > 4)
					{
						vector <double> temporary;
						temporary.resize(0);
						temporary.push_back(rGrid);
						temporary.push_back(zrotGrid);
						temporary.push_back(yrotGrid);
						temporary.push_back(offsetGrid);
						temporary.push_back(EnergyGrid);
						temporary.push_back(Fitness);
						temporary.push_back(newProb);

						if (temporary[4] < bestEnergy[4][4])//energy among 5 lowest, replace highest energy structure in top 5 with new structure, resort.
						{
							//??? remove 4th vector (row) of bestEnergy... need to learn... p. 491-Matt's book
							bestEnergy[4] = temporary;
							insertionSort(bestEnergy, 5);// resort the array
						}
					}

					if (counter > 4)
					{
						fBestGridLog << "top five energies" << endl;
						for (int a = 0; a < 5 ; a++)
						{
							for (int b = 0; b < 7; b++)
							{
								fBestGridLog << bestEnergy[a][b] << " ";
							}
							fBestGridLog << endl;
						}
						fBestGridLog << endl;
					}



					undoBackboneGrid(prot, zrotGrid, rGrid, yrotGrid, offsetGrid);

				}
			}
		}
	}

	fGridLog.close();
	fBestGridLog.close();

//creates a pdb file for the 5 lowest energy structures

//	int sNumber;
//	for (sNumber = 0; sNumber < 5; sNumber++)
//	{
//		rGrid = bestEnergy[sNumber][0];
//		zrotGrid = bestEnergy[sNumber][1];
//		yrotGrid = bestEnergy[sNumber][2];
//		offsetGrid = bestEnergy[sNumber][3];
//		createBackboneGrid(prot, zrotGrid, rGrid, yrotGrid, offsetGrid);
//		mapSequenceOntoProtein (prot, seqArray, 0); // assumes wild type is first!
//		if (sNumber == 0)
//		   pdbWriter(prot, bestFileOne);
//		if (sNumber == 1)
//		   pdbWriter(prot, bestFileTwo);
//		if (sNumber == 2)
//		   pdbWriter(prot, bestFileThree);
//		if (sNumber == 3)
//		   pdbWriter(prot, bestFileFour);
//		if (sNumber == 4)
//		   pdbWriter(prot, bestFileFive);
//		undoBackboneGrid(prot, zrotGrid, rGrid, yrotGrid, offsetGrid);
//	}

return 0;
}




void createBackboneGrid(protein* _prot, double _zrotGrid, double _rGrid, double _yrotGrid, double _offsetGrid)
{
	//cout << "create " << _zrotGrid << " " << _r << " " << _yrotGrid << " " << _offsetGrid << endl;
	_prot->rotate(0,Z_axis,_zrotGrid);
	_prot->rotate(1,Z_axis,_zrotGrid);
	_prot->translate(0,0,_rGrid,_offsetGrid);
	_prot->translate(1,0,_rGrid,_offsetGrid);
	_prot->rotate(0,Y_axis,_yrotGrid);
	_prot->rotate(1,Y_axis,_yrotGrid);
	_prot->rotate(1,Z_axis,180);
	return;
}


void undoBackboneGrid(protein* _prot, double _zrotGrid, double _rGrid, double _yrotGrid, double _offsetGrid)
{
	//cout << "undo " << _zrotGrid << " " << _rGrid << " " << _yrotGrid << " " << _offsetGrid << endl;
	_prot->rotate(1, Z_axis,180);
	_prot->rotate(0, Y_axis,-1*_yrotGrid);
	_prot->rotate(1, Y_axis,-1*_yrotGrid);
	_prot->translate(0, 0, -1*_rGrid, -1*_offsetGrid);
	_prot->translate(1, 0, -1*_rGrid, -1*_offsetGrid);
	_prot->rotate(1, Z_axis,-1*_zrotGrid);
	_prot->rotate(0, Z_axis,-1*_zrotGrid);
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



double calcEnsembleProb (protein* _prot, vector <vector <UInt> > _seqArray, vector <string> _typeSeq, vector <double> _oldEnergyArray, vector <double> &_newEnergyArray, double _temp, ran &_ranNumber, vector <vector <UInt> > _activePositions, double &_Fitness, double &_EnergyGrid, vector <double> _monEnergyArray)
{
	//initialize local variables

		//for in general
		double energy = 0.0;
		double e = 2.71828182846;
		//vector < double > referenceEnergy;
		_newEnergyArray.resize(0);
		double EnergyWildT = 0.0;


		// for probability function
		double newProb = 0.0;
		double num = 0.0;
		double den = 0.0;

		// for fitness function
		double PermSum = 0.0;

		double DisrFactor = 0.0;
		double EnergyDiff = 0.0;

	//calculates energy for each sequence
	for (UInt m = 0; m < _seqArray.size(); m++)
	{
		UInt wtIndex = 0;
		cout << "Sequence Number" << m << endl;


		if (_typeSeq[m] == "w")// for wild type sequence, optimize ALL sidechains
		{
			mapSequenceOntoProtein (_prot, _seqArray, m);
			_prot->optimizeRotamers();
			optimizeRotations(_prot, _activePositions, 2, 15, _ranNumber);
			_prot->saveCurrentState();
			wtIndex = m;
		}
		else  // only optimize all sidechains for the wildtype, else just optimize changed positions
		{
			_prot->undoState();
			_prot->saveCurrentState();
			vector < UIntVec > mutations(0);

			for (UInt i = 0; i < _seqArray[0].size(); i ++)
			{
				if (_seqArray[m][i] != _seqArray[wtIndex][i])//changes mutated positions for sequence m.
				{
					cout << "sequence " << m << " ... mutating " << i << " to " << _seqArray[m][i] << endl;
					_prot->mutate(_activePositions[i][0], _activePositions[i][1], _seqArray[m][i]);
					if (_seqArray[m][i] != 0 && _seqArray[m][i] != 7 && _seqArray[m][i] != 19) mutations.push_back(_activePositions[i]);
				}
			}
			if (mutations.size() != 0) // optimizes rotamer for that changed position.
			{
				_prot->optimizeRotamers(mutations);
				optimizeRotations(_prot, mutations, 2, 15, _ranNumber);
			}
		}

		//addition, to make pdb of two of the mutant sequences
		if (m == 6 || m == 7)
		{
			if (m == 6)
			{
				pdbWriter(_prot, "rotTest6.pdb");
			}
			if (m == 7)
			{
				pdbWriter(_prot, "rotTest7.pdb");
			}
		}

		if (_typeSeq[m]=="w") _prot->saveCurrentState();
		double intraEnergy = _prot->intraEnergy();//no H bond energy factor
		_newEnergyArray.push_back(intraEnergy);
		//double refEnergy = _prot->intraEnergy(0);//energy of just the chain itself
		//referenceEnergy.push_back(refEnergy);
	}

	cout << "calculated energy for sequences" << endl;


	ofstream fenergylog("energy.log", (ios::out | ios::app));//second argument so program does not overwrite itself after each iteration.
	ofstream fdimerlog ("dimer.log", (ios::out | ios::app));// list of energies of dimer

	//go thru all the sequences to calculate all the needed values for the different functions
	UInt n = 0;
	for (n = 0; n < _seqArray.size(); n++)
	{
		energy = _newEnergyArray[n] - 2*_monEnergyArray[n];// energy of chains associating minus energy of chains by themselves
		fenergylog << " " << energy;
		fdimerlog << " " << _newEnergyArray[n];

		// for (less sensitive) fitness function
		if (_typeSeq[n] == "w" || _typeSeq[n] == "f")
		{
			PermSum += energy;
		}
		if (_typeSeq[n] == "u")
		{
			EnergyDiff = (_newEnergyArray[0] - 2*_monEnergyArray[0]) - energy;//energy of WT - energy of mutant, in that order
			double DisrFactorContr = 0.0;
			DisrFactorContr = pow(e, ((-EnergyDiff)/(10*energy)));


			// let the disruptive factor contribution term be the minimum of 1 or pow(e, ((-EnergyDiff)/(10*energy))).
			if (DisrFactorContr >= 1.0)
			{
				DisrFactor += pow(e, ((-EnergyDiff)/(10*energy)));
			}
			else
			{
				DisrFactor += 1.0;
			}
		}


		// for (more sensitive) probability function (only for unfavorable sequences and wild type)

		if ( _typeSeq[n] == "w")
		{
			num = pow(e, (-1*energy / (0.00258*_temp)));
		}
		if (_typeSeq[n] == "w" || _typeSeq[n] == "u")
		{
			den += pow(e, (-1* energy / (0.00258*_temp)));
		}



	}
	fenergylog << endl;
	fenergylog.close();
	fdimerlog << endl;
	fdimerlog.close();

	if (den !=0)
	{
		newProb = num/den;

	}
	else
	{
		cout << "division by zero, quitting." << endl;
		newProb = 0.0;
	}

	EnergyWildT = _newEnergyArray[0] - 2*_monEnergyArray[0];

	//calculate (the less sensitive) fitness function for a collection of sequences

	_Fitness = ((PermSum)/n) - (1/n)*(EnergyWildT)*(DisrFactor);

	//determine the energy of the wild type sequence

	_EnergyGrid = EnergyWildT;



	return newProb;
}



void mapSequenceOntoProtein (protein* _prot, vector <vector <UInt> > _seqArray, UInt _m)
{
	for (UInt aaNum = 0; aaNum < _seqArray[_m].size(); aaNum++)//go thru all aa's in each sequence
	{
		_prot->mutate(0,aaNum,_seqArray [_m][aaNum]);//by mutating one chain the other symmetry linked one is mutated as well
	}
	return;
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

//this function orders a group of strings from string with smallest 4th element to string with largest 4th element.

void insertionSort (vector <vector<double> > &_bestEnergy, int _rows)
{
	for (int i = 1; i < _rows; i++)
	{
		if (_bestEnergy[i][4] < _bestEnergy [i-1][4])
		{
			vector <double> temp = _bestEnergy[i];
			int j = i;

			do
			{
				_bestEnergy [j] = _bestEnergy [j-1];
				--j;
			} while ((j>0) && (_bestEnergy[j-1][4] > temp[4]));

		_bestEnergy[j] = temp;
		}
	}
	int m;
	cout << "Ordered List of Numbers (smallest to largest)" << endl;
			for (m=0; m < 5; m++)
			{
				cout << _bestEnergy[m][4] << " ";
			}
			cout << endl;
}

