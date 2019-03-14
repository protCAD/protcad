#include <iostream>
#include <string>
#include <vector>
#include "ensemble.h"
#include "PDBInterface.h"
#include "ran.h"
#include "generalio.h"

#define PI 3.141596
#define NMAX 5000
#define CREATEBACKBONE createHomoTetramer
#define UNDOBACKBONE undoHomoTetramer
#define LARGE 1e20
#define TINY 1e-20
#define VecVecUIntVec vector < vector < vector < unsigned int >  > >
#define VecStrVec vector < vector < string > > 
#define VecUIntVec vector < vector < unsigned int > >
#define VecDblVec vector < vector < double > >
#define FIXEDCHAIN 4
enum aminoAcid {A, R, N, D, C, Q, E, G, H, I, L, K, M, F, P, S, T, W, Y, V};
protein* prot;
bool verbose = true;
// function declarations


VecDblVec getBackbone(string _backboneFile); 
VecDblVec getMCParams(string _mcParmasFile); 
VecDblVec monteCarlo(VecDblVec _mcParams, ran _ranNumber, VecDblVec _backboneParams);
void report(VecDblVec _backboneParams, string _pdbFile, int _seqIndex); 
void createHomoTetramer(VecDblVec _backboneParams); 
void undoHomoTetramer(VecDblVec _backboneParams); 
VecDblVec modifyBackbone(VecDblVec _currentBackbone, ran &number); 
double getEnergy(); 
double getEnergyWithM2Restraints(); 
double getHBondEnergy(dblVec _BA, dblVec _BC, double _eps, double _r0, double _theta0);
double hisAmdEnergy();
double expoCool(DblVec _args, UInt _iteration); 
void alignCB(UInt _chain, UInt _res);

int main (int argc, char* argv[])
{
	residue::setCutoffDistance(8.0);
	pmf::setScaleFactor(0.0);
	microEnvironment::setScaleFactor(0.0);
	amberVDW::setScaleFactor(1.0);
	amberVDW::setRadiusScaleFactor(1.0);
	amberVDW::setLinearRepulsionDampeningOff();
	amberElec::setScaleFactor(0.0);
	solvation::setItsScaleFactor(0.0);

	// read in PDB file
	string protFile = argv[1];
	PDBInterface* thePDB = new PDBInterface(protFile);
	ensemble* theEnsemble = thePDB->getEnsemblePointer();
	molecule* theMol = theEnsemble->getMoleculePointer(0);
	prot = static_cast<protein*>(theMol);

	// read in backbone params	
	string backboneFile = argv[2];
	VecDblVec backboneParams = getBackbone(backboneFile);

	// align histidine on axis

	UInt hisPos = 1000;

	for (UInt i = 0; i < prot->getNumResidues(0); i ++)
	{
		if (prot->getTypeFromResNum(0,i) == H)
		{
			hisPos = i;
		}
	}


	dblVec oneCA = prot->getCoords(0,0,"CA");
	dblVec sevCA = prot->getCoords(0,6,"CA");
	if (sevCA[2] > oneCA[2]) 
	{
		for (UInt i = 0; i < 4; i ++)
		{
			prot->rotate(i, Y_axis, 180.0);
		}
	}

	if (hisPos != 1000)
	{
		alignCB(0,hisPos);
		alignCB(1,hisPos);
		alignCB(2,hisPos);
		alignCB(3,hisPos);
	}
	
	// silence amantadine backbone atoms
	if (prot->getNumChains() > 4)
	{
		for (UInt i = 0; i < 4; i ++)
		{
			prot->makeAtomSilent(FIXEDCHAIN,0,i);
		}
	}
	
	// read in monte carlo scheme
	string mcFile = argv[3];
	VecDblVec mcParams = getMCParams(mcFile);

	UInt seed = 0;
	sscanf(argv[4], "%u", &seed);
	ran ranNumber;
	ranNumber.setSeed(seed);

	VecDblVec minBackboneParams = monteCarlo(mcParams, ranNumber, backboneParams);
	report(minBackboneParams, "final.pdb", 0);


	return 0;
}

void alignCB(UInt _chain, UInt _res)
{
	if (verbose) cout << " calling alignCB " << endl;
	dblVec hisCB = prot->getCoords(_chain, _res, "CB");
	dblVec tempX(3);
	tempX[0]=0.0;
	tempX[1]=-1.0;
	tempX[2]=0.0;
	hisCB[2]=0.0;

	double magCB = sqrt(CMath::dotProduct(hisCB,hisCB));
	double magX = sqrt(CMath::dotProduct(tempX,tempX));

	double hisCBAngle = acos(CMath::dotProduct(hisCB,tempX)/(magX*magCB));
	hisCBAngle = (180/PI)*hisCBAngle;
	if(verbose)cout << "his angle " << hisCBAngle << endl;
	prot->rotate(_chain, Z_axis, hisCBAngle);

	return;
}

void createHomoTetramer(VecDblVec _backboneParams)
{
	if (_backboneParams.size() != 7 )
	{
		cout << "ERROR - Backbone parameters not compatable with homotetramer" << endl;
		exit(1);
	}
	if (verbose) cout << "creating backbone ... " << endl;
	double radius = _backboneParams[0][0];
	double phase = _backboneParams[1][0];
	double crossY = _backboneParams[2][0];
	double crossX = _backboneParams[3][0];
	double offset = _backboneParams[4][0];
	double zmove = _backboneParams[5][0];
	double zrot = _backboneParams[6][0];
	for (UInt i = 0; i < 4; i ++)
	{	
		prot->rotate(i,Z_axis, phase);
		prot->translate(i,0.0, 0.0, offset);
		prot->rotate(i,Y_axis, crossY);
		prot->rotate(i,X_axis, crossX);
		prot->translate(i,0.0,0.0, -1*offset);
		prot->translate(i,0.0, radius, 0.0);
		prot->translate(i,0.0, 0.0, zmove);
		prot->rotate(i,Z_axis, zrot);
	}
	prot->rotate(0, Z_axis, 90.0);
	prot->rotate(1, Z_axis, 180.0);
	prot->rotate(2, Z_axis, 270.0);
	

	return;
}

void undoHomoTetramer(VecDblVec _backboneParams)
{
	if (_backboneParams.size() != 7 )
	{
		cout << "ERROR - Backbone parameters not compatible with homotetramer" << endl;
		exit(1);
	}
	if (verbose) cout << "undoing backbone ... " << endl;
	double radius = _backboneParams[0][0];
	double phase = _backboneParams[1][0];
 	double crossY = _backboneParams[2][0];
	double crossX = _backboneParams[3][0];
	double offset = _backboneParams[4][0];
	double zmove = _backboneParams[5][0];
	double zrot = _backboneParams[6][0];



	prot->rotate(0, Z_axis, -90.0);
	prot->rotate(1, Z_axis, -180.0);
	prot->rotate(2, Z_axis, -270.0);
	for (UInt i = 0; i < 4; i ++)
	{
		prot->rotate(i,Z_axis, -1*zrot);
		prot->translate(i,0.0, 0.0, -1*zmove);
    		prot->translate(i,0.0, -1*radius, 0.0);
		prot->translate(i,0.0,0.0,offset);
		prot->rotate(i,X_axis, -1*crossX);
    		prot->rotate(i,Y_axis, -1*crossY);
		prot->translate(i,0.0, 0.0, -1*offset);
    		prot->rotate(i,Z_axis, -1*phase);
	}
	return;
}

double expoCool(DblVec _args, UInt _i)
{
	if (verbose) cout << "calling expocool" << endl;
	double bi = 0;
	double b0 = _args[0];
	double bN = _args[1];
	double N = _args[2];

	bi = b0 * pow( bN/b0, (double)_i / N);
	return bi;
}

VecDblVec modifyBackbone(VecDblVec _currentBackbone, ran &ranNumber)
{
	if (verbose) cout << "calling modifyBackobne" << endl;
	UInt param = (UInt)(ranNumber.getNext()*(double)_currentBackbone.size());
	if (param > (_currentBackbone.size() -1)) param = _currentBackbone.size() -1;	

	double min = _currentBackbone[param][1];
	double max = _currentBackbone[param][2];
	double newVal = min + ranNumber.getNext()*(max - min);

	_currentBackbone[param][0] = newVal;
	return _currentBackbone;
}



double getEnergy()
{
	if (verbose) cout << "calling getEnergy" << endl;
	double energy = 0.0;
	for (UInt i = 0; i < prot->getNumChains(); i ++)
	{
		energy += prot->intraEnergy(i);
		for (UInt j = i+1; j < prot->getNumChains(); j ++)
		{
			energy += prot->intraEnergy(i,j);
		}
	}
	return energy;
}

double getEnergyWithM2Restraints()
{
	if (verbose) cout << "calling getEnergyWithM2Restraints" << endl;
	double energy = 0.0;
	energy += getEnergy();
	if (prot->getNumChains() > 4) energy += hisAmdEnergy();
	return energy;
}

double hisAmdEnergy()
{
	UInt hisNum = 1000;
	for (UInt i = 0; i < prot->getNumResidues(0); i ++)
	{
		if (prot->getTypeFromResNum(0,i) == H) hisNum = i;
	}	
	if (hisNum == 1000)
	{
		cout << "HIS not found in the sequence" << endl;
		return 0.0;
	}
	dblVec HisNE2 = prot->getCoords(0,hisNum,"NE2");
	dblVec AmdCB9 = prot->getCoords(4,0,"CB9");
	dblVec AmdNB0 = prot->getCoords(4,0,"NB0");

	dblVec BA = AmdNB0 - HisNE2;
	dblVec BC = AmdNB0 - AmdCB9;

	double angle = 120.0 * PI / 180.0;

	double HBenergy = getHBondEnergy(BA,BC,1.0,2.8,angle);
	cout << "HIS AMD HBond " << HBenergy << " kcals " << endl;
	return 4.0 * HBenergy;
}

double getHBondEnergy(dblVec _BA, dblVec _BC, double _eps, double _r0, double _theta0)
{
	double magBC = sqrt(CMath::dotProduct(_BC, _BC));
	double magBA = sqrt(CMath::dotProduct(_BA, _BA));

	double angle = acos( (CMath::dotProduct(_BA,_BC)) / (magBA * magBC) );

	double distratio = _r0/magBA;
	double angleFactor = pow(cos(angle - _theta0),2);
	double energy = angleFactor * _eps*(5*pow(distratio,12) - 6 * pow(distratio,10));

	return energy;
}



VecDblVec getBackbone(string _backboneFile)
{
	if (verbose) cout << "getBackbone reading in backbone file " << endl;
	ifstream inFile;
	inFile.open(_backboneFile.c_str());
	if (!inFile)
	{
		cout << "Unable to find backbone file " << inFile << endl;
		exit(1);
	}
	string currentLine;
	StrVec parsedStrings;
	VecDblVec backbone(0);
	while (getline(inFile, currentLine, '\n'))
	{
		parsedStrings = Parse::parse(currentLine);
		DblVec chainBackbone(0);
		for (UInt i = 0; i < parsedStrings.size(); i ++)
		{
			double param;
			sscanf(parsedStrings[i].c_str(), "%lf", &param);
			chainBackbone.push_back(param);
		}
		backbone.push_back(chainBackbone);
	}
	inFile.close();
	return backbone;
}
	
VecDblVec getMCParams(string _mcParamsFile)
{
	if (verbose) cout << "getMCParams reading in monte carlo paramter files" << endl;
	ifstream inFile;
	inFile.open(_mcParamsFile.c_str());
	if (!inFile)
	{
		cout  << "Unable to find monte carlo params file " << inFile << endl;
		exit(1);
	}
	string currentLine;
	StrVec parsedStrings;
	VecDblVec mcParams(0);
	while (getline(inFile, currentLine, '\n'))
	{
		parsedStrings = Parse::parse(currentLine);
		DblVec cycleParams(0);
		for (UInt i = 0; i < parsedStrings.size(); i ++)
		{
			double param;
			sscanf(parsedStrings[i].c_str(), "%lf", &param);
			cycleParams.push_back(param);
		}
		mcParams.push_back(cycleParams);
	}
	inFile.close();
	return mcParams;
}


void report(VecDblVec _backbone, string _pdbFile, int _index)
{
	if (verbose) cout << "calling report" << endl;
	cout << "params report" << endl;
	for (UInt i = 0; i < _backbone.size(); i ++)
	{
		cout << "backbone " << i << " " << endl;
		for (UInt j = 0; j < _backbone[i].size(); j ++)
		{
			cout << _backbone[i][j] << " ";
		}
		cout << endl;
	}
	cout << endl;
	if (_index != -1)
	{
		CREATEBACKBONE(_backbone);
		//prot->undoState();
		//prot->saveCurrentState();
		pdbWriter(prot, _pdbFile);
		UNDOBACKBONE(_backbone);
		string logfile = _pdbFile + ".out";
		ofstream fout(logfile.c_str());
		for (UInt i = 0; i < _backbone.size(); i ++)
		{
			fout << "backbone " << i << " " << endl;
			for (UInt j = 0; j < _backbone[i].size(); j ++)
			{
				fout << _backbone[i][j] << " ";
			}
			fout << endl;
		}
	}
}


VecDblVec monteCarlo(VecDblVec _mcParams, ran _ranNumber, VecDblVec _backboneParams)
{
	if (verbose) cout << "calling monteCarlo" << endl;
	double e =  2.71828182846;
	VecDblVec bestBackbone = _backboneParams;
	VecDblVec currentBackbone = _backboneParams;
	for (UInt cycle = 0; cycle < _mcParams.size(); cycle ++)
	{
		cout << "cycle " << cycle << endl;
		double bestScore = LARGE;
		double currentScore = bestScore;
		// _mParams is start, end, iterations
		DblVec coolfuncargs = _mcParams[cycle];
		double iterations = _mcParams[cycle][2];

		for ( double i = 0.0; i < iterations; i ++)
		{
			double beta = expoCool(coolfuncargs, (UInt)floor(i));
			cout << "Iteration:  " << i << " temp " << beta << endl;
			VecDblVec tmpBackbone = modifyBackbone(currentBackbone,_ranNumber);
			CREATEBACKBONE(tmpBackbone);
			double score;
			score = getEnergyWithM2Restraints(); 
			UNDOBACKBONE(tmpBackbone);
			if (score < bestScore)
			{
				bestScore = score;
				bestBackbone = tmpBackbone;
				cout << "best" << endl;
				report(bestBackbone, "best.pdb", 0);
			}
			if (score < currentScore)
			{
				currentScore = score;
				currentBackbone = tmpBackbone;
				cout << "ACCEPTED GOOD " << score <<  endl;
			}
			else
			{	
				double prob = pow(e, -1.0 * beta * (score - currentScore));
				if (_ranNumber.getNext() < prob)
				{
					currentScore = score;
					currentBackbone = tmpBackbone;
					cout << "ACCEPTED METROPOLIS " << score << endl;
				}
				else
				{
					cout << "REJECTED " << score << endl;
				}
			}
			cout << "current " << currentScore << endl;
			report(currentBackbone, "current.pdb", 0);
		}
	}
	return bestBackbone;
}	

